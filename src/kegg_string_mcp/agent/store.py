"""Append-only run store. Written by the pipeline, never by the model.

Everything the tools returned is appended here *before* the model is shown it, so
the store is the ground truth a citation is checked against. The model cannot
write to it, which is what makes validation meaningful: if the store said it, a
tool returned it.

JSONL, append-only, one file per run. No updates, no deletes -- a claim about
what the model was given must not be rewritable after the fact.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass
class RunStore:
    path: Path
    run_id: str
    _citable: set[str] = field(default_factory=set)
    _records: dict[str, dict[str, Any]] = field(default_factory=dict)
    _per_target: dict[str, set[str]] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def _append(self, kind: str, payload: dict[str, Any]) -> None:
        with self.path.open("a", encoding="utf-8") as fh:
            fh.write(json.dumps({"run_id": self.run_id, "at": _now(), "kind": kind, **payload}) + "\n")

    # -- writes -------------------------------------------------------------

    def tool_result(self, tool: str, arguments: dict[str, Any], result: dict[str, Any]) -> None:
        """Record a tool call and everything it returned. This is what defines the
        citable set: an ID is citable if and only if a tool actually returned it."""
        # Track which IDs came back for WHICH gene. A global membership check
        # cannot catch a citation that is real but attached to the wrong gene,
        # which is the more plausible-looking error.
        # Alias the target. `claimed_target` is what the caller asked to annotate,
        # but the model is invited to call tools with a different spelling of the
        # same gene -- katG, Rv1908c, mtu:Rv1908c. Keying only on the literal
        # argument made every ID from such a call report as cross_target, failing
        # a correctly-cited run.
        aliases = {str(arguments.get("gene", "")).strip().upper()}
        resolved = result.get("resolved", {})
        for key in ("kegg_gene_id", "string_id", "preferred_name"):
            value = str(resolved.get(key) or "").strip()
            if not value:
                continue
            aliases.add(value.upper())
            for separator in (":", "."):
                if separator in value:
                    aliases.add(value.split(separator, 1)[1].upper())
        aliases.discard("")

        # The resolved identifiers are part of what the tool returned, even though
        # they are not "records": kegg_gene_id and string_id name the very gene
        # being annotated. Omitting them made the validator flag a model for citing
        # the subject of its own query -- a false positive on correct behaviour,
        # which is the one failure a validator must never have.
        for key in ("kegg_gene_id", "string_id"):
            resolved_id = resolved.get(key)
            if resolved_id:
                self._citable.add(resolved_id)
                self._records.setdefault(resolved_id, {
                    "record_id": resolved_id, "type": "resolved_subject",
                    "name": result.get("resolved", {}).get("preferred_name", ""),
                })
                for alias in aliases:
                    self._per_target.setdefault(alias, set()).add(resolved_id)

        for record in result.get("records", []):
            rid = record.get("record_id")
            if not rid:
                continue
            self._citable.add(rid)
            existing = self._records.get(rid)
            if existing is None:
                self._records[rid] = record
            else:
                # Merge, do not overwrite. The same PMID returned for two genes was
                # keeping only the last record -- and `mentions` with it -- so a paper
                # discussing both genes was silently dropped from the first gene's
                # corpus, which is the field a downstream pipeline filters on.
                merged = set(existing.get("detail", {}).get("mentions", []))
                merged.update(record.get("detail", {}).get("mentions", []))
                if merged:
                    existing.setdefault("detail", {})["mentions"] = sorted(merged)
            for alias in aliases:
                self._per_target.setdefault(alias, set()).add(rid)
        self._append("tool_result", {"tool": tool, "arguments": arguments, "result": result})

    def derived(self, label: str, payload: dict[str, Any]) -> None:
        """Deterministic pipeline computation (set intersections, pathway sizes).
        Recorded separately from tool results so an auditor can tell which numbers
        came from an upstream API and which this code computed."""
        self._append("derived", {"label": label, "payload": payload})

    def decision(self, payload: dict[str, Any]) -> None:
        """One turn of the agent loop: what the model chose to do, and why it stopped."""
        self._append("decision", payload)

    def output(self, payload: dict[str, Any]) -> None:
        self._append("output", payload)

    # -- reads --------------------------------------------------------------

    @property
    def citable_ids(self) -> set[str]:
        return set(self._citable)

    @property
    def per_target(self) -> dict[str, set[str]]:
        return {k: set(v) for k, v in self._per_target.items()}

    @property
    def records(self) -> dict[str, dict[str, Any]]:
        """Every record the tools returned, keyed by ID. Quote validation needs the
        stored `quotable_text`, so it must read the pipeline's copy -- never the
        model's account of what it was shown."""
        return dict(self._records)

    def record(self, record_id: str) -> dict[str, Any] | None:
        return self._records.get(record_id)

    def corpus_manifest(self) -> list[dict[str, Any]]:
        """Papers this run found, deduped, for a downstream full-text pipeline.

        `mentions` is the genes actually present in the retrieved title/abstract,
        NOT the genes queried -- PubMed matches on metadata the model never sees,
        so a paper can be returned for a gene it never names. Filtering a corpus
        on the query term would quietly admit those.

        `in_pmc` is the handle that matters for full text: a DOI resolves to the
        publisher, usually paywalled and with terms forbidding bulk retrieval,
        whereas a PMCID is the licit route. Presence in PMC is still not
        sufficient -- the Open Access subset is a subset of PMC and per-article
        licences vary within it, which the consumer must check.
        """
        papers: dict[str, dict[str, Any]] = {}
        for record in self._records.values():
            if record.get("type") != "article":
                continue
            detail = record.get("detail", {})
            pmid = record["record_id"]
            entry = papers.setdefault(pmid, {
                "pmid": pmid,
                "doi": detail.get("doi", ""),
                "pmcid": detail.get("pmcid", ""),
                "in_pmc": detail.get("in_pmc", False),
                "title": detail.get("title", ""),
                "journal": detail.get("journal", ""),
                "year": detail.get("year", ""),
                "url": record.get("url", ""),
                "mentions": [],
            })
            entry["mentions"] = sorted(set(entry["mentions"]) | set(detail.get("mentions", [])))
        return [papers[k] for k in sorted(papers)]

    def replay(self) -> Iterator[dict[str, Any]]:
        if not self.path.exists():
            return iter(())
        with self.path.open(encoding="utf-8") as fh:
            for line in fh:
                if line.strip():
                    yield json.loads(line)
