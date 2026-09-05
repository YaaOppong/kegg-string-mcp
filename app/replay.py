"""Replay a pre-computed run and re-validate it. No network, no model, no API key.

The demo does not call anything. It loads a run captured earlier from the live
pipeline and re-runs `validate()` over it, so what a visitor sees is the current
validator's verdict on a real run rather than a stored screenshot of one.

That matters: the validator has been corrected repeatedly, several times to remove
false positives. Replaying a stored verdict would show conclusions the code no
longer reaches, and a demo that disagrees with the library it demonstrates is
worse than no demo.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from kegg_string_mcp.agent.store import RunStore
from kegg_string_mcp.agent.validate import ValidationReport, validate

RUNS_DIR = Path(__file__).resolve().parent.parent / "demo" / "runs"

# Which runs to offer, and in what order. Data, not presentation: it lives here so
# the ordering invariant -- that the first run offered is one where the checking
# catches something -- can be tested without installing a web framework.
#
# furA is first on purpose. A demo where everything passes proves nothing, so a
# visitor who changes nothing still sees the point.
ORDERED = ["furA", "gyrB", "katG", "gyrA", "ahpC", "rpoB", "pncA", "phoP",
           "katG-ahpC-epistasis", "phoP-phoR-epistasis"]

LABELS = {
    "furA": "furA — a transcriptional regulator KEGG has no pathway for",
    "gyrB": "gyrB — DNA gyrase subunit B, annotated mostly from literature",
    "katG": "katG — catalase-peroxidase; resistance-associated for isoniazid",
    "gyrA": "gyrA — DNA gyrase subunit A, absent from KEGG pathways",
    "ahpC": "ahpC — alkyl hydroperoxide reductase, no KEGG pathway",
    "rpoB": "rpoB — RNA polymerase beta subunit",
    "pncA": "pncA — pyrazinamidase, three KEGG pathways",
    "phoP": "phoP — carries a lineage marker; neither source describes its function",
    "katG-ahpC-epistasis": "katG + ahpC — looking for a mechanistic link between two genes",
    "phoP-phoR-epistasis": "phoP + phoR — a real interaction that may still be confounded",
}


@dataclass
class Replay:
    id: str
    mode: str
    target: str
    organism: str
    summary: str
    turns: list[dict[str, Any]]
    calls: list[dict[str, Any]]
    report: ValidationReport
    records: dict[str, dict[str, Any]]

    @property
    def failures(self) -> list[Any]:
        return self.report.unsupported + self.report.cross_target

    @property
    def clean(self) -> bool:
        return self.report.passed


def available(runs_dir: Path | None = None) -> list[str]:
    directory = runs_dir or RUNS_DIR
    return sorted(p.stem for p in directory.glob("*.json"))


def load(name: str, runs_dir: Path | None = None) -> Replay:
    directory = runs_dir or RUNS_DIR
    record = json.loads((directory / f"{name}.json").read_text(encoding="utf-8"))

    # Rebuild the store from the captured tool results, exactly as the pipeline
    # built it at run time, so the citable set is derived rather than trusted.
    store = RunStore(path=Path("/dev/null"), run_id=f"replay-{name}")
    for call in record["calls"]:
        store.tool_result(call["tool"], call["arguments"], call["result"])

    claimed = record["target"].strip().upper() if record["mode"] == "single" else None
    # store.notes as well as store.records: a summary may quote a tool note
    # verbatim, and checking only record text reported that as fabrication --
    # which the demo would have shown a visitor as a caught failure.
    report = validate(record["summary"], store.citable_ids, store.per_target,
                      claimed, records=store.records, notes=store.notes)

    return Replay(
        id=record["id"], mode=record["mode"], target=record["target"],
        organism=record.get("organism", "mtu"), summary=record["summary"],
        turns=record["turns"], calls=record["calls"], report=report,
        records=store.records,
    )


def tool_call_rows(replay: Replay) -> list[list[str]]:
    """One row per tool call, labelled with who asked for it.

    Not every call comes from the model. In epistasis mode the pipeline fetches
    each gene's pathways and partners itself, deterministically, before the model
    sees anything. That split -- what the pipeline computes versus what the model
    chooses -- is the point of the design, so it is shown rather than flattened.

    `origin` is recorded at capture time from the store's ordering, not inferred
    here: the model frequently requests a call the pipeline already made, and no
    amount of argument matching can tell those apart after the fact.
    """
    rows: list[list[str]] = []
    for call in replay.calls:
        result = call.get("result", {})
        records = result.get("record_ids", [])
        notes = result.get("notes", [])
        rows.append([
            call.get("origin", "model"),
            call["tool"],
            ", ".join(f"{k}={v}" for k, v in call["arguments"].items()),
            f"{len(records)} record(s)" + (f": {', '.join(records[:4])}" if records else ""),
            notes[0][:160] if notes else "",
        ])
    return rows


def citation_rows(replay: Replay) -> list[list[str]]:
    """Every identifier the summary cites, and whether it survives checking."""
    rows = []
    for citation in replay.report.citations:
        record = replay.records.get(citation.identifier, {})
        rows.append([
            citation.identifier,
            record.get("type", ""),
            {"verified": "verified",
             "unsupported": "UNSUPPORTED",
             "cross_target": "CROSS-TARGET"}.get(citation.status, citation.status),
            citation.detail or (record.get("name", "")[:70]),
        ])
    return rows


def quote_rows(replay: Replay) -> list[list[str]]:
    return [[q.record_id,
             "verified" if q.status == "verified" else "NOT IN SOURCE",
             q.quote[:120],
             "" if q.status == "verified" else (q.closest_span[:120] or "—")]
            for q in replay.report.quotes]
