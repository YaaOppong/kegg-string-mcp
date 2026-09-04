"""LangGraph orchestration for gene-pair evidence retrieval.

A straight chunk -> embed -> query pipeline does not need a graph framework, and
wrapping one in LangGraph would be scaffolding. What earns it is the **cycle**:

    retrieve -> judge sufficiency -> rewrite the query -> retrieve again

That is corrective retrieval. It has state (the queries already tried, the hits
accumulated), a branch (sufficient or not), and a loop with a bound. Those are the
things a graph framework exists for, and none of them are present in a one-shot
search.

The sufficiency judge is deterministic by default, and that is deliberate. The
repo already computes, per retrieved paper, which of the query's genes the text
actually **mentions** -- so "did this retrieval surface papers discussing both
genes of the pair?" is answerable without a model, without an API key, and
without a judgement anyone has to trust. A model-based judge can be injected, but
the default keeps the whole graph runnable and testable offline.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field
from typing import Any, TypedDict

from langgraph.graph import END, START, StateGraph

from kegg_string_mcp.retrieval.index import DEFAULT_K, Hit

MAX_ROUNDS = 3


class GraphState(TypedDict, total=False):
    genes: list[str]
    query: str
    tried: list[str]
    hits: list[Hit]
    rounds: int
    sufficient: bool
    reason: str
    no_new_hits: bool


def papers_naming_all(hits: list[Hit], genes: list[str]) -> list[Hit]:
    """Passages whose own text names EVERY gene in the pair.

    The union across hits is the wrong test and made the cycle unreachable: the
    corpus is built by querying each gene separately, so eight hits almost always
    mention both genes between them while no single paper discusses the pair.
    Evidence for a *pair* is one paper that names both.

    Scored on `genes_named` -- what the passage's text names, computed over the
    whole corpus -- not on `mentions`, which records only the terms of the query
    that fetched the paper. Judging on `mentions` made the loop reject evidence
    sitting in its own hits.
    """
    wanted = {g.lower() for g in genes}
    return [h for h in hits if wanted <= {m.lower() for m in h.genes_named}]


def genes_covered(hits: list[Hit], genes: list[str]) -> list[str]:
    """Genes named somewhere in the retrieved text, across all hits."""
    seen: set[str] = set()
    for hit in hits:
        seen.update(m.lower() for m in hit.genes_named)
    return [g for g in genes if g.lower() in seen]


@dataclass
class RetrievalGraph:
    """Corrective retrieval over any retriever exposing `.search(query, k)`."""

    retriever: Any
    k: int = DEFAULT_K
    max_rounds: int = MAX_ROUNDS
    judge: Callable[[list[Hit], list[str]], tuple[bool, str]] | None = None
    rewrite: Callable[[list[str], list[str]], str] | None = None
    trace: list[dict[str, Any]] = field(default_factory=list)

    # -- nodes --------------------------------------------------------------

    def _retrieve(self, state: GraphState) -> GraphState:
        hits = self.retriever.search(state["query"], k=self.k)
        # Accumulate across rounds: a rewritten query that finds the missing gene
        # should add to the evidence, not replace what round one already found.
        merged = {h.passage_id: h for h in state.get("hits", [])}
        for hit in hits:
            merged.setdefault(hit.passage_id, hit)
        self.trace.append({"node": "retrieve", "round": state.get("rounds", 0) + 1,
                           "query": state["query"], "new_hits": len(hits),
                           "total_hits": len(merged)})
        added = len(merged) - len(state.get("hits", []))
        return {"hits": list(merged.values()),
                "tried": [*state.get("tried", []), state["query"]],
                "rounds": state.get("rounds", 0) + 1,
                "no_new_hits": state.get("rounds", 0) > 0 and added == 0}

    def _judge(self, state: GraphState) -> GraphState:
        judge = self.judge or self._default_judge
        ok, reason = judge(state["hits"], state["genes"])
        self.trace.append({"node": "judge", "round": state["rounds"],
                           "sufficient": ok, "reason": reason})
        return {"sufficient": ok, "reason": reason}

    def _rewrite(self, state: GraphState) -> GraphState:
        covered = genes_covered(state["hits"], state["genes"])
        missing = [g for g in state["genes"] if g not in covered]
        if self.rewrite is not None:
            query = self.rewrite(state["genes"], missing)
        else:
            query = self._default_rewrite(state["genes"], missing, state["rounds"])
        self.trace.append({"node": "rewrite", "round": state["rounds"],
                           "missing": missing, "new_query": query})
        return {"query": query}

    # -- defaults, deterministic and offline --------------------------------

    @staticmethod
    def _default_judge(hits: list[Hit], genes: list[str]) -> tuple[bool, str]:
        joint = papers_naming_all(hits, genes)
        if joint:
            return True, (f"{len(joint)} passage(s) name all of {', '.join(genes)} "
                          f"(e.g. PMID {joint[0].pmid})")
        covered = genes_covered(hits, genes)
        missing = [g for g in genes if g not in covered]
        if missing:
            return False, f"nothing retrieved names {', '.join(missing)}"
        return False, (f"all of {', '.join(genes)} appear, but never in the same paper -- "
                       f"no joint evidence for the pair")

    @staticmethod
    def _default_rewrite(genes: list[str], missing: list[str], round_number: int = 1) -> str:
        """Narrow onto what is missing, and vary between rounds.

        A rewrite that returns the same string spins the loop for nothing -- the
        first version did exactly that, issuing an identical query twice and adding
        zero hits. Each round reframes: co-occurrence first, then mechanism.
        """
        focus = missing or genes
        joined = " ".join(focus)
        angles = [
            f"{joined} interaction regulation Mycobacterium tuberculosis",
            f"{joined} co-occurrence combined mutation epistasis compensatory",
            f"{joined} pathway mechanism functional relationship",
        ]
        return angles[min(round_number - 1, len(angles) - 1)]

    # -- graph --------------------------------------------------------------

    def build(self):
        graph = StateGraph(GraphState)
        graph.add_node("retrieve", self._retrieve)
        graph.add_node("judge", self._judge)
        graph.add_node("rewrite", self._rewrite)
        graph.add_edge(START, "retrieve")
        graph.add_edge("retrieve", "judge")
        graph.add_conditional_edges("judge", self._route,
                                    {"rewrite": "rewrite", END: END})
        graph.add_edge("rewrite", "retrieve")   # the cycle
        return graph.compile()

    def _route(self, state: GraphState) -> str:
        if state.get("sufficient"):
            return END
        if state.get("no_new_hits"):
            self.trace.append({"node": "route", "round": state["rounds"],
                               "stopped": "a rewrite returned no passages the previous "
                                          "rounds had not already found"})
            return END
        if state.get("rounds", 0) >= self.max_rounds:
            # Bounded: give up and say so, rather than looping on a corpus that
            # simply does not contain the answer. "Nothing found" is a result.
            self.trace.append({"node": "route", "round": state["rounds"],
                               "stopped": "max_rounds reached"})
            return END
        return "rewrite"

    def run(self, genes: list[str], query: str | None = None) -> dict[str, Any]:
        self.trace = []
        opening = query or (f"What is the relationship between {' and '.join(genes)} "
                            f"in Mycobacterium tuberculosis?")
        final = self.build().invoke({"genes": genes, "query": opening,
                                     "tried": [], "hits": [], "rounds": 0})
        return {"genes": genes, "rounds": final["rounds"], "queries": final["tried"],
                "sufficient": final.get("sufficient", False),
                "reason": final.get("reason", ""),
                "hits": final["hits"], "trace": list(self.trace)}
