"""Rule annotation: read a rule learner's output and say what is already known.

The upstream analysis emits rules over per-locus variant states predicting a
resistance phenotype. This package resolves the features those rules are written
in, then -- in later stages -- classifies each rule by what the structured
sources already account for: a known resistance locus, a known compensatory one,
an alternative route taken when the canonical locus is at reference, or nothing
known at all.

The division of labour is the one the rest of the package uses. What the
catalogue, the annotation and the interaction data can settle is settled by code.
The model is reserved for rules those sources leave open, and everything it says
is checked against what the tools returned.
"""

from kegg_string_mcp.rules.annotation import Annotation, Gene, Intergenic
from kegg_string_mcp.rules.annotation import parse as parse_annotation
from kegg_string_mcp.rules.catalogue import ABSENT, ANCHOR, ASSESSED_NEGATIVE, Catalogue
from kegg_string_mcp.rules.evidence import Sources, classify_all, rename_map, summarise
from kegg_string_mcp.rules.features import (
                                            CODING,
                                            INTERGENIC,
                                            UNRESOLVED,
                                            Coverage,
                                            Feature,
                                            Resolver,
                                            resolve_all,
)
from kegg_string_mcp.rules.parse import Condition, Rule, vocabulary
from kegg_string_mcp.rules.parse import parse as parse_rules
from kegg_string_mcp.rules.signature import Link, RuleSignature, classify

__all__ = [
                                            "ABSENT",
                                            "ANCHOR",
                                            "ASSESSED_NEGATIVE",
                                            "CODING",
                                            "INTERGENIC",
                                            "UNRESOLVED",
                                            "Annotation",
                                            "Catalogue",
                                            "Condition",
                                            "Coverage",
                                            "Feature",
                                            "Gene",
                                            "Intergenic",
                                            "Link",
                                            "Resolver",
                                            "Rule",
                                            "RuleSignature",
                                            "Sources",
                                            "classify",
                                            "classify_all",
                                            "parse_annotation",
                                            "parse_rules",
                                            "rename_map",
                                            "resolve_all",
                                            "summarise",
                                            "vocabulary",
]
