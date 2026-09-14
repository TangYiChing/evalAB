"""Tier 1 (sequence-only) antibody pre-screening pipeline.

See docs/antibody-prescreen-agent-plan.md for the design and rationale.
Entry points: screen_candidate() and screen_batch() in .fusion.
"""

from .fusion import CandidateResult, format_report, screen_batch, screen_candidate

__all__ = ["screen_candidate", "screen_batch", "format_report", "CandidateResult"]
