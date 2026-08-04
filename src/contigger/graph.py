"""Graph boundary for a later conservative implementation.

See DESIGN.md sections 14 and 15. No traversal or simplification is implemented
until ambiguity preservation rules have executable tests.
"""

from contigger.exceptions import FeatureNotImplementedError
from contigger.models import MergeComponent, Relationship


def build_components(relationships: tuple[Relationship, ...]) -> tuple[MergeComponent, ...]:
    """Reserve graph construction without returning plausible biological results."""
    raise FeatureNotImplementedError("overlap graph construction is not implemented")
