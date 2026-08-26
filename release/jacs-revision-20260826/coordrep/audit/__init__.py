"""Release-scale, source-faithful CoordRep audit records.

The audit layer is intentionally distinct from a fully resolved
``CoordRep-State``.  It gives every in-domain metal site a typed structural
record and preserves partial or ambiguous source information without forcing
an unsupported chemical assignment.
"""

from .records import (
    SCHEMA_VERSION,
    CoordinationIncidence,
    DonorGroupRelation,
    DonorPairRelation,
    DonorGroup,
    EntryRecord,
    ExternalCenterRelation,
    ExternalDonorRelation,
    MetalRelation,
    MetalSiteRecord,
    canonical_json,
    record_digest,
    validate_entry_record,
)

__all__ = [
    "SCHEMA_VERSION",
    "CoordinationIncidence",
    "DonorGroupRelation",
    "DonorPairRelation",
    "DonorGroup",
    "EntryRecord",
    "ExternalCenterRelation",
    "ExternalDonorRelation",
    "MetalRelation",
    "MetalSiteRecord",
    "canonical_json",
    "record_digest",
    "validate_entry_record",
]
