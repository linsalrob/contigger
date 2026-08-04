"""Domain-specific exceptions exposed by Contigger."""


class ContiggerError(Exception):
    """Base class for expected, user-facing Contigger errors."""


class ConfigurationError(ContiggerError):
    """Raised when run configuration is inconsistent or unsupported."""


class InputValidationError(ContiggerError):
    """Raised when a manifest or biological input is invalid."""


class FastaFormatError(InputValidationError):
    """Raised when a FASTA stream violates the supported syntax."""


class ManifestError(InputValidationError):
    """Raised when a sample manifest is invalid."""


class ExternalToolError(ContiggerError):
    """Raised when an external tool is missing or fails."""


class FeatureNotImplementedError(ContiggerError):
    """Raised when a deliberately unsupported biological operation is requested."""
