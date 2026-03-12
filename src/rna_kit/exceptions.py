from __future__ import annotations


class RNAAssessmentError(Exception):
    """Base class for project-specific exceptions."""


class InvalidResidueRangeError(RNAAssessmentError):
    """Raised when a residue range specification is malformed."""


class NormalizationError(RNAAssessmentError):
    """Raised when a PDB file cannot be normalized safely."""


class SequenceMismatchError(RNAAssessmentError):
    """Raised when two indexed structures do not describe the same sequence."""


class ToolNotAvailableError(RNAAssessmentError):
    """Raised when an optional third-party tool cannot be resolved or executed."""


class ToolResolutionError(ToolNotAvailableError):
    """Raised when an optional third-party tool cannot be located."""


class ToolExecutionError(ToolNotAvailableError):
    """Raised when an optional third-party tool is found but execution fails."""


class MetricCalculationError(RNAAssessmentError):
    """Raised when a metric cannot be computed from the provided structures."""


class ManifestFormatError(RNAAssessmentError):
    """Raised when a benchmark manifest is missing required fields or uses an unsupported format."""


class ReportGenerationError(RNAAssessmentError):
    """Raised when a report or visualization artifact cannot be generated."""


class SchemaValidationError(RNAAssessmentError):
    """Raised when a result document does not satisfy the exported schema contract."""
