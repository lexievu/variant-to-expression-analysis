"""Custom exceptions for the variant-to-expression pipeline.

Using exceptions instead of ``sys.exit()`` makes every pipeline function
testable and composable — callers can catch specific errors rather than
having the process killed out from under them.
"""


class PipelineInputError(Exception):
    """Raised when a required pipeline input is missing or invalid.

    Examples: missing VCF file, missing API key, malformed gene ID.
    """
