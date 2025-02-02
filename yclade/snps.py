"""Tools for working with SNPs"""

import logging
from yclade.types import SnpResults, SnpAliases


def parse_snp_results(snp_string: str) -> SnpResults:
    """Parse a string of comma separated SNPs into SnpResults.

    Args:
        snp_string: A string of comma separated SNPs.

    The SNPs can be separated by commas and, optionally, spaces.
    All SNPs must have the form snp+ (for a positively tested SNP)
    or snp- (for a negatively tested SNP), otherwise they are ignored.
    """
    snps = [snp.strip() for snp in snp_string.split(",")]
    positive_snps = {snp.rstrip("+") for snp in snps if snp.endswith("+")}
    negative_snps = {snp.rstrip("-") for snp in snps if snp.endswith("-")}
    snp_results = SnpResults(positive=positive_snps, negative=negative_snps)
    return warn_and_remove_duplicates(snp_results)


def warn_and_remove_duplicates(snp_results: SnpResults) -> SnpResults:
    """Warn about and remove duplicate SNPs from the results."""
    overlap = snp_results.positive & snp_results.negative
    if overlap:
        logger = logging.getLogger("yclade")
        logger.warning("Ignoring SNPs with conflicting test results: %s", overlap)
        return SnpResults(
            positive=snp_results.positive - overlap,
            negative=snp_results.negative - overlap,
        )
    return snp_results


def normalize_snp_results(
    snp_results: SnpResults, snp_aliases: SnpAliases
) -> SnpResults:
    """Normalize the SNP results by casting aliases to the canonical form."""
    snp_results = SnpResults(
        positive={snp_aliases.get(snp, snp) for snp in snp_results.positive},
        negative={snp_aliases.get(snp, snp) for snp in snp_results.negative},
    )
    return warn_and_remove_duplicates(snp_results)
