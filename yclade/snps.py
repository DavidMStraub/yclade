"""Tools for working with SNPs"""

from yclade.types import SnpResults


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
    return SnpResults(positive=positive_snps, negative=negative_snps)
