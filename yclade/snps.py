"""Tools for working with SNPs"""

from yclade.types import SnpResults


def parse_snp_results(snp_string: str) -> SnpResults:
    """Parse a string of comma separate SNPs into SnpResults."""
    snps = [snp.strip() for snp in snp_string.split(",")]
    positive_snps = {snp.rstrip("+") for snp in snps if snp.endswith("+")}
    negative_snps = {snp.rstrip("-") for snp in snps if snp.endswith("-")}
    return SnpResults(positive=positive_snps, negative=negative_snps)
