"""Tests for the yclade.snps module."""

import yclade.snps


def test_parse_snp_results_with_white_space():
    """Test the parse_snp_results function with white space present."""
    snp_string = "A+, Bb-, Ccc-, D+"
    results = yclade.snps.parse_snp_results(snp_string)
    assert results.positive == {"A", "D"}
    assert results.negative == {"Bb", "Ccc"}


def test_parse_snp_results_without_white_space():
    """Test the parse_snp_results function without white space present."""
    snp_string = "A+,Bb-,Ccc-,D+"
    results = yclade.snps.parse_snp_results(snp_string)
    assert results.positive == {"A", "D"}
    assert results.negative == {"Bb", "Ccc"}


def test_parse_snp_results_empty():
    """Test the parse_snp_results function with an empty string."""
    snp_string = ""
    results = yclade.snps.parse_snp_results(snp_string)
    assert results.positive == set()
    assert results.negative == set()

def test_parse_snp_results_invalid():
    """Test the parse_snp_results function with an invalid string."""
    snp_string = "invalid"
    results = yclade.snps.parse_snp_results(snp_string)
    assert results.positive == set()
    assert results.negative == set()
