# yclade

A Python library to find the position on the human Y chromosome tree given a set of single nucleotide polymorphisms (SNPs).

## Installation

```bash
pip install yclade
```

## Usage

```python
import yclade

snp_string = "M215+, BY61636-, FTF15749-, TY15744-"
yclade.find_clade(snp_string)
```

## Command-line usage

```bash
python -m yclade "M215+, BY61636-, FTF15749-, TY15744-"
```

You can also use a file with the SNPs:

```bash
python -m yclade -f snps.txt
```

## Credits

The Y tree data is from [YFull](https://www.yfull.com/) and [is shared](https://github.com/YFullTeam/YTree) under the [Creative Commons Attribution 4.0 International Public License](https://github.com/YFullTeam/YTree/blob/master/LICENSE.md).

The code has been tested on SNP data obtained from [YSEQ](https://www.yseq.net/) test results.

Inspiration has also been drawn from [clade-finder](https://github.com/hprovyn/clade-finder).
