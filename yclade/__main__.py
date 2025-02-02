"""Command line interface for yclade."""

import click
import logging

import yclade


@click.command()
@click.argument("snp_string", required=False)
@click.option(
    "--version", "-V", default=None, help="The YFull version to use (optional)."
)
@click.option(
    "--data-dir", default=None, help="The directory to store the YFull data (optional)."
)
@click.option("--file", "-f", type=click.File("r"))
@click.option("--verbose", "-v", is_flag=True, help="Print verbose output.")
def main(snp_string, version, data_dir, file, verbose):
    """Find the clade that matches the given SNP string."""
    logger = logging.getLogger("yclade")
    console_handler = logging.StreamHandler()
    formatter = logging.Formatter("%(levelname)s:%(name)s:%(message)s")
    console_handler.setFormatter(formatter)
    if verbose:
        logger.setLevel(logging.DEBUG)
        console_handler.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.WARNING)
        console_handler.setLevel(logging.DEBUG)
    logger.addHandler(console_handler)
    if not snp_string and not file:
        raise click.UsageError("You must provide either a SNP string or a file.")
    if file:
        snp_string = file.read().strip()
    clade_info = yclade.find_clade(snp_string, version=version, data_dir=data_dir)
    for info in clade_info:
        click.echo(info)


if __name__ == "__main__":
    main()
