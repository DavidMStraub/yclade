"""Command line interface for yclade."""

import click

import yclade


@click.command()
@click.argument("snp_string", required=False)
@click.option("--version", default=None, help="The YFull version to use (optional).")
@click.option(
    "--data-dir", default=None, help="The directory to store the YFull data (optional)."
)
@click.option("--file", "-f", type=click.File("r"))
def main(snp_string, version, data_dir, file):
    """Find the clade that matches the given SNP string."""
    if not snp_string and not file:
        raise click.UsageError("You must provide either a SNP string or a file.")
    if file:
        snp_string = file.read().strip()
    clade_info = yclade.find_clade(snp_string, version=version, data_dir=data_dir)
    for info in clade_info:
        click.echo(info)


if __name__ == "__main__":
    main()
