import sys
import click
from cli.prepare import prepare
from cli.apamix import apamix
from cli.extract import extract

@click.group()
def cli():
    """
    A pure python package for analysing alternative polyadenylation at single cell levels.
    Current version: 1.0.0

    """
    pass

cli.add_command(prepare)
cli.add_command(apamix)
cli.add_command(extract)