import click
from ortho_seq_code.orthogonal_polynomial import cli as orthogonal_polynomial


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    pass


cli.add_command(orthogonal_polynomial, name="orthogonal-polynomial")
