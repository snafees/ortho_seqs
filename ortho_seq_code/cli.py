import click
from ortho_seq_code.orthogonal_polynomial import cli as orthogonal_polynomial
from ortho_seq_code.gui.gui import gui_run

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group()
# @click.pass_context
def cli():
    pass


cli.add_command(orthogonal_polynomial, name="orthogonal-polynomial")
# cli.add_command(gui_run, name="gui")


if __name__ == "__main__":
    cli()
