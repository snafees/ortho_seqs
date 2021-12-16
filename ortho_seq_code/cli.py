import click
from ortho_seq_code.orthogonal_polynomial import cli as orthogonal_polynomial
from ortho_seq_code.gui.gui import gui_run
from ortho_seq_code.plotclass import plotclass


CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(invoke_without_command=True, context_settings=CONTEXT_SETTINGS)
@click.pass_context
def cli(ctx):
    if ctx.invoked_subcommand is None:
        gui_run()
    else:
        pass


cli.add_command(orthogonal_polynomial, name="orthogonal-polynomial")

cli.add_command(plotclass, name="rf1d")

if __name__ == "__main__":
    cli()
