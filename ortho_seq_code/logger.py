import sys
import os


class Logger(object):
    def __init__(self, out_dir=None, name="cli_output"):
        self.terminal = sys.stdout
        self.log = open(str(os.path.join(out_dir, name + ".txt")), "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass
