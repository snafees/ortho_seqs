import sys
import os


class Logger(object):
    def __init__(self, out_dir=None):
        self.terminal = sys.stdout
        print(out_dir)
        self.log = open(str(os.path.join(out_dir, "cli_output.txt")), "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        # this flush method is needed for python 3 compatibility.
        # this handles the flush command by doing nothing.
        # you might want to specify some extra behavior here.
        pass
