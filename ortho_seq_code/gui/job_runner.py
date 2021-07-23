from PyQt5.QtWidgets import QWidget
import os

from ortho_seq_code.orthogonal_polynomial import orthogonal_polynomial
from ortho_seq_code.gui.worker import Worker


class JobRunner(QWidget):
    def __init__(self, parent, threadpool):
        super(QWidget, self).__init__(parent)
        self.parent = parent
        self.threadpool = threadpool

        self.result_images = []

    def job_func(self, progress_callback):

        orthogonal_polynomial(
            self.filename,
            self.pheno_file,
            self.molecule,
            2,
            4,
            1,
            "first",
            False,
            os.getcwd(),
        )

        # Initialize required restoration instances

        return None

    def thread_complete(self):
        print("Thread COMPLETE")

    def launch(self):
        self.filename = self.parent.upload_button_1.text()
        self.pheno_file = self.parent.upload_button_2.text()
        self.molecule = self.parent.styleComboBox1.currentText()


        # Pass the function to execute
        worker = Worker(
            self.job_func
        )  # Any other args, kwargs are passed to the run function

        worker.signals.finished.connect(self.thread_complete)

        # Execute
        self.threadpool.start(worker)
