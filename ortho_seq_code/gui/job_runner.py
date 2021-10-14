from PyQt5.QtWidgets import QWidget

from ortho_seq_code.orthogonal_polynomial import orthogonal_polynomial
from ortho_seq_code.gui.worker import Worker


class JobRunner(QWidget):
    def __init__(self, parent, threadpool):
        super(QWidget, self).__init__(parent)
        self.parent = parent
        self.threadpool = threadpool

    def job_func(self, progress_callback):
        orthogonal_polynomial(
            self.filename,
            self.pheno_file,
            self.molecule,
            self.poly_order,
            self.precomputed,
            self.out_dir,
            self.alphbt_input,
            self.min_pct
        )

    def thread_complete(self):
        print("Thread COMPLETE")

    def launch(self):
        self.filename = self.parent.upload_button_1.text()
        self.pheno_file = self.parent.upload_button_2.text()
        self.molecule = self.parent.styleComboBox1.currentText()
        self.poly_order = self.parent.poly
        self.precomputed,
        self.out_dir,
        self.alphbt_input,
        self.min_pct

        worker = Worker(self.job_func)

        worker.signals.finished.connect(self.thread_complete)

        # Execute
        self.threadpool.start(worker)
