from PyQt5.QtWidgets import QWidget

from ortho_seq_code.orthogonal_polynomial import orthogonal_polynomial
from ortho_seq_code.plotclass import rf1d
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
            out_dir="../results_ortho_seq_testing/",
            self.alphbt_input,
            min_pct=75,
            pheno_name=None,
        )

    def thread_complete(self):
        print("Thread COMPLETE")

    def launch(self):
        self.filename = self.parent.upload_button_1.text()
        self.pheno_file = self.parent.upload_button_2.text()
        self.molecule = self.parent.molecule_combobox.currentText()
        self.poly_order = self.parent.poly_order_combobox.currentText()
        # self.precomputed = self.parent.precomputed_combobox.currentText() == "Yes"
        self.precomputed = self.parent.upload_button_3.text()
        self.alphabet_input = self.parent.alphabet_text.text()

        worker = Worker(self.job_func)

        worker.signals.finished.connect(self.thread_complete)

        # Execute
        self.threadpool.start(worker)
