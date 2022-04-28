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
            filename=self.filename,
            pheno_file=self.pheno_file,
            molecule=self.molecule,
            poly_order=self.poly_order,
            precomputed=self.precomputed,
            alphbt_input=self.alphbt_input,
            out_dir="../results_ortho_seq_testing/",
            min_pct=self.min_pct,
            pheno_name=self.pheno_name,
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
        self.alphbt_input = self.parent.alphabet_text.text()
        self.min_pct = int(self.parent.pct_text.text())
        self.pheno_name = self.parent.pheno_text.text()

        worker = Worker(self.job_func)

        worker.signals.finished.connect(self.thread_complete)

        # Execute
        self.threadpool.start(worker)
