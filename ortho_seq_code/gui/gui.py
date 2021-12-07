from PyQt5.QtCore import QThreadPool
from PyQt5.QtWidgets import (
    QApplication,
    QLabel,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QComboBox,
    QPushButton,
    QHBoxLayout,
    QCheckBox,
    QFileDialog,
)
import sys
import click


from ortho_seq_code.gui.job_runner import JobRunner


class MainWidget(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.widget_layout = QVBoxLayout()  # vertical layout

        # Title
        label = QLabel(
            "User interface to compute tensor-based orthogonal polynomials for sequence data"
        )
        self.widget_layout.addWidget(label)

        # Upload file buttons
        label2 = QLabel("Upload sequence file:")
        self.upload_button_1 = QPushButton("seq_file")
        self.upload_button_1.clicked.connect(
            lambda: self.upload_button_1.setText(self.openFileNamesDialog())
        )
        label3 = QLabel("Upload phenotype file:")
        self.upload_button_2 = QPushButton("pheno_file")
        self.upload_button_2.clicked.connect(
            lambda: self.upload_button_2.setText(self.openFileNamesDialog())
        )

        upload_layout1 = QHBoxLayout()
        upload_layout1.addWidget(label2)
        upload_layout2 = QHBoxLayout()
        upload_layout2.addWidget(label3)
        upload_layout1.addWidget(self.upload_button_1)
        upload_layout2.addWidget(self.upload_button_2)
        self.widget_layout.addLayout(upload_layout1)
        self.widget_layout.addLayout(upload_layout2)

        # molecule combobox
        upload_ComboBox1_layout = QHBoxLayout()
        self.molecule_combobox = QComboBox()
        self.molecule_combobox.addItems(
            ["DNA", "protein", "dna_n", "protein_n", "protein_pnp"]
        )
        styleLabel1 = QLabel("&Molecule:")
        styleLabel1.setBuddy(self.molecule_combobox)
        upload_ComboBox1_layout.addWidget(self.molecule_combobox)
        self.widget_layout.addWidget(styleLabel1)  # for molecule
        self.widget_layout.addWidget(self.molecule_combobox)
        self.widget_layout.addLayout(upload_ComboBox1_layout)

        # poly_order combobox
        upload_ComboBox2_layout = QHBoxLayout()
        self.poly_order_combobox = QComboBox()
        self.poly_order_combobox.addItems(
            ["first", "second"]
        )  # only doing first and second order polynomials so far
        styleLabel2 = QLabel("&poly_order:")
        styleLabel2.setBuddy(self.poly_order_combobox)
        upload_ComboBox2_layout.addWidget(self.poly_order_combobox)
        self.widget_layout.addWidget(styleLabel2)  # for poly_order
        self.widget_layout.addWidget(self.poly_order_combobox)
        self.widget_layout.addLayout(upload_ComboBox2_layout)

        # Precomputed checkboxes
        #upload_layout3 = QHBoxLayout()
        # label4 = QLabel("Precomputed?")
        #upload_layout3.addWidget(label4)
        # self.precomputed_combobox = QComboBox()
        # self.precomputed_combobox.addItems(["No", "Yes"])
        # upload_layout3.addWidget(self.precomputed_combobox)
        # self.widget_layout.addLayout(upload_layout3)  # for checkbox

        # precomputed dir path
        # Upload file buttons
        upload_layout3 = QHBoxLayout()
        label4 = QLabel("Select directory with precomputed sequence file:")
        upload_layout3.addWidget(label4)
        self.upload_button_3 = QPushButton("precomputed_dir")
        self.upload_button_3.clicked.connect(
            lambda: self.upload_button_3.setText(
                self.getExistingDirectory(self, "Select Directory")
            )
        )
        # self.upload_button_3.clicked.connect(
        #     lambda: self.upload_button_3.setText(self.openFileNamesDialog())
        # )

        # upload_layout3.addWidget(label4)
        # upload_layout3.addWidget(self.upload_button_3)
        self.widget_layout.addLayout(upload_layout3)


        # RUN button
        start_button = QPushButton("RUN")
        self.threadpool = QThreadPool(self)
        self.job_runner = JobRunner(self, self.threadpool)
        start_button.clicked.connect(self.job_runner.launch)
        self.widget_layout.addWidget(start_button)

        self.setLayout(self.widget_layout)

    def openFileNamesDialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        files, _ = QFileDialog.getOpenFileNames(
            self, "Open File(s)", "", "All Files (*)", options=options
        )

        if files:
            return files[0]


class MainWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setCentralWidget(MainWidget())
        self.show()


# @click.command(help="GUI option instead of command line input")


@click.command()
# @click.option("--gui", help="to run gui")
def gui_run():
    app = QApplication(sys.argv)
    window = MainWindow()

    sys.exit(app.exec())  # instead of just app.exec()


# if __name__ == "__main__":
#     gui_run()
