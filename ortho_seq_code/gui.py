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
)
import sys


class MainWidget(QWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.widget_layout = QVBoxLayout()  # vertical layout
        label = QLabel(
            "ortho_seqs GUI: User interface to compute tensor-based orthogonal polynomials for sequence data"
        )
        label2 = QLabel("Upload sequence file:")
        upload_button_1 = QPushButton("seq_file")
        label3 = QLabel("Upload phenotype file:")
        upload_button_2 = QPushButton("pheno_file")
        start_button = QPushButton("RUN")

        upload_layout1 = QHBoxLayout()
        upload_layout1.addWidget(label2)
        upload_layout2 = QHBoxLayout()
        upload_layout2.addWidget(label3)
        upload_layout1.addWidget(upload_button_1)
        upload_layout2.addWidget(upload_button_2)

        upload_ComboBox1_layout = QHBoxLayout()
        styleComboBox1 = QComboBox()
        styleComboBox1.addItems(["DNA", "protein", "dna_n", "protein_n", "protein_pnp"])
        styleLabel1 = QLabel("&Molecule:")
        styleLabel1.setBuddy(styleComboBox1)
        upload_ComboBox1_layout.addWidget(styleComboBox1)

        upload_ComboBox2_layout = QHBoxLayout()
        styleComboBox2 = QComboBox()
        styleComboBox2.addItems(
            ["first", "second"]
        )  # only doing first and second order polynomials so far
        styleLabel2 = QLabel("&poly_order:")
        styleLabel2.setBuddy(styleComboBox2)
        upload_ComboBox2_layout.addWidget(styleComboBox2)

        label4 = QLabel("Precomputed?")
        upload_layout3 = QHBoxLayout()
        upload_layout3.addWidget(label4)
        checkbox1 = QCheckBox("Yes")
        checkbox2 = QCheckBox("No")

        self.widget_layout.addWidget(label)
        self.widget_layout.addLayout(upload_layout1)
        self.widget_layout.addLayout(upload_layout2)

        self.widget_layout.addWidget(styleLabel1)  # for moledule
        self.widget_layout.addWidget(styleComboBox1)
        self.widget_layout.addLayout(upload_ComboBox1_layout)

        self.widget_layout.addWidget(styleLabel2)  # for poly_order
        self.widget_layout.addWidget(styleComboBox2)
        self.widget_layout.addLayout(upload_ComboBox2_layout)

        self.widget_layout.addLayout(upload_layout3)  # for checkbox
        self.widget_layout.addWidget(checkbox1)
        self.widget_layout.addWidget(checkbox2)

        self.widget_layout.addWidget(start_button)
        self.setLayout(self.widget_layout)


class MainWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.setCentralWidget(MainWidget())
        self.originalPalette = QApplication.palette()
        self.show()


def main():
    app = QApplication(sys.argv)
    window = MainWindow()

    sys.exit(app.exec())  # instead of just app.exec()


if __name__ == "__main__":
    main()
