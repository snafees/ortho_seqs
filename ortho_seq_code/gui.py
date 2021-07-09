from PyQt5.QtWidgets import (
    QApplication,
    QLabel,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QComboBox,
    QPushButton,
    QHBoxLayout,
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

        upload_layout = QHBoxLayout()
        upload_layout.addWidget(label2)
        upload_layout.addWidget(label3)
        upload_layout.addWidget(upload_button_1)
        upload_layout.addWidget(upload_button_2)
        styleComboBox = QComboBox()
        styleComboBox.addItems(["DNA", "protein", "dna_n", "protein_n", "protein_pnp"])
        styleLabel = QLabel("&Molecule:")
        styleLabel.setBuddy(styleComboBox)

        self.widget_layout.addWidget(label)
        self.widget_layout.addLayout(upload_layout)
        self.widget_layout.addWidget(styleLabel)
        self.widget_layout.addWidget(styleComboBox)
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
