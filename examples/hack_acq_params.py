import sys
import h5py
from PyQt5.QtWidgets import (
    QApplication,
    QWidget,
    QVBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QMessageBox,
    QHBoxLayout,
)
from pyspecdata import search_filename


class EditAcqParams(QWidget):
    def __init__(self, hdf_filename, nodename):
        super().__init__()
        self.hdf_filename = hdf_filename
        self.nodename = nodename
        self.acq_params = {}

        # Load HDF5 data
        self.load_hdf5_data()

        # Define property names and labels
        self.property_names = [
            "coherence_pathway",
            "chemical",
            "concentration",
            "postproc_type",
        ]
        self.property_labels = [
            "Coherence Pathway",
            "Chemical Name",
            "Concentration",
            "Postproc Type",
        ]

        # Create UI dynamically
        self.init_ui()

    def load_hdf5_data(self):
        with h5py.File(self.hdf_filename, "r") as hdf_file:
            if self.nodename in hdf_file:
                node = hdf_file[self.nodename]
                if "other_info" in node and "acq_params" in node["other_info"]:
                    self.acq_params = node["other_info"]["acq_params"]
                    print("the string attributes will appear here:")
                    for k,v in self.acq_params.attrs.items():
                        print(k,v)
                else:
                    raise ValueError(
                        "acq_params not found in"
                        f" '{self.nodename}/other_info'."
                    )
            else:
                raise ValueError(
                    f"Expno '{self.nodename}' not found in the HDF5 file."
                )

    def init_ui(self):
        layout = QVBoxLayout()

        # Store input fields for easy access later
        self.input_fields = {}

        # Dynamically create labeled text fields
        for prop_name, label_text in zip(
            self.property_names, self.property_labels
        ):
            h_layout = QHBoxLayout()
            label = QLabel(label_text)
            text_field = QLineEdit(str(self.acq_params.get(prop_name, "")))
            h_layout.addWidget(label)
            h_layout.addWidget(text_field)
            layout.addLayout(h_layout)
            self.input_fields[prop_name] = text_field

        # Save button
        save_button = QPushButton("Save Changes")
        save_button.clicked.connect(self.save_changes)
        layout.addWidget(save_button)

        self.setLayout(layout)
        self.setWindowTitle("Edit Acq Params")
        self.show()

    def save_changes(self):
        try:
            with h5py.File(self.hdf_filename, "a") as hdf_file:
                node = hdf_file[self.nodename]
                acq_params = node["other_info"]["acq_params"].attrs

                for prop_name in self.property_names:
                    acq_params[prop_name] = self.input_fields[prop_name].text()

                QMessageBox.information(
                    self, "Success", "Values saved successfully!"
                )
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))


def main():
    if len(sys.argv) != 4:
        print("Usage: python hack_acq_params.py <hdf_filename> <nodename>")
        sys.exit(1)

    hdf_filename = search_filename(
        sys.argv[2], unique=True, exp_type=sys.argv[3]
    )
    nodename = sys.argv[1]

    app = QApplication(sys.argv)
    editor = EditAcqParams(hdf_filename, nodename)
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
