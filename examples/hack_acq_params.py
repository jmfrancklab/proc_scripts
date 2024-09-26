"""
Hack H5 Properties
==================

There are various labeling properties that are stored in the HDF5 file.
Sometimes we mess up, and we need to hack these properties ex-post-facto.

Note that is intended *only* for hacking labeling/categorizing properties, and
should *never* be used to hack information about pulse program parameters, etc,
that were actually sent to the spectrometer.

The GPT conversation is `here <https://chatgpt.com/c/66f5aef7-8f68-800b-8a13-0b133ede1ca5>`_.
"""
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
        self.hdf_file = None
        self.node = None

        # Define property names and labels
        self.property_names = [
            "date",
            "coherence_pathway",
            "chemical",
            "concentration",
            "postproc_type",
        ]
        self.property_labels = [
            "Date",
            "Coherence Pathway",
            "Chemical Name",
            "Concentration",
            "Postproc Type",
        ]

    def __enter__(self):
        # Open HDF5 file and store reference
        self.hdf_file = h5py.File(self.hdf_filename, "r+")
        if self.nodename in self.hdf_file:
            self.node = self.hdf_file[self.nodename]
        else:
            raise ValueError(
                f"Expno '{self.nodename}' not found in the HDF5 file."
            )
        self.read_acq_params()
        self.init_ui()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        # Close the HDF5 file when exiting
        if self.hdf_file:
            self.hdf_file.close()

    def read_acq_params(self):
        if (
            "other_info" in self.node
            and "acq_params" in self.node["other_info"]
        ):
            self.acq_params = self.node["other_info"]["acq_params"].attrs
            print("The string attributes will appear here:")
            for k, v in self.acq_params.items():
                # Decode byte strings to regular strings, if necessary
                if isinstance(v, bytes):
                    v = v.decode("utf-8")
                print(k, v)
        else:
            raise ValueError(
                f"acq_params not found in '{self.nodename}/other_info'."
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

            # Retrieve the existing v from the HDF5 acq_params and use it to
            # pre-fill the text fields
            v = self.acq_params.get(prop_name, "")
            if isinstance(v, bytes):
                v = v.decode("utf-8")
            print("attempting", prop_name, "which is", v)
            text_field = QLineEdit(str(v))

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
            acq_params = self.node["other_info"]["acq_params"].attrs

            # Update the HDF5 file with new values from the text fields
            for prop_name in self.property_names:
                value = self.input_fields[prop_name].text()
                # Convert string values to byte strings if they are not already
                if isinstance(value, str):
                    value = value.encode("utf-8")
                acq_params[prop_name] = value

            QMessageBox.information(
                self, "Success", "Values saved successfully!"
            )
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))


def main():
    if len(sys.argv) != 4:
        print(
            "Usage: python hack_acq_params.py <nodename> <filename_pattern>"
            " <exp_type>"
        )
        sys.exit(1)

    hdf_filename = search_filename(
        sys.argv[2], unique=True, exp_type=sys.argv[3]
    )
    nodename = sys.argv[1]

    app = QApplication(sys.argv)
    with EditAcqParams(hdf_filename, nodename) as editor:
        print(editor)
        pass
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
