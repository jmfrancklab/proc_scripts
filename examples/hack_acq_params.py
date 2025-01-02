"""
Hack H5 Properties
==================

There are various labeling properties that are stored in the HDF5 file.
Sometimes we mess up, and we need to hack these properties ex-post-facto.

Note that is intended *only* for hacking labeling/categorizing properties, and
should *never* be used to hack information about pulse program parameters, etc,
that were actually sent to the spectrometer.

The GPT conversation is `here
<https://chatgpt.com/c/66f5aef7-8f68-800b-8a13-0b133ede1ca5>`_.
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


class NodeAsDict:
    def __init__(self, node):
        self.node = node

    def __getitem__(self, key):
        if key in self.node.attrs:
            value = self.node.attrs[key]
            if isinstance(value, bytes):
                return value.decode(
                    "utf-8"
                )  # Decode byte string to UTF-8 string
            return value
        elif key in self.node:
            return NodeAsDict(
                self.node[key]
            )  # Return a group as another NodeAsDict instance
        else:
            return None

    def __setitem__(self, key, value):
        if isinstance(value, (int, float)):  # Numbers stored as attributes
            self.node.attrs[key] = value
        elif isinstance(
            value, str
        ):  # Strings stored as byte strings (UTF-8 encoded)
            self.node.attrs[key] = value.encode("utf-8")
        elif isinstance(value, dict):  # Dictionaries stored as groups
            if key not in self.node:
                self.node.create_group(key)  # Create group if it doesn't exist
            group = self.node[key]
            node_dict = NodeAsDict(group)
            for subkey, subvalue in value.items():
                node_dict[
                    subkey
                ] = subvalue  # Recursively assign subkeys and values
        else:
            raise TypeError(
                f"Unsupported value type for key '{key}': {type(value)}"
            )

    def __delitem__(self, key):
        if key in self.node.attrs:
            del self.node.attrs[key]
        elif key in self.node:
            del self.node[key]
        else:
            raise KeyError(f"Key '{key}' not found")

    def __contains__(self, key):
        return key in self.node.attrs or key in self.node

    def get(self, key, default=None):
        """Get a value from the attributes or group, with UTF-8 decoding and
        default handling."""
        if key in self.node.attrs:
            value = self.node.attrs[key]
            if isinstance(value, bytes):
                return value.decode(
                    "utf-8"
                )  # Decode byte string to UTF-8 string
            return value
        elif key in self.node:
            return NodeAsDict(
                self.node[key]
            )  # Return a group as another NodeAsDict instance
        else:
            return default

    def keys(self):
        """Yield sorted attribute keys and group names together
        alphabetically."""
        all_keys = list(self.node.attrs.keys()) + list(self.node.keys())
        for key in sorted(all_keys):
            yield key

    def items(self):
        """Yield key-value pairs using __getitem__ for both attributes and
        groups."""
        for key in self.keys():
            yield key, self[key]

    def __str__(self):
        """Return a string representation similar to a dictionary."""
        return (
            "{"
            + ", ".join(f"{repr(k)}: {repr(v)}" for k, v in self.items())
            + "}"
        )

    def __repr__(self):
        return f"<NodeAsDict for HDF5 node '{self.node.name}'>"


class EditAcqParams(QWidget):
    def __init__(self, hdf_filename, nodename):
        super().__init__()
        self.hdf_filename = hdf_filename
        self.nodename = nodename
        self.hdf_file = None
        self.node = None
        self.acq_params = None

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
            # Set acq_params as a NodeAsDict object
            self.acq_params = NodeAsDict(self.node["other_info"]["acq_params"])
            print("The string attributes will appear here:")
            for k, v in self.acq_params.items():
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

            # Retrieve the existing value from NodeAsDict with default handling
            v = self.acq_params.get(prop_name, "")
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
        self.hdf_file = h5py.File(self.hdf_filename, "r+")
        if self.nodename in self.hdf_file:
            self.node = self.hdf_file[self.nodename]
        else:
            raise ValueError(
                f"Expno '{self.nodename}' not found in the HDF5 file."
            )
        self.read_acq_params()
        try:
            # Update the HDF5 file with new values from the text fields
            for prop_name in self.property_names:
                value = self.input_fields[prop_name].text()
                self.acq_params[
                    prop_name
                ] = value  # Set the value using NodeAsDict

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
