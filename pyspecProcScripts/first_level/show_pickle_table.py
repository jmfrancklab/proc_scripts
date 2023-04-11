import pickle
import os


def setup_pickle_table(pickle_file, whatprops):
    # {{{ make a pickle file, all_concs.pickle that has a
    #     dictionary of the concentrations of all samples
    #     that we are looking at
    #     the label should have the batch date and
    #     residue, so should be good identifier
    if os.path.exists(pickle_file):
        with open(pickle_file, "rb") as fp:
            pickle_vars = pickle.load(fp)
    else:
        pickle_vars = {}
    # }}}
    # {{{ set up properties we will be documenting
    for thisvar in whatprops:
        if thisvar not in pickle_vars.keys():
            pickle_vars[thisvar] = {}
    # }}}
    return pickle_vars
def show_pickle_table(pickle_file):
    """show a pickle file with scalar variables in a
    latex table"""
    with open(pickle_file, "rb") as fp:
        variables = pickle.load(fp)
    all_props = []
    all_keys = set()
    for k, v in variables.items():
        if type(v) is dict:
            all_props.append(k)
            all_keys |= set(v.keys())
        else:
            raise ValueError(
                f"pickle file {os.getcwd()} {pickle_file} is now expected to contain only dictionaries, which are expected to provide different properties (dict name) pertaining to different samples (dict keys), with the samples shared between dicts.  This doesn't appear to be in that format.  If you have pre-9/1/22 pickle files, you may need to delete them and start over"
            )
    all_keys = list(all_keys)
    all_keys.sort()
    print(r"\par\begin{tabular}{cc}")
    print(r"\textbf{label} & ", end="")
    print(" & ".join([r"\textbf{" + j + r"}" for j in all_props]))
    print(r"\\\hline\hline")
    for k in all_keys:
        print(f"{k}", end="")
        for thisprop in all_props:
            if k in variables[thisprop].keys():
                print(f" & {variables[thisprop][k]}", end="")
            else:
                print(f" & ", end="")
        print("\\\\")
    print(r"\end{tabular}\par")
