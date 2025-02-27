import pickle


def show_pickle_table(pickle_file):
    """show a pickle file with scalar variables in a
    latex table

    Useful function when comparing multiple data sets and
    takes advantage of the capability of :func:`~pySpecProcScripts.QESR`
    function to output pickle files.
    A latex table is generated pulling from the pickle file called.
    """
    with open(pickle_file, "rb") as fp:
        variables = pickle.load(fp)
    print(r"\par\begin{tabular}{cc}")
    print(r"\textbf{label} & \textbf{concentration}\\\hline\hline")
    for k, v in variables.items():
        print(f"{k} & {v} \\\\")
    print(r"\end{tabular}\par")
