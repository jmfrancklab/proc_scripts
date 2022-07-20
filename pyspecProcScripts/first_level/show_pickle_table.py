import pickle
def show_pickle_table(pickle_file):
    """show a pickle file with scalar variables in a
    latex table"""
    with open(pickle_file, "rb") as fp:
        variables = pickle.load(fp)
    print(r'\par\begin{tabular}{cc}')
    print(r'\textbf{label} & \textbf{concentration}\\\hline\hline')
    for k, v in variables.items():
        print(f"{k} & {v} \\\\")
    print(r'\end{tabular}\par')
