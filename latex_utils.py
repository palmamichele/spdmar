
import pandas as pd
import sys 
import os 

def csv_to_latex_table(
    csv_file,
    output_file="table.tex",
    sample_sizes=(200, 500, 1000),
    methods=("RM", "PGD", "VAR"),
    include_covg=False
):
    # Read CSV
    df = pd.read_csv(csv_file, sep=None, engine="python")

    # First two columns: index and q_ij
    index_col = df.columns[0]
    q_col = df.columns[1]

    # Number of columns per method
    metrics = ["hatq_ij", "RMSE"]
    if include_covg:
        metrics.append("COVG")

    n_metrics = len(metrics)
    n_methods = len(methods)
    n_samples = len(sample_sizes)

    # Build LaTeX column specification
    col_spec = "r S[table-format=-1.3]"

    for _ in range(n_samples * n_methods * n_metrics):
        col_spec += " S[table-format=-3.3]"

    lines = []

    lines.append(r"\begin{table}")
    lines.append(r"\centering")
    lines.append(r"\footnotesize")
    lines.append(r"\setlength{\tabcolsep}{3pt}")
    lines.append(r"\resizebox{\linewidth}{!}{%")
    lines.append(r"\begin{tabular}{" + col_spec + "}")
    lines.append(r"\toprule")

    # First header: sample sizes
    header1 = " & & "
    header1 += " & ".join(
        [rf"\multicolumn{{{n_methods * n_metrics}}}{{c}}{{{n}}}"
         for n in sample_sizes]
    )
    lines.append(header1 + r" \\")

    # Rules below sample sizes
    for i in range(n_samples):
        start = 3 + i * n_methods * n_metrics
        end = start + n_methods * n_metrics - 1
        lines.append(rf"\cmidrule(lr){{{start}-{end}}}")

    # Second header: methods
    header2 = " & & "
    header2 += " & ".join(
        [rf"\multicolumn{{{n_metrics}}}{{c}}{{{method}}}"
         for _ in sample_sizes
         for method in methods]
    )
    lines.append(header2 + r" \\")

    # Rules below methods
    for i in range(n_samples * n_methods):
        start = 3 + i * n_metrics
        end = start + n_metrics - 1
        lines.append(rf"\cmidrule(lr){{{start}-{end}}}")

    # Third header: metrics
    metric_labels = {
        "hatq_ij": r"$\hat q_{ij}$",
        "RMSE": "RMSE",
        "COVG": "COVG"
    }

    header3 = " & {$q_{ij}$}"
    for _ in sample_sizes:
        for _ in methods:
            for metric in metrics:
                header3 += " & {" + metric_labels[metric] + "}"

    lines.append(header3 + r" \\")
    lines.append(r"\midrule")

    # Data rows
    for _, row in df.iterrows():
        values = [
            str(int(row[index_col])),
            f"{row[q_col]:.3f}"
        ]

        for n in sample_sizes:
            for method in methods:
                for metric in metrics:
                    col = f"{n} {metric} {method}"
                    values.append(f"{row[col]:.3f}")

        lines.append(" & ".join(values) + r" \\")

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}}")
    lines.append(r"\end{table}")

    # Save LaTeX
    latex = "\n".join(lines)

    with open(output_file, "w", encoding="utf-8") as f:
        f.write(latex)

    return latex




if __name__ == "__main__":

    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print(
            "Usage: python latex_utils.py "
            "<inputPathToCSV> [<outputPathTex OR .>]"
        )
        sys.exit(1)

    input_path = sys.argv[1]

    output_path = sys.argv[2] if len(sys.argv) == 3 else "table.tex"

    if output_path == ".":
        output_path = "table.tex"

    use_covg = (
        input("Use COVG field? (y/n): ").strip().lower()
        in ("y", "yes", "1")
    )

    try:
        csv_to_latex_table(
            input_path,
            output_file=output_path,
            sample_sizes=(200, 500, 1000),
            methods=("RM", "PGD", "VAR"),
            include_covg=use_covg
        )

        print(f"LaTeX table successfully saved to: {output_path}")

    except (FileNotFoundError, ValueError, KeyError) as e:
        print(f"Error: {e}")
        sys.exit(1)
