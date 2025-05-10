import os
import pandas as pd

quast_dirs = ["spades_noerror","spades_error","flye_noerror","flye_error"]

# Metrics of interest and their keys in QUAST's report.tsv
metrics = {
    "Assembly": "Assembly",
    "Total length": "Total length (>= 0 bp)",
    "# contigs": "# contigs (>= 0 bp)",
    "GC (%)": "GC (%)",
    "N50": "N50",
    "N90": "N90",
    "L50": "L50",
    "Largest contig": "Largest contig",
    "Genome fraction (%)": "Genome fraction (%)",
    "Duplication ratio": "Duplication ratio",
    "# misassemblies": "# misassemblies",
    "# mismatches per 100 kbp": "# mismatches per 100 kbp",
    "# indels per 100 kbp": "# indels per 100 kbp"
}

summary_rows = []

for subdir in quast_dirs:
    report_path = os.path.join(subdir, "report.tsv")
    if os.path.exists(report_path):
        df = pd.read_csv(report_path, sep='\t', header=None, names=["metric", "value"])
        df = df.set_index("metric")["value"]
        row = {"Assembly": subdir}
        for key, metric in metrics.items():
            row[key] = df.get(metric, "NA")
        summary_rows.append(row)

# Scan all folders in the current directory
#for subdir in sorted(os.listdir(".")):
#    if os.path.isdir(subdir):
#        report_path = os.path.join(subdir, "report.tsv")
#        if os.path.exists(report_path):
#            df = pd.read_csv(report_path, sep='\t', header=None, names=["metric", "value"])
#            df = df.set_index("metric")["value"]
#            row = {"Assembly": subdir}
#            for key, metric in metrics.items():
#                row[key] = df.get(metric, "NA")
#            summary_rows.append(row)

# Write summary to TSV
summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv("quast_summary_spades_flye.tsv", sep="\t", index=False)
print("[✓] QUAST summary saved to quast_summary_spades_flye.tsv")
