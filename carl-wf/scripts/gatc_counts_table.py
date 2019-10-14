import re

import pandas as pd

HEADER = [
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "count",
    "covered",
    "bases",
    "fraction",
]


def main():
    df = pd.concat((parse_counts_file(file_name) for file_name in snakemake.input), axis=1)
    df.rename_axis("GATC_site").to_csv(snakemake.output[0], sep="\t")


def parse_counts_file(file_name):
    sample = re.findall(r"carl-wf/(?P<sample>.*)/", file_name)[0]
    return (
        pd.read_csv(file_name, sep="\t", header=None, names=HEADER, usecols=["name", "count"])
        .set_index("name")
        .loc[:, "count"]
        .rename(sample)
    )


if __name__ == "__main__":
    main()
