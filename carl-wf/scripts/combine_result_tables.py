import re
import pandas as pd

def main():
    with pd.ExcelWriter(snakemake.output[0]) as writer:
        for file_name in snakemake.input.files:
            sheet_name = make_sheet_name(file_name)
            df = pd.read_csv(file_name, sep="\t", index_col=0).query("padj <= 0.05")
            if snakemake.input.get("annot", False):
                df = df.join(read_annot(snakemake.input.annot)).set_index("gene_symbol", append=True)
            df.to_excel(writer, sheet_name=sheet_name)


def make_sheet_name(file_name):
    res = re.findall(r".*_(?P<tissue>.*)_(?P<driver>.*)\.(?P<ext>\d+)\.tsv$", file_name)[0]
    return "{}_{}_{}".format(*res)

def read_annot(file_name):
    return (
        pd.read_csv(file_name, sep="\t", usecols=["gene_symbol", "primary_FBgn"])
        .set_index("primary_FBgn")
        .squeeze()
        .rename_axis("FBgn")
    )


if __name__ == "__main__":
    main()