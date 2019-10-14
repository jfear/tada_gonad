"""Map GATC to FBgn.

Use bedtools to find which genes each GATC overlapped. Make a simple table of
gatc_id to FBgn.

"""
import re

from pybedtools import BedTool


def main():
    gatc = BedTool(snakemake.input.gatc)
    gtf = BedTool(snakemake.input.gtf)

    with open(snakemake.output[0], "w") as fo:
        fo.write("gatc_id\tFBgn\n")
        for record in gatc.intersect(gtf, wo=True):
            if record[8] == "gene":
                fo.write("\t".join(gatc_to_fbgn(record)) + "\n")


def gatc_to_fbgn(record):
    gatc_id = record.name
    fbgn = re.findall(f"(FBgn\d+)", record[14])[0]
    return gatc_id, fbgn


if __name__ == "__main__":
    main()
