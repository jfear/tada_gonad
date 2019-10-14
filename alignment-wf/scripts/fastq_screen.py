import os
from snakemake.shell import shell
import sys
import tempfile

# Pull in parameters
extra = snakemake.params.get("extra", "")
aligner = snakemake.params.get("aligner", "bowtie2")
subset = snakemake.params.get("subset", 100000)

# Make temporary DATABASE file
tmp = tempfile.NamedTemporaryFile(delete=False).name
with open(tmp, "w") as fout:
    for k, v in snakemake.input.items():
        if k != "fastq":
            label = k
            index = v.replace(".1.bt2", "")
            fout.write("\t".join(["DATABASE", label, index, aligner.upper()]) + "\n")
    config_file = tmp

# fastq_screen hard-codes filenames according to this prefix. We will send
# hard-coded output to a temp dir, and then move them later.
tempdir = tempfile.mkdtemp()

# Note that we assume only R1 is coming in.
prefix = os.path.basename(snakemake.input.fastq.split(".fastq")[0])

shell(
    "fastq_screen --outdir {tempdir} "
    "--force "
    "--aligner {aligner} "
    "--conf {config_file} "
    "--subset {subset} "
    "--threads {snakemake.threads} "
    "{extra} "
    "{snakemake.input.fastq} "
    "> {snakemake.log} 2>&1"
)

# Move output to the filenames specified by the rule
shell("cp {tempdir}/{prefix}_screen.txt {snakemake.output.txt}")

# Check for the output of the --tag option to fastq_screen
if os.path.isfile("{tempdir}/{prefix}.tagged.fastq.gz"):
    shell("cp {tempdir}/{prefix}.tagged.fastq.gz {snakemake.output.txt}.tagged.fastq.gz")

# Check for the output of the --filter XXXXXX option to fastq_screen
if os.path.isfile("{tempdir}/{prefix}.tagged_filter.fastq.gz"):
    shell(
        "cp {tempdir}/{prefix}.tagged_filter.fastq.gz {snakemake.output.txt}.tagged_filter.fastq.gz"
    )

# Clean up temp
shell("rm -r {tempdir}")
shell("rm {tmp}")
