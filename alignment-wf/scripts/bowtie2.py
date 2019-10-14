from snakemake.shell import shell

def main():
    prefix = snakemake.input.index.replace(".1.bt2", "")
    sam = snakemake.output[0].replace('.bam', '.sam')

    shell(
        "bowtie2 "
        "-x {prefix} "
        "-U {snakemake.input.fastq} "
        '--no-unal '  # NOTE: suppress unaligned reads
        "--threads {snakemake.threads} "
        "-S {sam} "
        "> {snakemake.log} 2>&1"
    )

    shell(
        "samtools view -Sb {sam} "
        "| samtools sort - -o {snakemake.output[0]} -O BAM "
        "&& rm {sam}"
    )

if __name__ == "__main__":
    main()