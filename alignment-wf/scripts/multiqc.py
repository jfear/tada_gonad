from snakemake.shell import shell

def main():
    analysis_directory = "../output/alignment-wf"
    outdir = "../output/alignment-wf" 
    basename = "multiqc"
    shell(
        'LC_ALL=en_US.UTF.8 LC_LANG=en_US.UTF-8 '
        'multiqc '
        '--quiet '
        '--outdir {outdir} '
        '--force '
        '--filename {basename} '
        '--config {snakemake.input.config} '
        '{analysis_directory} '
        '&> {snakemake.log} '
    )

if __name__ == "__main__":
    main()