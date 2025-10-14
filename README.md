# snakemake pipeline to assemble yeast and filamentous fungi using PacBio HiFi technology

## Requirements
	- PacBio HiFi reads
	- Illumina reads
	- Reference genomes (ie. contigs or chromosome-scale NCBI, Ensembl or JGI-Mycocosm)
	- snakemake installed (conda)

## Usage
	- snakemake --use-conda --cores 32 --printshellcmds --latency-wait 30 --keep-going --rerun-incomplete --configfile config/config.yaml
