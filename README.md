
# Whole-Genome Sequencing Genotyping Pipeline 



## Pipeline

The pipeline consists of tools such as:
 - Trimmomatic 
 - BWA 
 - SAMtools 
 - Freebayes



## Setup

### Pipeline

1. Install Python if you don't have it from before, and a cool python package - `Ruffus` (http://www.ruffus.org.uk/). 
Running jobs on a cluster (PBS, Slurm, etc) requires `drmaa` package. 
You might also need following packages: optparse, logging, shutil

2. Clone the pipeline repository

The pipeline is ready now, but you will need all of its components to perform the analysis.

### Pipeline components

Get a Docker image ....

OR

Install following tools manually:
2. Trimmomatic (http://www.usadellab.org/cms/?page=trimmomatic)
3. BWA (http://sourceforge.net/projects/bio-bwa/files)
4. SAMtools (https://github.com/samtools/samtools)
5. Freebayes

Optional tool for some QC tasks:
8. QualiMap (http://qualimap.bioinfo.cipf.es)

### Reference data

1. Download reference genome, for instance from GATK resource bundle:
https://www.broadinstitute.org/gatk/download

2. The reference genome will also require creating an index if you are planning to run alignment:
`bwa index -a bwtsw human_g1k_v37.fasta` 
 

## Usage

The NGS pipeline is run using `genotyping_pipeline.py` script:

* Running the script

    You can run the script using `python <PIPELINE_HOME>/genotyping_pipeline.py`.
    A list of possible options will be presented. The only required option is `--run_folder RUN_FOLDER`, 
    which specifies the location of the run folder dumped by Illumina seqeuncer.
    
    Important part of the pipeline is the settings file which contains paths to resources 
    (e.g. adapters, reference genome, annotation databases), docker settings, 
    and docker containers to use. See an exemplary file for all required options 
    in <PIPELINE_HOME>/pipeline_settings.cfg.
    If the settings file is not given as argument (--settings), it is expected in the RUN_FOLDER/settings.cfg
  
    If you want to follow progress of the script, use the verbose option (`-vvvvv`).
    In order to use multithreading, use the `-j` option (`-j 12`).

* Outputs

    The script will create RUN_ID folder in the scratch-root directory (given in settings). 
    Inside there will be several directories: 
    	SAMPLE_ID/ - one dir per sample, named after samples found in the RUN_FOLDER 
    	fastqs/    - fastq files
    	drmaa/     - slurm scripts created automatically (for debugging purposes)
    	qc/        - qc output

    After finishing, the sample directories will contain:
    	- a BAM file (`SAMPLE_ID.bam`)
    	- VCF file (`SAMPLE_ID.vcf`) with variants 
    	- quality control files (e.g. coverage statistics)
   
    Directly in the RUN_ID directory additional files will be created:
    	- multisample VCF with snps and indels for all samples (`RUN_ID.multisample.vcf`)

* Typical usage

    For running the genotyping analysis using 12 concurrent tasks:

	genotyping_pipeline.py --run_folder /incoming/RUN_XXXX_YYYY_ZZZZZ \
						    --settings /my_dir/pipeline_settings/my_settings.cfg \
							--target complete_run \
							-vvv -j 12 \
							--log_file /my_dir/pipeline_run_XXX.log





