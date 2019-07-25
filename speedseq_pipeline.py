#!/usr/bin/env python
"""

    speedseq_pipeline.py
                        [--settings PATH] (by default ./settings.cfg)
                        [--log_file PATH]
                        [--verbose]
                        [--target_tasks]  (by default the last task in the pipeline)
                        [--jobs N]        (by default 1)
                        [--just_print]
                        [--flowchart]
                        [--key_legend_in_graph]
                        [--forced_tasks]

"""
import sys
import os
import glob
import time


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Configuration / options

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


if __name__ == '__main__':
    
    from ruffus.proxy_logger import *
    from ruffus import *
    #import drmaa
    
    
    import pipeline.config
   
    # parse command-line args, and check for mandatory options
    options, helpstr = pipeline.config.parse_cli_args()
    pipeline.config.check_mandatory_options(options, 
                                        mandatory_options = [], 
                                        help_text=helpstr)
    
    # configure logging
    logger = pipeline.config.setup_logging(options.log_file, options.verbose)
    # allow logging across Ruffus pipeline
    def get_logger(logger_name, args):
        return logger

    (logger_proxy, logging_mutex) = \
        make_shared_logger_and_proxy(get_logger, logger.name, {})


    # Get pipeline settings from a config file  
    from pipeline.config import PipelineConfig
    cfg = PipelineConfig.getInstance() 
    cfg.set_logger(logger)
    cfg.set_num_jobs(options.jobs)
    cfg.load_settings_from_file(options.pipeline_settings if options.pipeline_settings != None 
                                                        else "settings.cfg")

    # init drmaa
    #cfg.drmaa_session = drmaa.Session()
    #cfg.drmaa_session.initialize()




#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#  Pipeline

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


from pipeline.utils import run_cmd, run_piped_command
from pipeline.tasks import *

from ruffus import *


#
# Prepare directory for every sample and link the input fastq files
# Expected format:
#    /path/to/file/[SAMPLE_ID]_R[12].fastq.gz
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"
#
@jobs_limit(1)    # to avoid problems with simultanous creation of the same sample dir
@subdivide(cfg.input_fastqs,
           formatter('(?P<PATH>.+)/(?P<SAMPLE_ID>[^/]+)_R[12].fastq\.gz$'), 
           os.path.join(cfg.runs_scratch_dir, '{SAMPLE_ID[0]}/{basename[0]}{ext[0]}'))
def link_fastqs(fastq_in, fastq_out):
    """Make working directory for every sample and link fastq files in"""
    if not os.path.exists(os.path.dirname(fastq_out)):
        os.mkdir(os.path.dirname(fastq_out))
    if not os.path.exists(fastq_out):
        os.symlink(fastq_in, fastq_out) 


#####################3
#
# Mapping
#
################

 
MAPPING_JOB_LIMIT=3

@jobs_limit(MAPPING_JOB_LIMIT)
@collate(link_fastqs, 
        formatter("(.+)/(?P<SAMPLE_ID>[^/]+)_R[12].fastq\.gz$"),
        ("{subpath[0][0]}/{SAMPLE_ID[0]}.bam", 
         "{subpath[0][0]}/{SAMPLE_ID[0]}.splitters.bam", 
         "{subpath[0][0]}/{SAMPLE_ID[0]}.discordants.bam"
        ), "{SAMPLE_ID[0]}")
#        regex(r'(.+)/([^/]+)_R[12].fastq\.gz$'),  
#        (r'\1/\2.bam', r'\1/\2.splitters.bam', r'\1/\2.discordants.bam'))
def align_reads(fastqs, bams, sample_id):
    read_group = "@RG\\tID:{id}\\tSM:{sm}\\tLB:{lb}\\tPL:{pl}\
                 ".format(id=sample_id, sm=sample_id, lb=sample_id, pl="ILLUMINA")
    bam_prefix = bams[0][:-len('.bam')]
    speedseq_align(bam_prefix, read_group, cfg.reference, fastqs[0], fastqs[1], threads=int(cfg.num_jobs/MAPPING_JOB_LIMIT))


####################
# 
# Calling
#
################


@merge(align_reads, os.path.join(cfg.runs_scratch_dir, 'multisample.vcf.gz'))
def call_variants(bams, vcf):
    out_prefix = vcf[:-len(".vcf.gz")]
    concordant_bams = [bam for (bam, _, _) in bams]
    speedseq_var(out_prefix, cfg.reference, concordant_bams, 
                 cfg.speedseq_include_bed, threads=cfg.num_jobs)

@merge(align_reads, os.path.join(cfg.runs_scratch_dir, 'multisample.sv.lumpy.vcf.gz'))
def call_svs(bams, vcf):
    out_prefix = vcf[:-len(".vcf.gz")]
    concordant_bams = [bam for (bam, _ , _ ) in bams]
    splitters_bams  = [bam for ( _ ,bam, _ ) in bams]
    discordant_bams = [bam for ( _ , _ ,bam) in bams]
    
    speedseq_sv(out_prefix, cfg.reference, 
                concordant_bams, splitters_bams, discordant_bams,
                exclude_bed=cfg.speedseq_lumpy_exclude_bed, genotype=False, 
                threads=cfg.num_jobs)

@merge(align_reads, os.path.join(cfg.runs_scratch_dir, 'multisample.sv.smoove.vcf.gz'))
def call_svs_smoove(bams, vcf):
    bams = [bam for (bam, _ , _ ) in bams]
    
    cmd = "docker run --rm \
             -v {workdir}:/work \
             -v {refdir}:/reference \
             -v {workdir}:{workdir} \
             -v {refdir}:{refdir} \
             brentp/smoove smoove".format(workdir=cfg.runs_scratch_dir,
                                    refdir=cfg.reference_root)
                                    
    exclude_bed = cfg.speedseq_lumpy_exclude_bed
    args = "call --genotype -x -d --name multisample.sv \
            --outdir {outdir} --fasta {ref} \
            {exclude} -p {threads} \
            {bams}".format(outdir=cfg.runs_scratch_dir,
                           ref=cfg.reference,
                           exclude="" if exclude_bed is None else "--exclude "+exclude_bed,
                           threads=cfg.num_jobs,
                           bams = ' '.join(bams))
            
    run_cmd(cmd+" {args}", args)

    # rename the output, delete temp output folder

@merge(align_reads, os.path.join(cfg.runs_scratch_dir, 'sv_manta', 'results', 'variants','diploidSV.vcf.gz'))
def call_svs_manta(bams, vcf):
    
    args = "--referenceFasta={ref} --runDir={rundir}\
           ".format(ref=cfg.reference,
                    rundir=os.path.join(cfg.runs_scratch_dir, 'sv_manta'))
    bams = [bam for (bam, _ , _ ) in bams]
    for bam in bams:
        args += " --bam={}".format(bam)
                        
    run_cmd(PipelineConfig.getInstance().manta, args)
    run_cmd(os.path.join(cfg.runs_scratch_dir, "sv_manta", "runWorkflow.py {args}"),
                         "-m local -g 32 -j {threads} --quiet".format(threads=cfg.num_jobs))


#######################
#
# Variant postprocessing
#
####################




@transform(call_variants, suffix(".vcf.gz"), ".norm.vcf.gz")
def normalize_variants(vcf, norm_vcf):
    vt_normalize(vcf, norm_vcf, cfg.reference)


@transform(normalize_variants, suffix(".vcf.gz"), ".innerAC9.vcf.gz")
def filter_inner_AC(normalized_vcf, filtered_vcf):
    AC_cutoff = 9
    tmp_vcf = vep_vcf.split(".")[:-1]
    args = "filter \
            -e 'INFO/AC[0]>{AC_cutoff}' \
            -o {outvcf} {invcf} \
           ".format(AC_cutoff=AC_cutoff,
                    outvcf = tmp_vcf,
                    invcf  = normalized_vcf)
    
    run_cmd(cfg.bcftools, args, None)
    
    
    bgzip_and_tabix(tmp_vcf)
    os.remove(tmp_vcf)



@transform(filter_inner_AC, suffix(".vcf.gz"), ".vep.vcf.gz")
def annotate_filtered_vcf(vcf, vep_vcf):
    vep_annotate_vcf(vcf, vep_vcf, threads = cfg.num_jobs)


@transform(call_svs_manta, formatter(), os.path.join(cfg.runs_scratch_dir, "multisample.sv.manta.vep.vcf.gz"))
def annotate_manta_diploidSVs(vcf, vep_vcf):
    vep_annotate_SVs(vcf, vep_vcf, threads = cfg.num_jobs)

@transform(call_svs_smoove, suffix(".vcf.gz"), ".vep.vcf.gz")
def annotate_smoove_svs(vcf, vep_vcf):
    vep_annotate_SVs(vcf, vep_vcf, threads = cfg.num_jobs)


##################
#
#  QC 
#
#############

   
@follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc')), mkdir(os.path.join(cfg.runs_scratch_dir,'qc','read_qc')))
@transform(link_fastqs, formatter('.+/(?P<SAMPLE_ID>[^/]+)\.fastq\.gz$'), 
           os.path.join(cfg.runs_scratch_dir,'qc','read_qc/')+'{SAMPLE_ID[0]}_fastqc.html')
def qc_raw_fastq(input_fastq, report):
    """ Generate FastQC report for raw FASTQs """
    produce_fastqc_report(input_fastq, os.path.dirname(report))


@follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc')))
@transform(call_variants, formatter(), os.path.join('{subpath[0][0]}','qc','multisample.variant_stats'))
def qc_multisample_vcf(vcf, output):
    """ Generate variant QC table for all samples """    
    args = "stats -F {ref} -s - {vcf} > {out}\
            ".format(ref=cfg.reference, 
                     vcf=vcf,
                     out=output)
    
    run_cmd(cfg.bcftools, args)
    




####################
#
#   Cleanup intermediate files
#
##################


def cleanup_files():
    pass


@posttask(cleanup_files)
@follows(qc_multisample_vcf)
def complete_run():
    pass





#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Main logic


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888
if __name__ == '__main__':
    if options.just_print:
        pipeline_printout(sys.stdout, options.target_tasks, options.forced_tasks,
                            gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                            verbose=options.verbose, #verbose_abbreviated_path=0,
                            checksum_level = 0)

    elif options.flowchart:
        pipeline_printout_graph (   open(options.flowchart, "w"),
                                    # use flowchart file name extension to decide flowchart format
                                    #   e.g. svg, jpg etc.
                                    os.path.splitext(options.flowchart)[1][1:],
                                    options.target_tasks,
                                    options.forced_tasks,
                                        gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                                    no_key_legend   = not options.key_legend_in_graph)
    else:        
        pipeline_run(options.target_tasks, options.forced_tasks,
                            multithread     = options.jobs,
                            logger          = cfg.logger,
                            verbose         = options.verbose,
                            gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                            checksum_level  = 0)
    
        
    #cfg.drmaa_session.exit()
    
