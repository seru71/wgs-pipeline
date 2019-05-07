#!/usr/bin/env python
"""

    sv_pipeline.py
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
# CNVnator
#
@transform(cfg.input_bams, suffix(".bam"), ".bam.sv.cnvnator.calls")
def call_sv_cnvnator(inbam, calls):
    """ Call CNVs using CNVnator """
    cnvnator_sv(inbam, calls, cfg.reference_chr_dir)

@transform(call_sv_cnvnator, suffix(".bam.sv.cnvnator.calls"), add_inputs(r"\1.bam"), ".bam.sv.cnvnator.calls.mqs")
def calc_mean_mq_for_calls(inputs, mqs):
    """ calculates mean MQ for each CNVnator call"""
    calls, bam = inputs
    get_mean_MQ_for_regions(bam, calls, mqs)

@transform(calc_mean_mq_for_calls, suffix(".calls.mqs"), add_inputs(r"\1.calls"), ".bed")
def convert_cnvnator2bed(inputs, bed):
    """ Convert CNVnator call to BED format """
    mqs, calls = inputs
    cnvnator_calls2bed(calls, bed, mqs)


@merge(convert_cnvnator2bed, os.path.join(cfg.runs_scratch_dir, 'multisample.sv.cnvnator.AC.bed'))
def get_cnvnator_AC(beds, ac_bed):
    """ Calculate how common BED entries are in the population of samples """
    args="multiinter -i %s" % ' '.join(beds)
    run_piped_command(cfg.bedtools, args, None,
                      "cut {args}", "-f1-4 > %s" % ac_bed, None)

@transform(convert_cnvnator2bed, suffix(".bed"), add_inputs(get_cnvnator_AC), ".AC.bed")
def add_AC_to_cnvnator_bed(inputs, outbed):
    """ Add AC columns to CNVnator BED """
    bed, ac_bed = inputs
    run_cmd(cfg.bedtools, "map -a {bed} -b {ac} -c 4 -o mean,min,max > {out}\
                          ".format(bed=bed, ac=ac_bed, out=outbed))


@transform(add_AC_to_cnvnator_bed, suffix(".bed"), ".bed.vep")
def annotate_cnvnator_bed(cnvnator_bed, vep_table):
    """ Annotate CNVnator BED using VEP """
    vep_annotate_bed(cnvnator_bed, vep_table)
    
    
@transform(annotate_cnvnator_bed, suffix(".bed.vep"), add_inputs(r"\1.bed"), ".vep.tsv")
def join_cnvnator_annotations_and_calls(inputs, out_table):
    """ Concatenate VEP output with variant call details """
    vep_table, cnvnator_bed = inputs
    with open(vep_table, 'r') as vep, \
         open(cnvnator_bed, 'r') as bed, \
         open(out_table, 'w') as table:
        
        last_vep_record_id = None
        for vl in vep:
            if vl.startswith("##"):
                table.write(vl)
            elif vl.startswith("#Uploaded_variation"):
                table.write('\t'.join(vl.strip().split('\t')[1:] + ["Length","Depth_Ratio","PV1","PV2","PV3","PV4","Q0","MQ_avg","meanAC","minAC","maxAC\n"]))
            else:
                vep_records = vl.strip().split('\t')
                if vep_records[0] != last_vep_record_id:
                    last_vep_record_id = vep_records[0]
                    bed_records = bed.readline().split('\t')
                
                try:
                    table.write('\t'.join(vep_records[1:] + [bed_records[4]] + bed_records[6:]))
                except IndexError:
                    print vep_records
                    print bed_records
    

@merge(cfg.input_bams, os.path.join(cfg.runs_scratch_dir, 'multisample.sv.speedseq.vcf.gz'))
def call_sv_speedseq(bams, vcf):
    out_prefix = vcf[:-len(".vcf.gz")]
    concordant_bams = [bam for (bam, _ , _ ) in bams]
    splitters_bams  = [bam for ( _ ,bam, _ ) in bams]
    discordant_bams = [bam for ( _ , _ ,bam) in bams]
    
    speedseq_sv(out_prefix, cfg.reference, 
                concordant_bams, splitters_bams, discordant_bams,
                exclude_bed=cfg.speedseq_lumpy_exclude_bed, threads=cfg.num_jobs)



##################
#
#  QC variants
#
#############






####################
#
#   Cleanup intermediate files
#
##################


def cleanup_files():
    pass


@posttask(cleanup_files)
@follows(call_sv_cnvnator)
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
    
