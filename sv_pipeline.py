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
from pipeline.tasks import speedseq_sv, cnvnator_sv, cnvnator_calls2bed

from ruffus import *


@transform(cfg.input_bams, suffix(".bam"), ".bam.sv.cnvnator.calls")
def call_sv_cnvnator(inbam, calls):
    """ call CNVs using CNVnator """
    cnvnator_sv(inbam, calls, cfg.reference_chr_dir)

@transform(call_sv_cnvnator, suffix(".bam.sv.cnvnator.calls"), add_inputs(r"\1.bam"), ".bam.sv.cnvnator.bed")
def convert_cnvnator2bed(inputs, bed):
    """ convert CNVnator call to BED format """
    calls, bam = inputs
    cnvnator_calls2bed(calls, bed, bam)


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
    
