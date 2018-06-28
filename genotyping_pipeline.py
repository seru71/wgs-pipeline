#!/usr/bin/env python
"""

    genotyping_pipeline.py
                        --run_folder PATH
                        [--settings PATH] (by default RUN_FOLDER/settings.cfg)
                        [--log_file PATH]
                        [--verbose]
                        [--target_tasks]  (by default the last task in the pipeline)
                        [--jobs N]        (by default 1)
                        [--just_print]
                        [--flowchart]
                        [--key_legend_in_graph]
                        [--forced_tasks]
                        [--run_on_bcl_tile TILE_REGEX]

"""
import sys
import os
import glob


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   options


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


if __name__ == '__main__':
    from optparse import OptionParser
    import StringIO

    parser = OptionParser(version="%prog 1.0", usage = "\n\n    %prog --run_folder PATH_TO_RUN_FOLDER [--settings pipeline_settings.cfg] [--target_task TASK] [more_options]")
    
    parser.add_option("-r", "--run_folder", dest="run_folder",
                        metavar="FILE",
                        type="string",
                        help="Path to the input run folder.")                  
    
                                
    #
    #   general options: verbosity / logging
    #
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="count", 
                      help="Print more verbose messages for each additional verbose level.")
    parser.add_option("-L", "--log_file", dest="log_file",
                      metavar="FILE",
                      type="string",
                      help="Name and path of log file")


    #
    #   pipeline
    #
    parser.add_option("-s", "--settings", dest="pipeline_settings",
                        metavar="FILE",
                        type="string",
                        help="File containing all the settings for the analysis. By default settings.cfg in the run_folder.")                  
    parser.add_option("-t", "--target_tasks", dest="target_tasks",
                        action="append",
                        metavar="JOBNAME",
                        type="string",
                        help="Target task(s) of pipeline.")
    parser.add_option("-j", "--jobs", dest="jobs",
                        metavar="N",
                        type="int",
                        help="Allow N jobs (commands) to run simultaneously.")
    parser.add_option("-n", "--just_print", dest="just_print",
                        action="store_true", 
                        help="Don't actually run any commands; just print the pipeline.")
    parser.add_option("--flowchart", dest="flowchart",
                        metavar="FILE",
                        type="string",
                        help="Don't actually run any commands; just print the pipeline "
                             "as a flowchart.")

    #
    #   Less common pipeline options
    #
    parser.add_option("--key_legend_in_graph", dest="key_legend_in_graph",
                        action="store_true",
                        help="Print out legend and key for dependency graph.")
    parser.add_option("--forced_tasks", dest="forced_tasks",
                        action="append",
                        metavar="JOBNAME",
                        type="string",
                        help="Pipeline task(s) which will be included even if they are up to date.")
    parser.add_option("--rebuild_mode", dest="rebuild_mode",
                        action="store_false", 
                        help="gnu_make_maximal_rebuild_mode")
    parser.add_option("--run_on_bcl_tile", dest="run_on_bcl_tile",
                        type="string",                        
                        help="Use only this tile when doing bcl2fastq conversion. For testing purposes.")
    
    parser.set_defaults(pipeline_settings=None, 
                        jobs=1, verbose=0, 
                        target_tasks=list(), forced_tasks=list(), 
                        just_print=False, key_legend_in_graph=False,
                        rebuild_mode=True, run_on_bcl_tile=None)
    

    # get help string
    f =StringIO.StringIO()
    parser.print_help(f)
    helpstr = f.getvalue()
    (options, remaining_args) = parser.parse_args()


    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Change this if necessary                  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    #
    #   Add names of mandatory options,
    #       strings corresponding to the "dest" parameter
    #       in the options defined above
    #
    mandatory_options = ['run_folder']

    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Change this if necessary                  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


    def check_mandatory_options (options, mandatory_options, helpstr):
        """
        Check if specified mandatory options have been defined
        """
        missing_options = []
        for o in mandatory_options:
            if not getattr(options, o):
                missing_options.append("--" + o)
    
        if not len(missing_options):
            return
    
        raise Exception("Missing mandatory parameter%s: %s.\n\n%s\n\n" %
                        ("s" if len(missing_options) > 1 else "",
                         ", ".join(missing_options),
                         helpstr))
        
    check_mandatory_options(options, mandatory_options, helpstr)
    
    
    #
    # check presence of the run folder, and sample sheet file
    #
    if not os.path.exists(options.run_folder) or not os.path.exists(os.path.join(options.run_folder,'SampleSheet.csv')):
        raise Exception("Missing sample sheet file: %s.\n" % os.path.join(options.run_folder,'SampleSheet.csv'))
            
    


    #vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
    #                                             #
    #   Get pipeline settings from a config file  #
    #                                             #
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    import ConfigParser

    config = ConfigParser.ConfigParser()
    try:
        if options.pipeline_settings == None:
            options.pipeline_settings = os.path.join(options.run_folder, 'settings.cfg')
        config.read(options.pipeline_settings)
    except FileNotFound:
        raise Exception('Provided settings file [%s] does not exist or cannot be read.' % options.pipeline_settings)


    # Root dirs
    reference_root = config.get('Paths','reference-root')
    
    run_id = os.path.dirname(options.run_folder)    # get ID of the run and use it to create scratch results folder
    runs_scratch_dir = os.path.join(config.get('Paths','scratch-root'), run_id)
      
    # optional results and fastq archive dirs  
    results_archive = None
    try:
        results_archive = config.get('Paths','results-archive')
    except ConfigParser.NoOptionError:
        print 'No results-archive provided. Results will not be archived outside of the current (working) directory.'
    
    fastq_archive = None
    try:
        fastq_archive = config.get('Paths','fastq-archive')
    except ConfigParser.NoOptionError:
        print 'No fastq-archive provided. Fastq files will not be archived outside of the current (working) directory.'

    
    # optional /tmp dir
    tmp_dir = '/tmp'
    try:
        tmp_dir = config.get('Paths','tmp-dir')
    except ConfigParser.NoOptionError:
        print 'No custom tmp-dir provided. /tmp will be used.'
    
    
    # reference files
    reference = os.path.join(reference_root, config.get('Resources','reference-genome'))
    capture = os.path.join(reference_root, config.get('Resources','capture-regions-bed'))
    capture_qualimap = os.path.join(reference_root, config.get('Resources','capture-regions-bed-for-qualimap'))
    capture_plus = os.path.join(reference_root, config.get('Resources', 'capture-plus-regions-bed'))
    gene_coordinates = os.path.join(reference_root, config.get('Resources', 'gene-coordinates'))
    
    adapters = os.path.join(reference_root, config.get('Resources', 'adapters-fasta'))
    
    # tools
    bcl2fastq = config.get('Tools','bcl2fastq')
    trimmomatic = config.get('Tools', 'trimmomatic') 
    bwa = config.get('Tools','bwa')
    samtools = config.get('Tools','samtools')
    freebayes = config.get('Tools','freebayes')
    qualimap = config.get('Tools','qualimap')


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Logger


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888



if __name__ == '__main__':
    import logging
    import logging.handlers

    MESSAGE = 15
    logging.addLevelName(MESSAGE, "MESSAGE")

    def setup_std_logging (logger, log_file, verbose):
        """
        set up logging using programme options
        """
        class debug_filter(logging.Filter):
            """
            Ignore INFO messages
            """
            def filter(self, record):
                return logging.INFO != record.levelno

        class NullHandler(logging.Handler):
            """
            for when there is no logging
            """
            def emit(self, record):
                pass

        # We are interesting in all messages
        logger.setLevel(logging.DEBUG)
        has_handler = False

        # log to file if that is specified
        if log_file:
            handler = logging.FileHandler(log_file, delay=False)
            handler.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)6s - %(message)s"))
            handler.setLevel(MESSAGE)
            logger.addHandler(handler)
            has_handler = True

        # log to stderr if verbose
        if verbose:
            stderrhandler = logging.StreamHandler(sys.stderr)
            stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
            stderrhandler.setLevel(logging.DEBUG)
            if log_file:
                stderrhandler.addFilter(debug_filter())
            logger.addHandler(stderrhandler)
            has_handler = True

        # no logging
        if not has_handler:
            logger.addHandler(NullHandler())


    #
    #   set up log
    #
    module_name = "ngspipe"
    logger = logging.getLogger(module_name)
    setup_std_logging(logger, options.log_file, options.verbose)

    #
    #   Allow logging across Ruffus pipeline
    #
    def get_logger (logger_name, args):
        return logger

    from ruffus.proxy_logger import *
    (logger_proxy,
     logging_mutex) = make_shared_logger_and_proxy (get_logger,
                                                    module_name,
                                                    {})






#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#  Common functions 


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

import drmaa
drmaa_session = drmaa.Session()
drmaa_session.initialize()

from ruffus.drmaa_wrapper import run_job, error_drmaa_job

"""
cmd is given in a form:
    
    command {args}
    interpreter {interpreter_args} command {atgs}

The strings {args} and {interpreter_args} are replaced with args and interpreter_args values.
Examples of correct commands:
    cmd = "bcl2fastq {args}"
    cmd = "bcl2fastq -p param {args}"
    cmd = "bcl2fastq -p2 param {args} -p2 param2"
    cmd = "java {interpreter_args} -jar myjarfile.jar {args} -extras extra_arg
    cmd = "java -XmX4G {interpreter_args} -jar myjarfile.jar {args} -extras extra_arg
"""
def run_cmd(cmd, args, interpreter_args=None, run_locally=True,
            cpus=1, mem_per_cpu=1024, walltime='24:00:00', 
            retain_job_scripts = True, job_script_dir = os.path.join(runs_scratch_dir, "drmaa")):
    
    full_cmd = cmd.format(args=args, 
                          interpreter_args = interpreter_args if interpreter_args!=None else "")

    stdout, stderr = "", ""
    job_options = "--ntasks=1 \
                   --cpus-per-task={cpus} \
                   --mem-per-cpu={mem} \
                   --time={time} \
                  ".format(cpus=cpus, mem=int(1.2*mem_per_cpu), time=walltime)
    #print full_cmd                   
    try:
        stdout, stderr = run_job(full_cmd.strip(), 
                                 job_other_options=job_options,
                                 run_locally = run_locally, 
                                 retain_job_scripts = retain_job_scripts, job_script_directory = job_script_dir,
                                 logger=logger, working_directory=os.getcwd(),
                                 drmaa_session = drmaa_session)
    except error_drmaa_job as err:
        raise Exception("\n".join(map(str, ["Failed to run:", cmd, err, stdout, stderr])))



""" 
Only default job scheduling params of run_command available when executing via SLURM.
"""
def run_piped_command(*args):
    run_locally=True
    retain_job_scripts = True
    job_script_dir = os.path.join(runs_scratch_dir, "drmaa")	
    cpus=1
    mem_per_cpu=1024
    walltime="24:00:00"
 
    stdout, stderr = "", ""
    job_options = "--ntasks=1 \
                   --cpus-per-task={cpus} \
                   --mem-per-cpu={mem} \
                   --time={time} \
                  ".format(cpus=cpus, mem=int(1.2*mem_per_cpu), time=walltime)
	
    full_cmd = expand_piped_command(*args)
	
    try:
        stdout, stderr = run_job(full_cmd.strip(), 
                                 job_other_options=job_options,
                                 run_locally = run_locally, 
                                 retain_job_scripts = retain_job_scripts, job_script_directory = job_script_dir,
                                 logger=logger, working_directory=os.getcwd(),
                                 drmaa_session = drmaa_session)
    except error_drmaa_job as err:
        raise Exception("\n".join(map(str, ["Failed to run:", full_cmd, err, stdout, stderr])))
	
def expand_piped_command(cmd, cmd_args, interpreter_args=None, *args):
	expanded_cmd = cmd.format(args=cmd_args, interpreter_args = interpreter_args if interpreter_args!=None else "")
	expanded_cmd += (" | "+expand_piped_command(*args)) if len(args) > 0 else ""
	return expanded_cmd


def produce_fastqc_report(fastq_file, output_dir=None):
    args = fastq_file
    args += (' -o '+output_dir) if output_dir != None else ''
    run_cmd(fastqc, args)



                            
def index_bam(bam):
    """Use samtools to create an index for the bam file"""
    run_cmd(samtools, "index %s" % bam)
                          
def bam_quality_score_distribution(bam,qs,pdf):
    """Calculates quality score distribution histograms"""
    run_cmd(picard, "QualityScoreDistribution \
                    CHART_OUTPUT={chart} \
                    OUTPUT={out} \
                    INPUT={bam} \
                    VALIDATION_STRINGENCY=SILENT \
                    ".format(chart=pdf, out=qs, bam=bam),
            interpreter_args="")

def bam_alignment_metrics(bam,metrics):
    """Collects alignment metrics for a bam file"""
    run_cmd(picard, "CollectAlignmentSummaryMetrics \
                    REFERENCE_SEQUENCE={ref} \
                    OUTPUT={out} \
                    INPUT={bam} \
                    VALIDATION_STRINGENCY=SILENT \
                    ".format(ref=reference, out=metrics, bam=bam),
            interpreter_args="")
    
def bam_target_coverage_metrics(input_bam, output):
    """ Calculates and outputs bam coverage statistics """
    run_cmd(gatk, "-R {reference} \
                    -T DepthOfCoverage \
                    -o {output} \
                    -I {input} \
                    -L {capture} \
                    -ct 8 -ct 20 -ct 30 \
                    --omitDepthOutputAtEachBase --omitLocusTable \
                    ".format(reference=reference,
                             output=output,
                             input=input_bam,
                             capture=capture), 
            interpreter_args="-Xmx4g")

def bam_gene_coverage_metrics(input_bam, output):
    """ Calculates and outputs bam coverage statistics """
    run_cmd(gatk, "-R {reference} \
                    -T DepthOfCoverage \
                    -o {output} \
                    -I {input} \
                    -L {capture} \
                    -geneList {genes} \
                    -ct 5 -ct 10 -ct 20 \
                    --omitDepthOutputAtEachBase --omitLocusTable \
                    ".format(reference=reference,
                             output=output,
                             input=input_bam,
                             capture=capture,
                             genes=gene_coordinates), 
            interpreter_args="-Xmx4g")

def qualimap_bam(input_bam, output_dir):
    """ Generates Qualimap bam QC report """
    # create necessary folders first
    #if not os.path.exists('qc'): os.mkdir('qc')
    #if not os.path.exists('qc/qualimap/'): os.mkdir('qc/qualimap')
    if not os.path.exists(output_dir): os.mkdir(output_dir)
    run_cmd(qualimap, "bamqc -bam {bam} \
                        -c -outformat PDF \
                        -gff {target} \
                        -gd HUMAN -os \
                        -outdir {dir} &> {dir}/qualimap.err \
                        ".format(bam=input_bam,
                                target=capture_qualimap,
                                dir=output_dir),
            interpreter_args="")

def get_sample_ids():
    """ Provides meaningful result only after HaplotypeCaller step"""
    files = glob.glob(os.path.join(runs_scratch_dir,'*','*.gvcf'))
    return [ os.path.splitext(os.path.basename(f))[0] for f in files ]

def get_num_files():
    """ Provides meaningful result only after HaplotypeCaller step"""
    return len(get_sample_ids())
    


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#   Pipeline


#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

from ruffus import *


#
#
# Prepare FASTQ
# 

@follows(mkdir(runs_scratch_dir), mkdir(os.path.join(runs_scratch_dir,'fastqs')))
@files(options.run_folder, os.path.join(runs_scratch_dir,'fastqs','completed'))
@posttask(touch_file(os.path.join(runs_scratch_dir,'fastqs','completed')))
def bcl2fastq_conversion(run_directory, completed_flag):
    """ Run bcl2fastq conversion and create fastq files in the run directory"""
    out_dir = os.path.join(runs_scratch_dir,'fastqs')
    interop_dir = os.path.join(out_dir,'InterOp')

    # r, w, d, and p specify numbers of threads to be used for each of the concurrent subtasks of the conversion (see bcl2fastq manual) 
    args = "-R {indir} -o {outdir} --interop-dir={interopdir} -r1 -w1 -d2 -p4 \
            ".format(indir=run_directory, outdir=out_dir, interopdir=interop_dir)
    if options.run_on_bcl_tile != None:
        args += " --tiles %s" % options.run_on_bcl_tile
    run_cmd(bcl2fastq, args, cpus=8, mem_per_cpu=2048)
    


@active_if(fastq_archive != None)
@transform(bcl2fastq_conversion, formatter(".+/(?P<RUN_ID>[^/]+)/fastqs/completed"), str(fastq_archive)+"/{RUN_ID[0]}")
def archive_fastqs(completed_flag, archive_dir):
    """ Archive fastqs """    
    fq_dir = os.path.dirname(completed_flag)

# uncomment if archive should not be overwritten (risk of creating many archives of the same run)
#    if os.path.exists(archive_dir):
#	import time
#	archive_dir += "_archived_"+str(time.strftime("%Y%m%d_%H%M%S"))

    import shutil
    shutil.move(fq_dir, archive_dir)
    os.mkdir(fq_dir)
    for f in glob.glob(os.path.join(archive_dir,"*.fastq.gz")):
        os.symlink(f, os.path.join(fq_dir,os.path.basename(f)))


#
# Prepare directory for every sample and link the input fastq files
# Expected format:
#    /path/to/file/[SAMPLE_ID]_S[1-9]\d?_L\d\d\d_R[12]_001.fastq.gz
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"
#
@jobs_limit(1)    # to avoid problems with simultanous creation of the same sample dir
@follows(archive_fastqs)
@subdivide(os.path.join(runs_scratch_dir,'fastqs','*.fastq.gz'),
           formatter('(?P<PATH>.+)/(?P<SAMPLE_ID>[^/]+)_S[1-9]\d?_L\d\d\d_R[12]_001\.fastq\.gz$'), 
           '{subpath[0][1]}/{SAMPLE_ID[0]}/{basename[0]}{ext[0]}')
def link_fastqs(fastq_in, fastq_out):
    """Make working directory for every sample and link fastq files in"""
    if not os.path.exists(os.path.dirname(fastq_out)):
        os.mkdir(os.path.dirname(fastq_out))
    if not os.path.exists(fastq_out):
        os.symlink(fastq_in, fastq_out) 

    
    
#
# Input FASTQ filenames are expected to have following format:
#    [SAMPLE_ID]_[S_NUM]_[LANE_ID]_[R1|R2]_001.fastq.gz
# In this step, the two FASTQ files matching on the [SAMPLE_ID]_[S_ID]_[LANE_ID] will be trimmed together (R1 and R2). 
# The output will be written to two FASTQ files
#    [SAMPLE_ID]_[LANE_ID].fq1.gz
#    [SAMPLE_ID]_[LANE_ID].fq2.gz
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"
#
@collate(link_fastqs, regex(r'(.+)/([^/]+)_S[1-9]\d?_(L\d\d\d)_R[12]_001\.fastq\.gz$'),  r'\1/\2_\3.fq1.gz')
def trim_reads(inputs, output):
    outfq1 = output
    outfq2 = output.replace('fq1.gz','fq2.gz')
    unpaired = [outfq1.replace('fq1.gz','fq1_unpaired.gz'), outfq2.replace('fq2.gz','fq2_unpaired.gz')]               
    # logfile = output.replace('fq1.gz','trimmomatic.log')
    # -trimlog {log} \
    # log=logfile
    args = "PE -phred33 -threads 1 \
            {in1} {in2} {out1} {unpaired1} {out2} {unpaired2} \
            ILLUMINACLIP:{adapter}:2:30:10 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36".format(in1=inputs[0], in2=inputs[1],
                              out1=outfq1, out2=outfq2,
                              unpaired1=unpaired[0], unpaired2=unpaired[1],
                              adapter=adapters)
    max_mem = 2048
    run_cmd(trimmomatic, args, interpreter_args="-Xmx"+str(max_mem)+"m", cpus=1, mem_per_cpu=max_mem)


#
#
# Align reads and create raw BAM files (one per sample)
# 



def bwa_map_and_sort(output_bam, ref_genome, fq1, fq2=None, read_group=None, threads=1):
	
	bwa_args = "mem -t {threads} {rg} {ref} {fq1} \
	            ".format(threads=threads, 
                        rg="-R '%s'" % read_group if read_group!=None else "", 
                        ref=ref_genome, fq1=fq1)
	if fq2 != None:
		bwa_args += fq2

	samtools_args = "sort -o {out}".format(out=output_bam)

	run_piped_command(bwa, bwa_args, None,
	                  samtools, samtools_args, None)


"""
Could it be used to merge LANE bams for each sample?


def map_reads(fastq_list, ref_genome, output_bam, read_groups=None):
    
    # If no read groups is provided, we could make up default ones based on fq filenames.
    # This would most likely result in unpaired FQs ending up in different read group ids, and in consequence, being genotyped separately
    if read_groups==None:
        s_ids = [ os.path.basename(s[0][:-len('.fq1.gz')] if isinstance(s, tuple) else s[:-len('fq.gz')]) for s in fastq_list ]
        read_groups = [ '@RG\tID:{sid}\tSM:{sid}\tLB:{sid}'.format(sid=s) for s in s_ids]
    
    tmp_bams = [ output_bam+str(i) for i in range(0, len(fastq_list)) ]
    for i in range(0, len(fastq_list)):
		if isinstance(fastq_list[i], tuple):
			bwa_map_and_sort(tmp_bams[i], ref_genome, fastq_list[i][0], fastq_list[i][1], read_groups[i])
		else:
			bwa_map_and_sort(tmp_bams[i], ref_genome, fastq_list[i], read_groups[i])   
    
    merge_bams(output_bam, *tmp_bams)
    
    for f in tmp_bams:
		  os.remove(f)

"""

#
# FASTQ filenames are expected to have following format:
#    [SAMPLE_ID]_[LANE_ID].fq[1|2].gz
# In this step, the fq1 file coming from trim_reads is matched with the fq2 file and mapped together. 
# The output will be written to BAM file:
#    [SAMPLE_ID]_[LANE_ID].bam
#
#@transform(trim_reads, suffix('.fq1.gz'), add_inputs(r'\1.fq2.gz'), '.bam', )
@transform(trim_reads, formatter("(.+)/(?P<SAMPLE_ID>[^/]+)_(?P<LANE_ID>L\d\d\d).fq[12]\.gz$", 
        "{subpath[0][0]}/{SAMPLE_ID[0]}_{LANE_ID[0]}.bam", "{SAMPLE_ID[0]}", "{LANE_ID[0]}")
def align_reads(fastqs, bam, sample_id, lane_id):
    threads = 2
    
    print sample_id, lane_id
    
    # construct read group information from fastq file name (assuming [SAMPLE_ID]_[LANE_ID].fq[1|2].gz format)
    sample_lane = os.path.basename(fastqs[0])[0:-len(".fq1.gz")]
    sample = "_".join(sample_lane.split("_")[0:-1])
    lane = sample_lane.split("_")[-1]
    read_group = "@RG\\tID:{id}\\tSM:{sm}\\tLB:{lb}\\tPL:{pl}\\tPU:{pu} \
                 ".format(id=sample_lane, sm=sample, lb=sample, pl="ILLUMINA", pu=lane)
                 
    #args = "mem -t {threads} -R {rg} {ref} {fq1} {fq2} \
	#    ".format(threads=threads, rg=read_group, ref=reference, 
    #                 fq1=fastqs[0], fq2=fastqs[1])
    #iargs = "samtools view -b -o {bam} -".format(bam=bam)

    #run_cmd(bwa, args, interpreter_args=iargs, cpus=threads, mem_per_cpu=8192/threads)
    
    bwa_map_and_sort(bam, reference, fastqs[0], fastqs[1], read_group=read_group, threads=1)

#
# BAM filenames are expected to have following format:
#    [SAMPLE_ID]_[LANE_ID].bam
# In this step, all BAM files matching on the SAMPLE_ID will be merged into one BAM file:
#    [SAMPLE_ID].bam
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"
#
@jobs_limit(4)    # to avoid filesystem problems 
@collate(align_reads, regex(r"(.+)/([^/]+).+\.bam$"),  r'\1/\2.bam')
def merge_lanes(lane_bams, out_bam):
    args = "MergeSamFiles O={bam} \
            ASSUME_SORTED=false \
            MAX_RECORDS_IN_RAM=2000000 \
            USE_THREADING=true \
            ".format(bam=out_bam)
    # include all bam files as args
    for bam in lane_bams:
        args += " I={bam}".format(bam=bam)
        
    run_cmd(picard, args, interpreter_args="-Xmx8g", cpus=4, mem_per_cpu=2048)
    

def clean_fastqs_and_lane_bams():
    """ Remove the trimmed fastq files, and SAM files. Links to original fastqs are kept """
    for f in glob.glob(os.path.join(runs_scratch_dir,'*','*.fq[12]*.gz')):
        os.remove(f)
    for f in glob.glob(os.path.join(runs_scratch_dir,'*','*_L\d\d\d.bam')):
        os.remove(f)


@posttask(clean_fastqs_and_lane_bams)
@transform(merge_lanes, suffix(".bam"), '.bam.bai')
def index(bam, output):
    """Create raw bam index"""
    index_bam(bam)



#
#
# QC the raw bam files
#

"""
@follows(index)
@transform(merge_lanes, suffix(".bam"), '.quality_score')
def qc_raw_bam_quality_score_distribution(input_bam, output):
    """docstring for metrics1"""
    bam_quality_score_distribution(input_bam, output, output + '.pdf')

@follows(index)
@transform(merge_lanes, suffix(".bam"), '.metrics')
def qc_raw_bam_alignment_metrics(input_bam, output):
    """docstring for metrics1"""
    bam_alignment_metrics(input_bam, output)
   """ 
@follows(index)
@transform(merge_lanes, suffix(".bam"), '.target_coverage.sample_summary', r'\1.target_coverage')
def qc_raw_bam_target_coverage_metrics(input_bam, output, output_format):
    bam_target_coverage_metrics(input_bam, output_format)

@follows(index)
@transform(merge_lanes, suffix('.bam'), '.gene_coverage.sample_summary', r'\1.gene_coverage')
def qc_raw_bam_gene_coverage_metrics(input_bam, output, output_format):
    bam_gene_coverage_metrics(input_bam, output_format)

@follows(index, mkdir(os.path.join(runs_scratch_dir,'qc')), mkdir(os.path.join(runs_scratch_dir,'qc','qualimap')))
@transform(merge_lanes, formatter(".*/(?P<SAMPLE_ID>[^/]+).bam"), '{subpath[0][1]}/qc/qualimap/{SAMPLE_ID[0]}')
def qc_raw_bam_qualimap_report(input_bam, output_dir):
    qualimap_bam(input_bam, output_dir)

@follows(qc_raw_bam_target_coverage_metrics, qc_raw_bam_qualimap_report)
def raw_bam_qc():
    """ Aggregates raw bam quality control steps """
    pass





#
#
# Prepare the bam files for variant calling
#


@follows(index)
@transform(merge_lanes, suffix(".bam"), '.dedup.bam')
def remove_dups(bam, output):
    """Use Picard to mark duplicates"""
    args = "MarkDuplicates \
            TMP_DIR={tmp} \
            REMOVE_DUPLICATES=true \
            INPUT={bam} \
            OUTPUT={out} \
            METRICS_FILE={bam}.dup_metrics \
            VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR CREATE_INDEX=true \
            ".format(tmp=tmp_dir, 
                     bam=bam, 
                     out=output)
    run_cmd(picard, args, interpreter_args="-Xmx4g", mem_per_cpu=4096)



#
# Call variants
#


def call_variants_freebayes(bams_list, vcf, ref_genome, targets, bam_list_filename='/tmp/bam_list'):
    
    threads = 1
    mem = 4096
    
    with open(bam_list_filename,'w') as f:
        for bam in bams_list:
            f.write(bam + '\n')
    
    args = args = " -f {ref} -v {vcf} -L {bam_list} -t {target} --report-monomorphic \
        ".format(ref=ref_genome, vcf=vcf, bam_list=bam_list_filename, target=targets)
            
    run_cmd(freebayes, args, cpus=threads, mem_per_cpu=int(mem/threads))
    
    os.remove(bam_list_filename)


@merge(map_trimmed_reads, os.path.join(runs_scratch_dir, "multisample.vcf"))
def jointcall_variants(bams, vcf):
    """ Call variants using freebayes on trimmed (not merged) reads """
    call_variants_freebayes(bams, vcf, reference, capture_plus)




def split_snp_parameters():
    multisample_vcf = os.path.join(runs_scratch_dir, run_id+'.multisample.vcf')
    for s_id in get_sample_ids():
        yield [multisample_vcf, os.path.join(runs_scratch_dir, s_id, s_id + '.vcf'), s_id]


@follows(jointcall_variants)
@files(split_snp_parameters)
def split_snps(vcf, output, sample):
    """ Split variants by sample, and use sample-specific statistics to filter: AD and DP"""
    AD_threshold=5
    DP_threshold=8
    # what if there is multiple alt alleles??? is the first one most covered      
    args = "-T SelectVariants \
            -R {ref} \
            --variant {vcf} \
            -sn {sample} \
            -select 'vc.getGenotype(\\\"{sample}\\\").getAD().1 >= {ad_thr} && vc.getGenotype(\\\"{sample}\\\").getDP() >= {dp_thr}' \
            -o {out} \
            ".format(ref=reference,
                     vcf=vcf, 
                     sample=sample, 
                     out=output, 
                     ad_thr=AD_threshold, 
                     dp_thr=DP_threshold)
            
    run_cmd(gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx2g" % tmp_dir, mem_per_cpu=2048)


@follows(mkdir(os.path.join(runs_scratch_dir,'qc')))
@merge(split_snps, os.path.join(runs_scratch_dir,'qc','variant_qc'))
def variants_qc(vcf, output):
    """ Generate variant QC table for all samples """    
    args = "-T VariantEval \
            -R {ref} \
            -o {out} \
            -noST -noEV -EV CountVariants \
            ".format(ref=reference,
                     vcf=vcf, 
                     sample=sample, 
                     out=output, 
                     ad_thr=AD_threshold, 
                     dp_thr=DP_threshold)
    
    for vcf in vcfs:
        args += " --eval:{sample} {vcf}" % (os.path.basename(vcf), vcf)
        
    run_cmd(gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx2g" % tmp_dir, mem_per_cpu=2048)
    

def archive_results():
    # if optional results_archive was not provided - do nothing
    if results_archive == None: return
    arch_path = os.path.join(results_archive, run_id)
    if not os.path.exists(arch_path): 
        os.mkdir(arch_path)
        
    run_cmd("cp %s/*/*.bam %s" % (runs_scratch_dir,arch_path), "", run_locally=True)
    run_cmd("cp %s/*/*.bam.gene_coverage* %s" % (runs_scratch_dir,arch_path), "", run_locally=True)
    run_cmd("cp %s/*/*.vcf %s" % (runs_scratch_dir,arch_path), "", run_locally=True)
    run_cmd("cp -r %s/qc %s" % (runs_scratch_dir,arch_path), "", run_locally=True)


def cleanup_files():
    run_cmd("rm -rf {dir}/*/*.recal_data.csv {dir}/*/*.realign* {dir}/*/*.dedup* \
            {dir}/*.multisample.indel.model* {dir}/*.multisample.snp.model* \
            {dir}/*/*.log {dir}/*.multisample.recalibratedSNPS.rawIndels.vcf* \
            {dir}/*.multisample.recalibrated.vcf* \
            ".format(dir=runs_scratch_dir), "", run_locally=True)


@posttask(archive_results, cleanup_files)
@follows(raw_bam_qc, variants_qc)
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
                            logger          = stderr_logger,
                            verbose         = options.verbose,
                            gnu_make_maximal_rebuild_mode = options.rebuild_mode,
                            checksum_level  = 0)
    
        
    drmaa_session.exit()
    
