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

#   Configuration / options

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


if __name__ == '__main__':
    
    from ruffus.proxy_logger import *
    from ruffus import *
    import drmaa
    
    
    import pipeline.config
   
    # parse command-line args, and check for mandatory options
    options, helpstr = pipeline.config.parse_cli_args()
    pipeline.config.check_mandatory_options(options, 
                                        mandatory_options = ['run_folder'], 
                                        help_text=helpstr)
    
    # configure logging
    logger = pipeline.config.setup_logging(options.log_file, options.verbose)
    # allow logging across Ruffus pipeline
    def get_logger(logger_name, args):
        return logger

    (logger_proxy, logging_mutex) = \
        make_shared_logger_and_proxy(get_logger, logger.name, {})


    # Get pipeline settings from a config file  
    cfg = pipeline.config.PipelineConfig() 
    cfg.set_logger(logger)
    cfg.set_runfolder(options.run_folder)
    cfg.set_num_jobs(options.jobs)
    cfg.load_settings_from_file(options.pipeline_settings if options.pipeline_settings != None 
                                                        else os.path.join(options.run_folder, 'settings.cfg') )

    # init drmaa
    cfg.drmaa_session = drmaa.Session()
    cfg.drmaa_session.initialize()




#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#  Common functions 

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888


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
def run_cmd(cfg, cmd, args, interpreter_args=None, run_locally=True,
            cpus=1, mem_per_cpu=1024, walltime='24:00:00', 
            retain_job_scripts = True, job_script_dir = None):
    
    if job_script_dir is None:
        job_script_dir = os.path.join(cfg.runs_scratch_dir, "drmaa")

    full_cmd = "nice "+cmd.format(args=args, 
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
                                 logger=cfg.logger, working_directory=os.getcwd(),
                                 drmaa_session = cfg.drmaa_session)
    except error_drmaa_job as err:
        raise Exception("\n".join(map(str, ["Failed to run:", cmd, err, stdout, stderr])))



""" 
Only default job scheduling params of run_command available when executing via SLURM.
"""
def run_piped_command(cfg, *args):
    run_locally=True
    retain_job_scripts = True
    job_script_dir = os.path.join(cfg.runs_scratch_dir, "drmaa")	
    cpus=1
    mem_per_cpu=1024
    walltime="24:00:00"
 
    stdout, stderr = "", ""
    job_options = "--ntasks=1 \
                   --cpus-per-task={cpus} \
                   --mem-per-cpu={mem} \
                   --time={time} \
                  ".format(cpus=cpus, mem=int(1.2*mem_per_cpu), time=walltime)
	
    full_cmd = "nice " + expand_piped_command(*args)
    print full_cmd	
    try:
        stdout, stderr = run_job(full_cmd.strip(), 
                                 job_other_options=job_options,
                                 run_locally = run_locally, 
                                 retain_job_scripts = retain_job_scripts, job_script_directory = job_script_dir,
                                 logger=cfg.logger, working_directory=os.getcwd(),
                                 drmaa_session = cfg.drmaa_session)
    except error_drmaa_job as err:
        raise Exception("\n".join(map(str, ["Failed to run:", full_cmd, err, stdout, stderr])))
	
def expand_piped_command(cmd, cmd_args, interpreter_args=None, *args):
	expanded_cmd = cmd.format(args=cmd_args, interpreter_args = interpreter_args if interpreter_args!=None else "")
	expanded_cmd += (" | "+expand_piped_command(*args)) if len(args) > 0 else ""
	return expanded_cmd


def produce_fastqc_report(cfg, fastq_file, output_dir=None):
    args = fastq_file
    args += (' -o '+output_dir) if output_dir != None else ''
    run_cmd(cfg, cfg.fastqc, args)



                            
def index_bam(bam):
    """Use samtools to create an index for the bam file"""
    run_cmd(cfg, cfg.samtools, "index %s" % bam)
                          
def bam_quality_score_distribution(bam,qs,pdf):
    """Calculates quality score distribution histograms"""
    run_cmd(cfg, cfg.picard, "QualityScoreDistribution \
                    CHART_OUTPUT={chart} \
                    OUTPUT={out} \
                    INPUT={bam} \
                    VALIDATION_STRINGENCY=SILENT \
                    ".format(chart=pdf, out=qs, bam=bam),
            interpreter_args="")


def get_sample_ids():
    """ Provides meaningful result only after HaplotypeCaller step"""
    files = glob.glob(os.path.join(cfg.runs_scratch_dir,'*','*.gvcf'))
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

@follows(mkdir(cfg.runs_scratch_dir), mkdir(os.path.join(cfg.runs_scratch_dir,'fastqs')))
@files(cfg.run_folder, os.path.join(cfg.runs_scratch_dir,'fastqs','completed'))
@posttask(touch_file(os.path.join(cfg.runs_scratch_dir,'fastqs','completed')))
def bcl2fastq_conversion(run_directory, completed_flag):
    """ Run bcl2fastq conversion and create fastq files in the run directory"""
    out_dir = os.path.join(cfg.runs_scratch_dir,'fastqs')
    interop_dir = os.path.join(out_dir,'InterOp')

    # r, w, d, and p specify numbers of threads to be used for each of the concurrent subtasks of the conversion (see bcl2fastq manual) 
    #args = "-R {indir} -o {outdir} --interop-dir={interopdir} -r1 -w1 -d2 -p4 \
    args = "-R {indir} -o {outdir} --no-lane-splitting \
            ".format(indir=run_directory, outdir=out_dir, interopdir=interop_dir)
    if cfg.run_on_bcl_tile != None:
        args += " --tiles %s" % cfg.run_on_bcl_tile
    run_cmd(cfg, cfg.bcl2fastq, args, cpus=8, mem_per_cpu=2048)
    


@active_if(cfg.fastq_archive != None)
@transform(bcl2fastq_conversion, formatter(".+/(?P<RUN_ID>[^/]+)/fastqs/completed"), str(cfg.fastq_archive)+"/{RUN_ID[0]}")
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
@subdivide(os.path.join(cfg.runs_scratch_dir,'fastqs','*.fastq.gz'),
           formatter('(?P<PATH>.+)/(?P<SAMPLE_ID>[^/]+)_S[1-9]\d?_R[12]_001\.fastq\.gz$'), 
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
#@collate(link_fastqs, regex(r'(.+)/([^/]+)_S[1-9]\d?_R[12]_001\.fastq\.gz$'),  r'\1/\2.fq1.gz')
@collate(link_fastqs, regex(r'(.+)/([^/]+)_S[1-9]\d?_R[12]_001\.fastq\.gz$'),  (r'\1/\2.fq1.gz',r'\1/\2.fq2.gz'))
def trim_reads(inputs, outputs):
    unpaired = [outputs[0].replace('fq1.gz','fq1_unpaired.gz'), outputs[1].replace('fq2.gz','fq2_unpaired.gz')]               
    # logfile = output.replace('fq1.gz','trimmomatic.log')
    # -trimlog {log} \
    # log=logfile
    args = "PE -phred33 -threads 1 \
            {in1} {in2} {out1} {unpaired1} {out2} {unpaired2} \
            ILLUMINACLIP:{adapter}:2:30:10 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36".format(in1=inputs[0], in2=inputs[1],
                              out1=outputs[0], out2=outputs[1],
                              unpaired1=unpaired[0], unpaired2=unpaired[1],
                              adapter=cfg.adapters)
    max_mem = 2048
    run_cmd(cfg, cfg.trimmomatic, args, interpreter_args="-Xmx"+str(max_mem)+"m", cpus=1, mem_per_cpu=max_mem)


###########################
#
# Align/map reads and create BAM files (one per sample)
# 
#########################

def bwa_map_and_sort(output_bam, ref_genome, fq1, fq2=None, read_group=None, threads=1):
	
	bwa_args = "mem -t {threads} {rg} {ref} {fq1} \
	            ".format(threads=threads, 
                        rg="-R '%s'" % read_group if read_group!=None else "", 
                        ref=ref_genome, fq1=fq1)
	if fq2 != None:
		bwa_args += fq2

	samtools_args = "sort -o {out}".format(out=output_bam)

	run_piped_command(cfg, cfg.bwa, bwa_args, None,
	                       cfg.samtools, samtools_args, None)


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
#@transform(trim_reads, formatter("(.+)/(?P<SAMPLE_ID>[^/]+)_(?P<LANE_ID>L\d\d\d).fq[12]\.gz$"), 
#        "{subpath[0][0]}/{SAMPLE_ID[0]}_{LANE_ID[0]}.bam", "{SAMPLE_ID[0]}", "{LANE_ID[0]}")
@jobs_limit(4)
@transform(trim_reads, formatter("(.+)/(?P<SAMPLE_ID>[^/]+).fq[12]\.gz$"), 
        "{subpath[0][0]}/{SAMPLE_ID[0]}.bam", "{SAMPLE_ID[0]}")
def align_reads(fastqs, bam, sample_id):  
    # construct read group information from fastq file name (assuming [SAMPLE_ID]_[LANE_ID].fq[1|2].gz format)
    #sample_lane = os.path.basename(fastqs[0])[0:-len(".fq1.gz")]
    #sample = "_".join(sample_lane.split("_")[0:-1])
    #lane = sample_lane.split("_")[-1]
    sample=sample_id
    read_group = "@RG\\tID:{id}\\tSM:{sm}\\tLB:{lb}\\tPL:{pl}\
                 ".format(id=sample, sm=sample, lb=sample, pl="ILLUMINA")
                 
    #args = "mem -t {threads} -R {rg} {ref} {fq1} {fq2} \
	#    ".format(threads=threads, rg=read_group, ref=reference, 
    #                 fq1=fastqs[0], fq2=fastqs[1])
    #iargs = "samtools view -b -o {bam} -".format(bam=bam)

    #run_cmd(cfg, bwa, args, interpreter_args=iargs, cpus=threads, mem_per_cpu=8192/threads)
    
    bwa_map_and_sort(bam, cfg.reference, fastqs[0], fastqs[1], read_group=read_group, threads=4)


#
# TODO can we merge in the align reads function?
#


#
# BAM filenames are expected to have following format:
#    [SAMPLE_ID]_[LANE_ID].bam
# In this step, all BAM files matching on the SAMPLE_ID will be merged into one BAM file:
#    [SAMPLE_ID].bam
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"
#
"""
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
        
    run_cmd(cfg, cfg.picard, args, interpreter_args="-Xmx8g", cpus=4, mem_per_cpu=2048)
    

def clean_fastqs_and_lane_bams():
    "" Remove the trimmed fastq files, and SAM files. Links to original fastqs are kept ""
    for f in glob.glob(os.path.join(cfg.runs_scratch_dir,'*','*.fq[12]*.gz')):
        os.remove(f)
    for f in glob.glob(os.path.join(cfg.runs_scratch_dir,'*','*_L\d\d\d.bam')):
        os.remove(f)
"""

#@posttask(clean_fastqs_and_lane_bams)
@transform(align_reads, suffix(".bam"), '.bam.bai')
def index(bam, output):
    """Create raw bam index"""
    index_bam(bam)



############################
#
# QC the bam files
#
#####################

#
#TODO: requires update. Use different tools?
#

#@follows(index)
#@transform(merge_lanes, suffix(".bam"), '.quality_score')
#def qc_raw_bam_quality_score_distribution(input_bam, output):
    #"""docstring for metrics1"""
    #bam_quality_score_distribution(input_bam, output, output + '.pdf')

@follows(index)
@transform(align_reads, suffix(".bam"), '.bam.alignment_metrics')
def qc_bam_alignment_metrics(input_bam, output):
    """Collects alignment metrics for a bam file"""
    run_cmd(cfg, cfg.picard, "CollectAlignmentSummaryMetrics \
                    REFERENCE_SEQUENCE={ref} \
                    OUTPUT={out} \
                    INPUT={bam} \
                    VALIDATION_STRINGENCY=SILENT \
                    ".format(ref=cfg.reference, out=metrics, bam=bam), interpreter_args="-Xmx2g")
                    

@follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc')))                    
@merge(qc_bam_alignment_metrics,os.path.join(cfg.runs_scratch_dir, 'qc', 'all_samples.alignment_metrics'))
def qc_aggregate_alignment_metrics(inputs,output):

	import tsv
	import os
	list_of_names = []
	categories = ''
	num = 0

	writer = tsv.TsvWriter(open(output, "w"))

	for name in inputs:
		base = os.path.basename(name)
		list_of_names.append(base[:-len('.bam_alignment_metrics')])
		with open(name, 'r') as f:
			for line in f:
				line = line.strip('\n')
				if categories == '' and line.startswith('CATEGORY'):
					categories = line.split('\t')
                                        categories[0]='SAMPLE'
					writer.list_line(categories)
				
				elif line.startswith('PAIR'):
					splited = line.split('\t')
					splited.remove('PAIR')
					splited.insert(0,list_of_names[num])
					num += 1
					writer.list_line(splited)
	writer.close()
    
    
    
@follows(index)
@transform(align_reads, suffix(".bam"), '.target_coverage.sample_summary', r'\1.target_coverage')
def qc_bam_target_coverage_metrics(input_bam, output, output_format):
    """ Calculate coverage on target """
    run_cmd(cfg, cfg.gatk, "-R {reference} \
                    -T DepthOfCoverage \
                    -o {output} \
                    -I {input} \
                    -L {capture} \
                    -ct 8 -ct 20 -ct 30 \
                    --omitDepthOutputAtEachBase --omitLocusTable \
                    ".format(reference=cfg.reference,
                             output=output_format,
                             input=input_bam,
                             capture=cfg.capture),
            interpreter_args="-Xmx4g")

@follows(index)
#@transform(align_reads, suffix('.bam'), '.gene_coverage.sample_summary', r'\1.gene_coverage')
@merge(align_reads, os.path.join(cfg.runs_scratch_dir, 'qc', 'all_samples.coverage'))
def qc_bam_gene_coverage_metrics(input_bams, output):
    """Calculates and outputs bam coverage statistics """
    bam_list_file = os.path.join(cfg.tmp_dir, 'file_with_bam_lists.list')
    with open(bam_list_file,'w') as f:
        for bam_path in input_bams:
            f.write(bam_path+'\n')
        
    run_cmd(cfg, cfg.gatk, "-R {reference} \
                    -T DepthOfCoverage \
                    -o {output} \
                    -I {inputs} \
                    -L {capture} \
                    -geneList {genes} \
                    -ct 10 -ct 20 \
                    --omitDepthOutputAtEachBase --omitLocusTable \
                    ".format(reference=cfg.reference,
                             output=output,
                             inputs=bam_list_file,
                             capture=cfg.capture,
                             genes=cfg.gene_coordinates),
            interpreter_args="-Xmx4g")
            
           #os.remove(bam_list_file)

@transform(align_reads, formatter(".*/(?P<SAMPLE_ID>[^/]+).bam"), '{subpath[0][1]}/qc/qualimap/{SAMPLE_ID[0]}')
def qc_bam_qualimap_report(input_bam, output_dir):
    """ Generates Qualimap bam QC report """
    if not os.path.exists(output_dir): os.mkdir(output_dir)
    run_cmd(cfg, cfg.qualimap, "bamqc -bam {bam} \
                        -c -outformat PDF \
                        -gff {target} \
                        -gd HUMAN -os \
                        -outdir {dir} &> {dir}/qualimap.err \
                        ".format(bam=input_bam,
                                target=cfg.capture_qualimap,
                                dir=output_dir),
            interpreter_args="")


@follows(qc_bam_alignment_metrics, qc_bam_target_coverage_metrics, qc_bam_qualimap_report)
def bam_qc():
    """ Aggregates raw bam quality control steps """
    pass





# ###################
#
# Prepare the bam files for variant calling
#
#################

#
# TODO: can we mark/remove dups on the fly, e.g. with sambamba. Samtools requires fixmate which requires name-sorted order
#

@follows(index)
@transform(align_reads, suffix(".bam"), '.dedup.bam')
def mark_dups(bam, output):
    """Use Picard to mark duplicates"""
    args = "MarkDuplicates \
            TMP_DIR={tmp} \
            REMOVE_DUPLICATES=true \
            INPUT={bam} \
            OUTPUT={out} \
            METRICS_FILE={bam}.dup_metrics \
            VALIDATION_STRINGENCY=LENIENT VERBOSITY=ERROR CREATE_INDEX=true \
            ".format(tmp=cfg.tmp_dir, 
                     bam=bam, 
                     out=output)
    run_cmd(cfg, cfg.picard, args, interpreter_args="-Xmx2g", mem_per_cpu=4096)


#####################
#
# Call variants
#
################3#

def call_variants_freebayes(bams_list, vcf, ref_genome, targets, bam_list_filename='/tmp/bam_list'):
    
    threads = 1
    mem = 4096
    
    with open(bam_list_filename,'w') as f:
        for bam in bams_list:
            f.write(bam + '\n')
    
#    args = args = " -f {ref} -v {vcf} -L {bam_list} -t {target} --report-monomorphic \
#        ".format(ref=ref_genome, vcf=vcf, bam_list=bam_list_filename, target=targets)
 
    args = args = " -f {ref} -v {vcf} -L {bam_list} -t {target} \
                  --min-coverage {min_dp} --min-alternate-count {min_alt_dp} \
                  --min-mapping-quality {mq} --report-monomorphic \
                  ".format(ref=ref_genome, 
                           vcf=vcf, 
                           bam_list=bam_list_filename, 
                           target=cfg.capture_plus,
                           min_dp     = 20,
                           min_alt_dp = 4,
                           mq         = 30)           

    run_cmd(cfg, cfg.freebayes, args, cpus=threads, mem_per_cpu=int(mem/threads))
    
    os.remove(bam_list_filename)


@merge(mark_dups, os.path.join(cfg.runs_scratch_dir, "multisample.fb.vcf"))
def jointcall_variants(bams, vcf):
    """ Call variants using freebayes on trimmed (not merged) reads """
    call_variants_freebayes(bams, vcf, cfg.reference, cfg.capture_plus)



#
# TODO:
# pcr_indel_model based on settings
#
@jobs_limit(16)
@transform(mark_dups, suffix('.bam'), '.gvcf')
def call_haplotypes(bam, output_gvcf):
    """Perform variant calling using GATK HaplotypeCaller"""
    cmd = "-T HaplotypeCaller \
            -R {ref} \
            -I {bam} \
            -o {gvcf} \
            --emitRefConfidence GVCF \
            -L {target} \
            -pcr_indel_model {pcr} \
            ".format(ref=cfg.reference, 
                     bam=bam, 
                     gvcf=output_gvcf, 
                     target=cfg.capture_plus,
                     pcr='NONE')

    run_cmd(cfg, cfg.gatk, args, '-Djava.io.tmpdir=%s -Xmx4g' % cfg.tmp_dir)


@merge(call_haplotypes, os.path.join(cfg.runs_scratch_dir, 'multisample.gatk.vcf'))
def genotype_gvcfs(gvcfs, output):
    """Combine the per-sample GVCF files and genotype""" 
    args = "-T GenotypeGVCFs \
            -R {ref} \
            -o {out} \
            -dcov 3000 \
            -nt {threads}" % (ref=cfg.reference, out=output, threads=cfg.num_jobs)
       
    for gvcf in gvcfs:
        args += " --variant %s" % gvcf
    
    run_cmd(cfg, cfg.gatk, args, '-Xmx16g')





#
# TODO - how to split variants in elegantly?
#
def split_snp_parameters():
    multisample_vcf = os.path.join(cfg.runs_scratch_dir, cfg.run_id+'.multisample.vcf')
    for s_id in get_sample_ids():
        yield [multisample_vcf, os.path.join(cfg.runs_scratch_dir, s_id, s_id + '.vcf'), s_id]


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
            ".format(ref=cfg.reference,
                     vcf=vcf, 
                     sample=sample, 
                     out=output, 
                     ad_thr=AD_threshold, 
                     dp_thr=DP_threshold)
            
    run_cmd(cfg, cfg.gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx2g" % cfg.tmp_dir, mem_per_cpu=2048)


#
# TODO - update variant stats generation
#
@follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc')))
@merge(split_snps, os.path.join(cfg.runs_scratch_dir,'qc','variant_qc'))
def variants_qc(vcf, output):
    """ Generate variant QC table for all samples """    
    args = "-T VariantEval \
            -R {ref} \
            -o {out} \
            -noST -noEV -EV CountVariants \
            ".format(ref=cfg.reference,
                     vcf=vcf, 
                     sample=sample, 
                     out=output, 
                     ad_thr=AD_threshold, 
                     dp_thr=DP_threshold)
    
    for vcf in vcfs:
        args += " --eval:{sample} {vcf}" % (os.path.basename(vcf), vcf)
        
    run_cmd(cfg, cfg.gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx2g" % cfg.tmp_dir, mem_per_cpu=2048)
    

def archive_results():
    # if optional results_archive was not provided - do nothing
    if cfg.results_archive == None: return
    arch_path = os.path.join(cfg.results_archive, cfg.run_id)
    if not os.path.exists(arch_path): 
        os.mkdir(arch_path)
        
    run_cmd(cfg, "cp %s/*/*.bam %s" % (cfg.runs_scratch_dir,arch_path), "", run_locally=True)
    run_cmd(cfg, "cp %s/*/*.bam.gene_coverage* %s" % (cfg.runs_scratch_dir,arch_path), "", run_locally=True)
    run_cmd(cfg, "cp %s/*/*.vcf %s" % (cfg.runs_scratch_dir,arch_path), "", run_locally=True)
    run_cmd(cfg, "cp -r %s/qc %s" % (cfg.runs_scratch_dir,arch_path), "", run_locally=True)


def cleanup_files():
    run_cmd(cfg, "rm -rf {dir}/*/*.recal_data.csv {dir}/*/*.realign* {dir}/*/*.dedup* \
            {dir}/*.multisample.indel.model* {dir}/*.multisample.snp.model* \
            {dir}/*/*.log {dir}/*.multisample.recalibratedSNPS.rawIndels.vcf* \
            {dir}/*.multisample.recalibrated.vcf* \
            ".format(dir=cfg.runs_scratch_dir), "", run_locally=True)


@posttask(archive_results, cleanup_files)
@follows(bam_qc, variants_qc)
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
    
        
    cfg.drmaa_session.exit()
    
