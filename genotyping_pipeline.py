#!/usr/bin/env python
"""

    genotyping_pipeline.py
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
    import drmaa
    
    
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
    cfg = pipeline.config.PipelineConfig() 
    cfg.set_logger(logger)
    cfg.set_num_jobs(options.jobs)
    cfg.load_settings_from_file(options.pipeline_settings if options.pipeline_settings != None 
                                                        else "settings.cfg")

    # init drmaa
    cfg.drmaa_session = drmaa.Session()
    cfg.drmaa_session.initialize()




#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888

#  Common functions 

#88888888888888888888888888888888888888888888888888888888888888888888888888888888888888888




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


from pipeline.utils import run_cmd, run_piped_command
from pipeline.tasks import produce_fastqc_report, samtools_index, bwa_map_and_sort, bam_quality_score_distribution
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
           '{subpath[0][1]}/{SAMPLE_ID[0]}/{basename[0]}{ext[0]}')
def link_fastqs(fastq_in, fastq_out):
    """Make working directory for every sample and link fastq files in"""
    if not os.path.exists(os.path.dirname(fastq_out)):
        os.mkdir(os.path.dirname(fastq_out))
    if not os.path.exists(fastq_out):
        os.symlink(fastq_in, fastq_out) 

    
   
#
# Input FASTQ filenames are expected to have following format:
#    [SAMPLE_ID]_[R1|R2].fastq.gz
# In this step, the two FASTQ files matching on the [SAMPLE_ID] will be trimmed together (R1 and R2). 
# The output will be written to two FASTQ files
#    [SAMPLE_ID].fq1.gz
#    [SAMPLE_ID].fq2.gz
# SAMPLE_ID can contain all signs except path delimiter, i.e. "\"
#
@collate(link_fastqs, regex(r'(.+)/([^/]+)_R[12].fastq\.gz$'),  (r'\1/\2.fq1.gz',r'\1/\2.fq2.gz'))
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





##################################3
#
#   QC reads
#
############################


@follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc')), mkdir(os.path.join(cfg.runs_scratch_dir,'qc','read_qc')))
@transform(link_fastqs, formatter('.+/(?P<SAMPLE_ID>[^/]+)\.fastq\.gz$'), 
           os.path.join(cfg.runs_scratch_dir,'qc','read_qc/')+'{SAMPLE_ID[0]}_fastqc.html')
def qc_raw_fastq(input_fastq, report):
    """ Generate FastQC report for raw FASTQs """
    produce_fastqc_report(cfg, input_fastq, os.path.dirname(report))

 


@follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc')), mkdir(os.path.join(cfg.runs_scratch_dir,'qc','read_qc')))
@transform(trim_reads, formatter('.+/(?P<SAMPLE_ID>[^/]+)\.fq1\.gz$', '.+/(?P<SAMPLE_ID>[^/]+)\.fq2\.gz$'), 
	  [os.path.join(cfg.runs_scratch_dir,'qc','read_qc')+'/{SAMPLE_ID[0]}.fq1_fastqc.html',
           os.path.join(cfg.runs_scratch_dir,'qc','read_qc')+'/{SAMPLE_ID[1]}.fq2_fastqc.html'])
def qc_trimmed_fastq(input_fastqs, reports):
    """ Generate FastQC report for trimmed FASTQs """
    produce_fastqc_report(cfg, input_fastqs[0], os.path.dirname(reports[0]))
    produce_fastqc_report(cfg, input_fastqs[1], os.path.dirname(reports[1]))





###########################
#
# Align/map reads and create BAM files (one per sample)
# 
#########################


#
# FASTQ filenames are expected to have following format:
#    [SAMPLE_ID].fq[1|2].gz
# In this step, the fq1 file coming from trim_reads is matched with the fq2 file and mapped together. 
# The output will be written to BAM file:
#    [SAMPLE_ID].bam
#
@jobs_limit(4)
@transform(trim_reads, formatter("(.+)/(?P<SAMPLE_ID>[^/]+).fq[12]\.gz$"), 
        "{subpath[0][0]}/{SAMPLE_ID[0]}.bam", "{SAMPLE_ID[0]}")
def align_trimmed_reads(fastqs, bam, sample_id):  
    sample=sample_id
    read_group = "@RG\\tID:{id}\\tSM:{sm}\\tLB:{lb}\\tPL:{pl}\
                 ".format(id=sample, sm=sample, lb=sample, pl="ILLUMINA")
    bwa_map_and_sort(bam, cfg.reference, fastqs[0], fastqs[1], read_group=read_group, threads=4)


@transform(align_reads, suffix(".bam"), '.bam.bai')
def index(bam, output):
    """Create raw bam index"""
    samtools_index(bam)



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
def qc_bam_gene_coverage_metrics_multisample(input_bams, output):
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
    
    
           
@transform(align_reads, suffix('.bam'), '.gene_coverage.sample_gene_summary', r'\1.gene_coverage')
def qc_bam_gene_coverage_metrics_singlesample(input_bam, _, output_prefix):
    """Calculates and outputs bam coverage statistics """      
    run_cmd(cfg, cfg.gatk, "-R {reference} \
                    -T DepthOfCoverage \
                    -o {output} \
                    -I {bam} \
                    -L {capture} \
                    -geneList {genes} \
                    -ct 10 -ct 20 \
                    --omitDepthOutputAtEachBase --omitLocusTable \
                    ".format(reference=cfg.reference,
                             output=output_prefix,
                             bam=input_bam,
                             capture=cfg.capture,
                             genes=cfg.gene_coordinates),
            interpreter_args="-Xmx4g")



@merge(qc_bam_gene_coverage_metrics_singlesample, 
       os.path.join(cfg.runs_scratch_dir, 'qc', 'all_samples.gene_coverage_metrics'))
def qc_bam_aggregate_gene_coverage_metrics(gene_cov_files, output):
    """
    Join several gene_coverage tables, side-to-side.
    Basic assumption is that in all input tables genes are listed in the same order, and that files have equal number of lines/records
    """

    import shutil

    tmp_file_name = os.path.join(cfg.tmp_dir, 'all_samples.gene_coverage_metrics'+'.tmp.'+str(time.time()))

    # drop unecessary columns from first table
    with open(gene_cov_files[0], 'r') as first, \
         open(output, 'w') as merged_gene_cov:
        for l in first.xreadlines():
            l_split = l.split('\t')
            t = [l_split[e] for e in [0,4,8,9]]
            merged_gene_cov.write('\t'.join(t)+'\n')

    # glue rest of files to the right
    for i in range(1, len(gene_cov_files)):
        with open(gene_cov_files[i], 'r') as gene_cov_file, \
             open(output, 'r') as merged_gene_cov, \
             open(tmp_file_name, 'w') as tmp_file:
            
            for line in merged_gene_cov.xreadlines():
                
                new_cols_line_split = gene_cov_file.readline().split('\t')
                new_cols = [new_cols_line_split[e] for e in [4,8,9]]
                
                line = line.strip()
                line += '\t' + '\t'.join(new_cols) + '\n'
                
                tmp_file.write(line)
        
        # replace output with tmp
        shutil.move(tmp_file_name, output)


    


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


@follows(qc_aggregate_alignment_metrics, 
         qc_bam_gene_coverage_metrics_multisample, 
         qc_bam_qualimap_report)
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
#  Lookup the GATK pipeline
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
@transform(mark_dups, suffix('.bam'), '.g.vcf')
def call_haplotypes(bam, output_gvcf):
    """Perform variant calling using GATK HaplotypeCaller"""
    args = "-T HaplotypeCaller \
            -R {ref} \
            -I {bam} \
            -o {gvcf} \
            --emitRefConfidence GVCF \
            -L {target} \
            --pcr_indel_model {pcr} \
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
            -nt {threads} \
            ".format(ref=cfg.reference, out=output, threads=cfg.num_jobs)
       
    for gvcf in gvcfs:
        args += " --variant %s" % gvcf
    
    run_cmd(cfg, cfg.gatk, args, '-Xmx16g')





#
# TODO - how to split variants in elegantly?
#
def split_snp_parameters():
    multisample_vcf = os.path.join(cfg.runs_scratch_dir, "multisample.vcf")
    for s_id in get_sample_ids():
        yield [multisample_vcf, os.path.join(cfg.runs_scratch_dir, s_id, s_id + '.vcf'), s_id]


#@follows(jointcall_variants)
@follows(genotype_gvcfs)
@transform(mark_dups, suffix(".bam"), ".vcf", 
           os.path.join(cfg.runs_scratch_dir, 'multisample.gatk.vcf'))
def split_snps(bam, vcf, multisample_vcf):
    
    
    sample_id = os.path.basename(bam)[:-len('.dedup.bam')]
    
    """ Split variants by sample, and use sample-specific statistics to filter: AD and DP"""
    AD_threshold=5
    DP_threshold=8
    # what if there is multiple alt alleles??? is the first one most covered      
    args = "-T SelectVariants \
            -R {ref} \
            --variant {vcf} \
            -sn {sample} \
            -o {out} \
            ".format(ref=cfg.reference,
                     vcf=multisample_vcf, 
                     sample=sample_id, 
                     out=vcf)
#                     ad_thr=AD_threshold, 
#                     dp_thr=DP_threshold)
# -select 'vc.getGenotype(\\\"{sample}\\\").getAD().1 >= {ad_thr} && vc.getGenotype(\\\"{sample}\\\").getDP() >= {dp_thr}' \            

    run_cmd(cfg, cfg.gatk, args, interpreter_args="-Djava.io.tmpdir=%s -Xmx2g" % cfg.tmp_dir, mem_per_cpu=2048)



##################
#
#  QC variants
#
#############


@follows(mkdir(os.path.join(cfg.runs_scratch_dir,'qc')))
@transform(genotype_gvcfs, formatter(), os.path.join('{subpath[0][0]}','qc','multisample.variant_stats'))
def qc_multisample_vcf(vcf, output):
    """ Generate variant QC table for all samples """    
    args = "stats -F {ref} -s - {vcf} > {out}\
            ".format(ref=cfg.reference, 
                     vcf=vcf,
                     out=output)
    
    run_cmd(cfg, cfg.bcftools, args)
    



#####################
#
#   Annotate variants
#
##################

from pipeline import annotation

@follows(mkdir(os.path.join(cfg.runs_scratch_dir, 'annotation')))
@transform(split_snps, formatter('.*/(?P<SAMPLE_ID>[^/]+).dedup.vcf'), 
           '{subpath[0][1]}/annotation/{SAMPLE_ID[0]}.avinput')
def prepare_annovar_inputs(vcf, output):
    """ convert to annovar format """
    args = "{vcf} -format vcf4 -withzyg -includeinfo -outfile {out} \
        ".format(vcf=vcf, out=output)
    run_cmd(cfg, cfg.convert_to_annovar, args)
    
   
@transform(prepare_annovar_inputs, suffix('.avinput'), 
                                        ['.avinput.hg19_EUR.sites.2015_08_filtered', 
                                         '.avinput.hg19_EUR.sites.2015_08_dropped'])
def filter_common_1000genomes(avinput, outputs):
    """ filter common 1000 genomes variants """
    args = "-build hg19 -filter -dbtype {eur1kg} -maf {maf} \
            -outfile {output_prefix} {input_file} {annodb} \
            ".format(eur1kg=cfg.annovar_1000genomes_eur,
                    maf=cfg.annovar_1000genomes_eur_maf_cutoff, 
                    output_prefix=avinput,
                    input_file=avinput, 
                    annodb=cfg.annovar_human_db)
    
    run_cmd(cfg, cfg.annovar_annotate, args)


#
#
# Currently not used / tested
#
@transform(filter_common_1000genomes, suffix('.hg19_EUR.sites.2014_10_filtered'),
                                        ['.hg19_EUR.sites.2015_08_filtered.common_inhouse_filtered',
                                         '.hg19_EUR.sites.2015_08_filtered.common_inhouse_dropped'])
def filter_common_inhouse(inputs, outputs):
    """ filter variants found in the inhouse databases """
    filtered = inputs[0]              # use only the filtered input file, leave dropped
    annotation.filter_common_inhouse(cfg, filtered, outputs)



#
#
# TODO: Refactor name of the get_stats_on_prefiltered_variants to reflect what it does. 
# TODO: Check if table annotation can work straight on result of filter_common_inhouse/filter_common_1000genomes
#

#@transform(filter_common_inhouse, suffix('.common_inhouse_filtered'), 
#                                           ['.common_inhouse_filtered.variant_function',
#                                         '.common_inhouse_filtered.exonic_variant_function',
#                                         '.common_inhouse_filtered.variant_function.stats',
#                                         '.common_inhouse_filtered.exonic_variant_function.stats'])
@transform(filter_common_1000genomes, suffix('hg19_EUR.sites.2015_08_filtered'), 
                                           ['hg19_EUR.sites.2015_08_filtered.variant_function',
                                         'hg19_EUR.sites.2015_08_filtered.exonic_variant_function',
                                         'hg19_EUR.sites.2015_08_filtered.variant_function.stats',
                                         'hg19_EUR.sites.2015_08_filtered.exonic_variant_function.stats'])
def annotate_function_of_rare_variants(inputs, outputs):
    """ annotate functional change in rare variants """
    filtered = inputs[0]              # use only the filtered input file, leave dropped
    annotation.get_stats_on_prefiltered_variants(cfg, input_file=filtered, outputs=outputs[2:4], cleanup=False)


# with common_inhouse:      formatter(".*/(?P<SAMPLE_ID>[^/]+).avinput.hg19_EUR.sites.2015_08_filtered.common_inhouse_filtered.variant_function", None, None, None),
@transform(annotate_function_of_rare_variants, 
           formatter(".*/(?P<SAMPLE_ID>[^/]+).avinput.hg19_EUR.sites.2015_08_filtered.variant_function", None, None, None),
           ['{subpath[0][1]}/{SAMPLE_ID[0]}/{SAMPLE_ID[0]}.rare_coding_and_splicing.avinput', 
            '{subpath[0][1]}/{SAMPLE_ID[0]}/{SAMPLE_ID[0]}.rare_coding_and_splicing.multianno.csv'])
def produce_variant_annotation_table(inputs, outputs):
    """ produce a table of various annotations per variant """
    annotation.produce_variant_annotation_table(cfg, inputs[0], inputs[1], outputs[0], outputs[1])

#'{path[1]}/{basename[1]}.with_omim.csv'
@transform(produce_variant_annotation_table, suffix('.csv'), '.with_omim.tsv')
def include_omim_phenotype_annotation(inputs, output_table):
    """ add OMIM annotation to the table """
    annotation.include_omim_phenotype_annotation(cfg, inputs[1], output_table)
 

@transform(include_omim_phenotype_annotation, suffix('.tsv'), '.recessive.tsv')
def extract_recessive_disorder_candidates(input_file, output_file):
    """ extract a part of the annotated table that contains candidates for recessive disorders """ 
    annotation.extract_recessive_disorder_candidates(input_file, output_file)

##########################
#
# Annotated variants' QC 
#
######################

@merge(prepare_annovar_inputs, 
       [os.path.join(cfg.runs_scratch_dir, 'qc','multisample.hetz_per_chr'),
        os.path.join(cfg.runs_scratch_dir, 'qc','multisample.homz_per_chr')])
def count_hetz_and_homz_per_chr(inputs, tables):
    """ produce a table of sample vs chromosome counts of hetero- and homozygotic variants """
    annotation.count_hetz_and_homz_per_chr(inputs, tables)


@transform(prepare_annovar_inputs, suffix('.avinput'), ['.avinput.variant_function.stats','.avinput.exonic_variant_function.stats'])
def get_stats_on_raw_variants(input_file, outputs):
    """ annotate functional change in raw variants, get stats, and remove annotated files """
    annotation.get_stats_on_prefiltered_variants(cfg, input_file, outputs, cleanup=True)


@transform(filter_common_1000genomes, suffix('.hg19_EUR.sites.2015_08_filtered'), 
                                        ['.hg19_EUR.sites.2015_08_filtered.variant_function.stats',
                                         '.hg19_EUR.sites.2015_08_filtered.exonic_variant_function.stats'])
def get_stats_on_1kg_filtered_variants(inputs, outputs):
    """ annotate functional change in 1kg filtered variants, get stats, and remove annotated files """
    annotation.get_stats_on_prefiltered_variants(cfg, inputs[0], outputs, cleanup=True)

# equivalent of annotate_rare_variants (input and output files are the same)
@transform(filter_common_inhouse, suffix('.common_inhouse_filtered'), 
                                        ['.common_inhouse_filtered.variant_function.stats',
                                         '.common_inhouse_filtered.exonic_variant_function.stats'])
def get_stats_on_inhouse_filtered_variants(inputs, outputs):
    """ annotate functional change in inhouse-exomes filtered variants, get stats, and remove annotated files """
    annotation.get_stats_on_prefiltered_variants(cfg, inputs[0], outputs, cleanup=True)


#@merge([get_stats_on_raw_variants,
#        get_stats_on_1kg_filtered_variants,
#        get_stats_on_inhouse_filtered_variants],
#        'all_samples_exonic_variant_stats.tsv')
@merge([get_stats_on_raw_variants, 
        get_stats_on_1kg_filtered_variants, 
        get_stats_on_1kg_filtered_variants], 
        os.path.join(cfg.runs_scratch_dir, 'qc','multisample.annotated_variant_stats'))
def produce_variant_stats_table(infiles, table_file):
    """ produce a table of per-sample counts of different type of exonic variants """
    annotation.produce_variant_stats_table(infiles, table_file)


####################
#
#   Cleanup intermediate files
#
##################


def cleanup_files():
    run_cmd(cfg, "rm -rf {dir}/*/*.recal_data.csv {dir}/*/*.realign* {dir}/*/*.dedup* \
            {dir}/*.multisample.indel.model* {dir}/*.multisample.snp.model* \
            {dir}/*/*.log {dir}/*.multisample.recalibratedSNPS.rawIndels.vcf* \
            {dir}/*.multisample.recalibrated.vcf* \
            ".format(dir=cfg.runs_scratch_dir), "", run_locally=True)


@posttask(archive_results, cleanup_files)
@follows(bam_qc, qc_multisample_vcf)
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
    
