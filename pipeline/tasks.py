
#
#
# Functions implementing recurring pipeline tasks
#

import os
from utils import run_cmd, run_piped_command
from pipeline.config import PipelineConfig

cfg = PipelineConfig.getInstance()

#
# alignment
#
def bwa_map_and_sort(output_bam, ref_genome, fq1, fq2=None, read_group=None, threads=1):
	
	bwa_args = "mem -t {threads} {rg} {ref} {fq1} \
	            ".format(threads=threads, 
                        rg="-R '%s'" % read_group if read_group!=None else "", 
                        ref=ref_genome, fq1=fq1)
	if fq2 != None:
		bwa_args += fq2

	samtools_args = "sort -o {out}".format(out=output_bam)

	run_piped_command(cfg.bwa, bwa_args, None,
	                  cfg.samtools, samtools_args, None)
                           
                           
def speedseq_align(output_prefix, read_group, ref_genome, fq1, fq2=None, threads=8):
    args = "align -t {threads} -o {out} -T {out}.tmp -R '{rg}' {ref} {fq1} {fq2} \
            ".format(out=output_prefix, rg=read_group, 
                     ref=ref_genome, fq1=fq1, 
                     fq2="" if fq2 is None else fq2,
                     threads = threads)
    
    run_cmd(cfg.speedseq, args, None)



def samtools_index(bam):
    """Use samtools to create an index for the bam file"""
    run_cmd(PipelineConfig.getInstance().samtools, "index %s" % bam)



def recalibrateBQ(bam, metrics_table):
    """ GATK BQSR """
    run_cmd(cfg.gatk, "BaseRecalibrator \
                        -I {bam} \
                        -R {ref} \
                        --known-sites {sites} \
                        -O {table} \
                        ".format(bam=bam, 
                                ref=cfg.ref_genome,
                                sites=cfg.dbsnp_vcf,
                                table=metrics_table))
                                


#
# Calling
#

def speedseq_var(output_prefix, ref_genome, bams, include_bed=None, annotate=False, threads=16):
    args = "var -k -o {out} -T {out}.tmp {bed} {annotate} -t {threads} {ref}\
           ".format(out=output_prefix, 
                    annotate="-A " if annotate else "",
                    bed = "" if include_bed is None else "-w "+include_bed,
                    threads = threads,
                    ref = ref_genome)
                    
    for bam in bams:
        args += " "+bam
    
    run_cmd(PipelineConfig.getInstance().speedseq, args, None)


def speedseq_sv(output_prefix, ref_genome, 
                concordant_bams, splitters_bams, discordant_bams, 
                exclude_bed=None, genotype=True, annotate=False, 
                threads=16):
    
    bams = ','.join(concordant_bams)
    splitters = ','.join(splitters_bams)
    discordants = ','.join(discordant_bams)
    
    args = "sv -o {out} -R {ref} {exclude} {genotype} {annotate} \
               -B {bams} -S {splitters} -D {discords} \
               -t {threads}  \
           ".format(out=output_prefix, ref = ref_genome,
                    exclude = "" if exclude_bed is None else "-x "+exclude_bed,
                    genotype="-g " if genotype else "",
                    annotate="-A " if annotate else "",
                    bams = bams, splitters = splitters, discords = discordants,
                    threads = threads)
                        
    run_cmd(PipelineConfig.getInstance().speedseq, args, None)




def cnvnator_sv(bam, calls, genome_dir, bin_size=100, 
                chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 \
                             chr10 chr11 chr12 chr13 chr14 chr15 chr16 \
                            chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"):
    
    ''' executes 5-step cnvnator calling process. Expects that thisroot.sh is sourced in advance '''
    
    #source /tools/root_v6.16.00/bin/thisroot.sh

    rootfile = bam+'.cnvnator.root'
    errfile  = rootfile+'.err'
    bs = str(bin_size)
    
    args1="-root {rf} -tree {bam} -unique -lite -chrom {chrom} \
          >> {err} 2>&1".format(rf=rootfile, bam=bam, 
                                 chrom=chromosomes, err=errfile)
                                      
    args2="-root {rf} -his {binsize} -d {genome_dir} >> {err} 2>&1\
          ".format(rf=rootfile, binsize=bs, 
                   genome_dir=genome_dir, err=errfile)

    args3="-root {rf} -stat {binsize} >> {err} 2>&1\
          ".format(rf=rootfile, binsize=bs, err=errfile)
          
    args4="-root {rf} -partition {binsize} >> {err} 2>&1\
          ".format(rf=rootfile, binsize=bs, err=errfile)

    args5="-root {rf} -call {binsize} > {calls} 2>>{err}\
          ".format(rf=rootfile, binsize=bs, calls=calls, err=errfile)

    args = [ args1, args2, args3, args4, args5 ]
    for arg in args:
        run_cmd(PipelineConfig.getInstance().cnvnator, arg, None)

    # cleanup
    os.remove(rootfile)


def cnvnator_calls2bed(calls, bed, bam=None):
    ''' reformats cnvnator calls file to a BED including a column with mean MQ '''

    if bam: # calculating mean MQ for regions
        mqs_file = calls+".mqs"
        if os.path.exists(mqs_file): os.remove(mqs_file)
        with open(calls, 'r') as calls_file:
            for l in calls_file:
                reg = l.strip().split('\t')[1]
                run_piped_command(PipelineConfig.getInstance().samtools, "view {bam} {reg}".format(bam=bam, reg=reg), None,
                                  "awk {args}", "\'{sum+=$5}END{if (NR==0) {print \"NA\"} else {print sum/NR}}\' >> %s" % mqs_file)
    
    with open(calls, 'r') as calls_file, \
         open(bed, 'w') as bed:
         
        mqs = open(calls+'.mqs', 'r') if bam else None
        for l in calls_file:
             ls = l.strip().split('\t')
             
             c,s,e = ls[1].replace(':','\t').replace('-','\t').split('\t')
             cnvtype = ls[0].replace("deletion","DEL").replace("duplication","DUP")
             
             bed.write('\t'.join([c, s, e, ls[2], ls[3], cnvtype] + \
                                 ls[4:10] + \
                                 [mqs.readline() if mqs else "\n"]))
        
        mqs.close()
            




#
# QC
#

def produce_fastqc_report(fastq_file, output_dir=None):
    args = fastq_file
    args += (' -o '+output_dir) if output_dir != None else ''
    run_cmd(cfg.fastqc, args)

                                                      
def bam_quality_score_distribution(bam,qs,pdf):
    """Calculates quality score distribution histograms"""
    run_cmd(cfg.picard, "QualityScoreDistribution \
                    CHART_OUTPUT={chart} \
                    OUTPUT={out} \
                    INPUT={bam} \
                    VALIDATION_STRINGENCY=SILENT \
                    ".format(chart=pdf, out=qs, bam=bam),
            interpreter_args="")
