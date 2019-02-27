
#
#
# Functions implementing recurring pipeline tasks
#

from utils import run_cmd, run_piped_command

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

	run_piped_command(cfg, cfg.bwa, bwa_args, None,
	                       cfg.samtools, samtools_args, None)
                           
                           
def speedseq_align(output_prefix, read_group, ref_genome, fq1, fq2=None, threads=8):
    args = "align -t {threads} -o {out} -R {rg} {ref} {fq1} {fq2} \
            ".format(out=output_prefix, rg=read_group, 
                     ref=ref_genome, fq1=fq1, 
                     fq2="" if fq2 is None else fq2,
                     threads = threads)
        
    run_cmd(cfg, cfg.speedseq, args, None)



def samtools_index(bam):
    """Use samtools to create an index for the bam file"""
    run_cmd(cfg, cfg.samtools, "index %s" % bam)

#
# Calling
#

def speedseq_var(output_prefix, ref_genome, bams, include_bed=None, threads=16):
    args = "var -o {out} {bed} -t {threads} {ref}\
           ".format(out=output_prefix, 
                    bed = "" if include_bed is None else "-w "+include_bed,
                    threads = threads,
                    ref = ref_genome)
                    
    for bam in bams:
        args += " "+bam
    
    run_cmd(cfg, cfg.speedseq, args, None)
    

def speedseq_sv(output_prefix, ref_genome, 
                concordant_bams, splitters_bams, discordant_bams, 
                exclude_bed=None, genotype=True, annotate=True, 
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
                    bams = bams, splitters = splitters, discordants = discordants,
                    threads = threads)
                        
    run_cmd(cfg, cfg.speedseq, args, None)


#
# QC
#

def produce_fastqc_report(cfg, fastq_file, output_dir=None):
    args = fastq_file
    args += (' -o '+output_dir) if output_dir != None else ''
    run_cmd(cfg, cfg.fastqc, args)

                                                      
def bam_quality_score_distribution(bam,qs,pdf):
    """Calculates quality score distribution histograms"""
    run_cmd(cfg, cfg.picard, "QualityScoreDistribution \
                    CHART_OUTPUT={chart} \
                    OUTPUT={out} \
                    INPUT={bam} \
                    VALIDATION_STRINGENCY=SILENT \
                    ".format(chart=pdf, out=qs, bam=bam),
            interpreter_args="")
