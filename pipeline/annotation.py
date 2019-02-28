import os
from utils import get_cfg, run_cmd

cfg = config.PipelineConfig.getInstance()

def filter_common_inhouse(input_variants, outputs):

    inhouse_dbs = cfg.annovar_inhouse_dbs.split(';')

    # take only the filtered file, leave dropped. Copy it to temp file
    filtered = input_variants+'.tmp'     
    shutil.copyfile(input_variants, filtered)
    
    for inhouse_db in inhouse_dbs:
        if inhouse_db.strip() == '':
            continue
            
        args = "-build hg19 -filter -dbtype generic -genericdbfile {inhouse} \
                -outfile {outfile} {infile} {annodb} \
                ".format(inhouse=inhouse_db, 
                         outfile=filtered[:-len('.tmp')], 
                         infile=filtered, 
                         annodb=annovar_human_db)
                         
        run_cmd(cmd, cmd.annovar_annotate, args)
        
        for output_file in outputs: 
            os.rename(output_file.replace('common_inhouse','hg19_generic'), output_file)

        # copy the filtered output to input, for next iteration
        shutil.copyfile(outputs[0], filtered)
    
    # delete the temp file
    os.remove(filtered)

def get_stats_on_prefiltered_variants(input_file, outputs, cleanup=True):    

    args = "-buildver hg19 -outfile {outfile_prefix} {input_file} {annodb}".format(
        tool=cfg.annovar_annotate,
        outfile_prefix=input_file, 
        input_file=input_file, 
        annodb=cfg.annovar_human_db)
        
    run_cmd(cfg.annovar_annotate, args)

    # calculate stats on files created by annovar - output files without ".stats" suffix
    run_cmd("cut -f 1 {f} | sort | uniq -c > {f}.stats".format(f=outputs[0][:-len('.stats')]), "")
    run_cmd("cut -f 2 {f} | sort | uniq -c > {f}.stats".format(f=outputs[1][:-len('.stats')]), "")
    # remove the annovar files
    if cleanup:
        os.remove(outputs[0][:-len('.stats')])
        os.remove(outputs[1][:-len('.stats')])


def produce_variant_annotation_table(variants_file, exonic_variants_file, 
                                    coding_and_splicing_output_file, annotation_table_file):
          
    # create input file for the table_annotation script
    f_out = open(coding_and_splicing_output_file,'w')
    f = open(exonic_variants_file)
    for l in f.xreadlines():
        lsplit=l.split('\t')
        if lsplit[1] != 'synonymous SNV':
            f_out.write('\t'.join(lsplit[3:]))
    f.close()
    f = open(variants_file)
    for l in f.xreadlines():
        lsplit=l.split('\t')
        if 'splicing' in (lsplit[0]).split(';'):
            f_out.write('\t'.join(lsplit[2:]))
    f.close()
    f_out.close()
    
    # annotate all variants selected above
    
    args = "--protocol refGene,1000g2015aug_eur,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_afr,gnomad_genome,avsnp150,clinvar_20180603,ljb26_all \
            --operation g,f,f,f,f,f,f,f,f,f --arg \'--splicing 4\',,,,,,,,, --nastring NA --build hg19 -csvout --otherinfo --remove \
            --outfile {output_prefix} {avinput} {db}".format(
            output_prefix=coding_and_splicing_output_file, 
            avinput=coding_and_splicing_output_file, 
            db=cfg.annovar_human_db)
    
    run_cmd(cfg.table_annovar, args)
    
    # rename annovar output table
    os.rename(coding_and_splicing_output_file+".hg19_multianno.csv", annotation_table_file)


#
# include omim phenotypes

def _get_omim_gene_phenotype_map(omim_file):
    
    print("RUN get_omim_gene_phenotype_map")
    
    gene_col=6
    pht_col=12
    map_pht={}
    f = open(omim_file)
    for l in f.xreadlines():
        lsplit=l.split('|')
        
        # ignore lines with no phenotype
        if lsplit[pht_col-1].strip()=='':
            continue

        genes, phenotype = lsplit[gene_col-1], lsplit[pht_col-1]
        for gene in genes.split(','):
            gene = gene.strip()
            if gene == '': continue
            try:
                map_pht[gene] = map_pht[gene]+'|'+phenotype.strip()
            except KeyError:
                map_pht[gene] = phenotype.strip()
        
    f.close()
    return map_pht


def quote_aware_split(string, delim=',', quote='"'):
    """ Split outside of quotes (i.e. ignore delimiters within quotes."""
    out = []
    s = ''
    open_quote=False
    for c in string:
        if c == quote: 
            open_quote = not open_quote
        if c == delim and not open_quote:
            out += [s]
            s = ''
        else: 
            s += c
    return out + [s]


def parenthesis_aware_split(string, delim=',', open_par='(', close_par=')'):
    """ Split outside of parenthesis (i.e. ignore delimiters within parenthesis."""
    out = []
    s = ''
    open_parenthesis=0
    for c in string:
        if c == open_par: 
            open_parenthesis+=1
        if c == close_par and open_parenthesis > 0:
            open_parenthesis-=1
        if c == delim and open_parenthesis==0:
            out += [s]
            s = ''
        else: 
            s += c
    return out + [s]


def include_omim_phenotype_annotation(input_table, output_table, gene_column=7, omim_column=15, delim=','):
    
    if cfg.omim_gene_phenotype_map == None:
        cfg.omim_gene_phenotype_map = _get_omim_gene_phenotype_map(cfg.omim_gene_phenotype_map_file)

    """ include OMIM phenotype into the annotation table """
    table_in = open(input_table,'r')
    table_out = open(output_table,'w')

    # header
    header_in=quote_aware_split(table_in.readline(), delim)
    if omim_column <= 0:        
        omim_column=len(header_in)+1
    header_out = header_in[:omim_column-1] + ['omim_phenotype'] + header_in[omim_column-1:]
    table_out.write(delim.join(header_out))

    # the rest of the table
    for l in table_in.xreadlines():
        lsplit = quote_aware_split(l,delim)
        gene = lsplit[gene_column-1].strip('"')

        # the gene record can be a list (e.g. overlapping genes), so it needs to be split
        genes = [gene]
        if gene.find(',')>=0:
            genes = parenthesis_aware_split(gene, delim=',')
        genes = [parenthesis_aware_split(gene,delim=';') for gene in genes]
        genes = set([gene for sublist in genes for gene in sublist])  # unlist and get unique gene ids only
        
        for gene in genes:
            # if present, remove suffix in parenthesis
            if gene.find('(') >= 0: 
                gene = gene[:gene.find('(')]
            
            # put the variant record in the map
            try:
                omim_phenotype = cfg.omim_gene_phenotype_map[gene]
            except KeyError:
                omim_phenotype = 'NA'

        table_out.write(delim.join(lsplit[:omim_column-1] + ['"'+omim_phenotype+'"'] + lsplit[omim_column-1:]) )

    table_in.close()
    table_out.close()


def extract_recessive_disorder_candidates(input_file, output, 
                                            gene_column_name='Gene.refGene', 
                                            zygozity_column_name='Otherinfo', 
                                            delim=','):
    """ extract a part of the annotated table that contains candidates for recessive disorders """ 
    table_in = open(input_file,'r')
    
    # get right column indexes
    header = quote_aware_split(table_in.readline().strip(), delim)
    gene_col_index     = header.index(gene_column_name)
    zygozity_col_index = header.index(zygozity_column_name) 
    
    variant_records_per_gene={}
    for l in table_in.xreadlines():
        lsplit = quote_aware_split(l,delim)
        
        gene = lsplit[gene_col_index].strip('"')
        
        # the gene record can be a list (e.g. overlapping genes), so it needs to be split
        genes = [gene]
        if gene.find(',')>=0:
            genes = parenthesis_aware_split(gene, delim=',')
        genes = [parenthesis_aware_split(gene,delim=';') for gene in genes]
        genes = set([gene for sublist in genes for gene in sublist])  # unlist and get unique gene ids only
        
        for gene in genes:
            # if present, remove suffix in parenthesis
            if gene.find('(') >= 0: gene = gene[:gene.find('(')]
            # put the variant record in the map
            try:
                variant_records_per_gene[gene] += [l]
            except KeyError: 
                variant_records_per_gene[gene] = [l]
    
    table_in.close()
    
    # write the table with candidates for recessive inheritance model
    table_out = open(output,'w')
    table_out.write(delim.join(header)+'\n')
    
    # iterate over the genes and select...
    for gene in variant_records_per_gene:
        # ...these with 2 or more variants...
        if len(variant_records_per_gene[gene]) >= 2:
            for l in variant_records_per_gene[gene]: 
                table_out.write(l)
        # or homozygous variants
        else:
            lsplit = quote_aware_split(variant_records_per_gene[gene][0])
            if lsplit[zygozity_col_index].find('"hom\t') >= 0:
                table_out.write(variant_records_per_gene[gene][0])
    
    table_out.close()
    
    

#
# QC on variant level
##

def count_hetz_and_homz_per_chr(infiles, table_files):
    """ produce a table of sample vs chromosome counts of hetero- and homozygotic variants """
    # the variant calls are expected to appear sorted by chromosome name in the following order
    # any other order should trigger an exception
    chromosomes = ['1','2','3','4','5','6','7','8','9',
                   '10','11','12','13','14','15','16','17','18','19',
                   '20','21','22','X','Y','MT',
                   'GL000191.1', 'GL000228.1', 'GL000209.1',  
                   'GL000223.1','GL000222.1', 'GL000194.1']

    hetz = open(table_files[0], 'w')
    homz = open(table_files[1], 'w')
    hetz.write('sample\t'+'\t'.join(chromosomes)+'\n')
    homz.write('sample\t'+'\t'.join(chromosomes)+'\n')
    for fname in infiles:
        sample_id=os.path.basename(fname).split('.')[0]
        hetz.write(sample_id)
        homz.write(sample_id)

        het_cnt=0
        hom_cnt=0
        curr_chr=0
        f = open(fname)
        for l in f.xreadlines():
            lsplit=l.split('\t')
            while lsplit[0] != chromosomes[curr_chr] and curr_chr+1<len(chromosomes):
                hetz.write("\t" + str(het_cnt))
                homz.write("\t" + str(hom_cnt))
                het_cnt=0
                hom_cnt=0
                curr_chr+=1
            # if the call is in a chromosome that is not in the list (or calls could be out of order)
            if lsplit[0] != chromosomes[curr_chr] and curr_chr+1 >= len(chromosomes):
                raise Exception("Unrecognized or out-of-order chromosome: "+lsplit[0])
            
            if lsplit[5] == 'het': het_cnt+=1
            elif lsplit[5] == 'hom': hom_cnt+=1
            else:
                raise Exception("Variant is not a het nor a hom")

        f.close()

        # flush out the counts for the last/remaining chromosomes
        while curr_chr<len(chromosomes):
            hetz.write("\t" + str(het_cnt))
            homz.write("\t" + str(hom_cnt))
            het_cnt=0
            hom_cnt=0             
            curr_chr+=1
        hetz.write('\n')
        homz.write('\n') 

    hetz.close()
    homz.close()
    


def produce_variant_stats_table(infiles, table_file):
    """ produce a table of per-sample counts of different type of exonic variants """

    var_types=['splicing','UTR3','UTR5','intronic','intergenic','exonic']
    
    # split the input files per task
    sample_no = len(infiles)/3
    raw_variant_files = infiles[0:sample_no]
    kg1_filtered_variant_files = infiles[sample_no:sample_no*2]
    inhouse_filtered_variant_files = infiles[sample_no*2:sample_no*3]

    out = open(table_file,'w')
   
    import itertools
    #header = ['sample'] + ['raw_'+t for t in var_types] + ['rare_'+t for t in var_types] + ['raw_synonymous','rare_synonymous']
    filtering_stages = ['raw_','kg1_','inhouse_']
    header = ['sample'] + \
            [f+t for (f,t) in itertools.product(filtering_stages, var_types)] + \
            [s+'synonymous' for s in filtering_stages]
    out.write(('\t'.join(header))+'\n')
    for i in range(0,sample_no):
        out.write(os.path.basename(raw_variant_files[i][0]).split('.')[0])
        for fname in [raw_variant_files[i][0], kg1_filtered_variant_files[i][0], \
                        inhouse_filtered_variant_files[i][0]]: # exonic variant stats of raw variants and rare variants
            counts = dict.fromkeys(var_types,0)
            f=open(fname)        
            for l in f.xreadlines():
                for var_type in var_types:
                    if l.find(' '+var_type)>0:
                        counts[var_type] += int(l.split()[0])
                        break
            out.write('\t'+'\t'.join([str(counts[t]) for t in var_types]))
            f.close()

        for fname in [raw_variant_files[i][1], kg1_filtered_variant_files[i][1], \
                        inhouse_filtered_variant_files[i][1]]: # synonymous variants stats of raw and rare variants
            f=open(fname)
            found=False
            for l in f.xreadlines():
                if l.find(" synonymous")>0:
                    found=True
                    out.write('\t'+l.split()[0])
                    break     
            if not found: out.write('\t0')       
            f.close()
        out.write('\n')
    out.close()

