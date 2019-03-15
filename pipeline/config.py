#
#
#
# Pipeline configurator
#

import os, json
import glob
import sys
import logging
import logging.handlers
import StringIO


def check_mandatory_options (options, mandatory_options, help_text):
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
                     help_text))
                     
def parse_cli_args():
    
    from optparse import OptionParser

    parser = OptionParser(version="%prog 1.0", usage = "\n\n %prog [--settings pipeline_settings.cfg] [--target_task TASK] [more_options]")
                                
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
                        help="File containing all the settings for the analysis. By default settings.cfg in searched for in the current directory.")                  
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
    
    parser.set_defaults(pipeline_settings=None, 
                        jobs=1, verbose=0, 
                        target_tasks=list(), forced_tasks=list(), 
                        just_print=False, key_legend_in_graph=False,
                        rebuild_mode=True)
    

    # get help string
    f = StringIO.StringIO()
    parser.print_help(f)
    helpstr = f.getvalue()
    
    (options, remaining_args) = parser.parse_args()
 
    return options, helpstr


def _setup_std_logging (logger, log_file, verbose):
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

def setup_logging(log_file, verbosity):

    MESSAGE = 15
    logging.addLevelName(MESSAGE, "MESSAGE")
    
    module_name = "ngspipe"
    logger = logging.getLogger(module_name)
    _setup_std_logging(logger, log_file, verbosity)

    return logger


import ConfigParser

class PipelineConfig:    
    
    instance = None
    
    @staticmethod
    def getInstance():
        if PipelineConfig.instance == None:
            PipelineConfig.instance = PipelineConfig()
        return PipelineConfig.instance


    def __init__(self):
        
        self.instance = None
        
        self.drmaa_session = None
        self.logger = None
        
        # tool paths
        self.speedseq    = None
        self.trimmomatic = None
        self.fastqc      = None
        self.bwa         = None
        self.samtools    = None
        self.picard      = None
        self.gatk        = None
        self.freebayes   = None
        self.bcftools    = None
        self.qualimap    = None

        # annovar settings
        self.convert_to_annovar = None
        self.annovar_annotate   = None
        self.table_annovar      = None
        self.annovar_human_db                   = None
        self.annovar_1000genomes_eur            = None
        self.annovar_1000genomes_eur_maf_cutoff = None
        self.annovar_inhouse_dbs                = None
        self.omim_gene_phenotype_map            = None


        # run paths
        self.reference_root = None
        self.scratch_root = None
        self.input_fastqs = None
        self.runs_scratch_dir = None
        self.tmp_dir      = None

        self.adapters  = None
        self.reference = None
        self.speedseq_include_bed = None
        self.speedseq_lumpy_exclude_bed = None

        # ruffus settings
        self.target_tasks    = []
        self.log_file        = None
        self.num_jobs        = 1
        self.verbosity_level = 0
        self.dry_run         = False
        self.rebuild_mode    = False


    def set_logger(self, logger):
        self.logger = logger
   

    def set_num_jobs(self,jobs):
        if jobs is not None:
            self.num_jobs = jobs


    def set_input_fastqs(self, fastqs):
        self.input_fastqs = fastqs
    
    
    # helper method
    def _get_optional_param(self, cfg, section, param, default_value=None, log_msg=None):
        try:
            return cfg.get(section, param)
        except ConfigParser.NoOptionError:
            if log_msg: self.logger.info(log_msg)
            return default_value

     
    
    def load_settings_from_file(self, cfg_file):
        
        #
        #
        # TODO
        # Missing settings should not cause exceptions
        #
        #
        #
        
        
        
        """ preset settings using config file """

        if not os.path.exists(cfg_file): 
            raise Exception('Provided config file [%s] does not exist or cannot be read.' % cfg_file)

        config = ConfigParser.ConfigParser()
        config.read(cfg_file)
        
        
        self.input_fastqs = glob.glob(config.get('Paths','input-fastqs'))
        
        self.reference_root = config.get('Paths','reference-root')
        
        self.scratch_root = self._get_optional_param(config, 'Paths','scratch-root', \
                            os.getcwd(), 'Scratch-root setting is missing. Using current directory: %s' % os.getcwd())
                  
    
        # workdir
        self.runs_scratch_dir = self.scratch_root
        self.logger.info('Work directory: %s' % self.runs_scratch_dir)
                  
                          
        # optional /tmp dir
        self.tmp_dir = self._get_optional_param(config, 'Paths','tmp-dir', \
                                                '/tmp', 'No tmp-dir provided. /tmp will be used.')
        
        # reference files
        self.reference        = os.path.join(self.reference_root, config.get('Resources', 'reference-genome'))
        self.gene_coordinates = os.path.join(self.reference_root, config.get('Resources', 'gene-coordinates'))       
        self.adapters         = os.path.join(self.reference_root, config.get('Resources', 'adapters-fasta'))
        
        # optional speedseq BED annotations
        beds = [self._get_optional_param(config, 'Resources', val, log_msg='Speedseq\'s BED not set: '+val) \
                 for val in ['speedseq-include-bed', 'speedseq-lumpy-exclude-bed'] ]
        self.speedseq_include_bed, self.speedseq_lumpy_exclude_bed = \
            [ os.path.join(self.reference_root, filename) if filename else None for filename in beds ]
        
        
        # tools
        self.speedseq    = config.get('Tools','speedseq')
        self.trimmomatic = config.get('Tools','trimmomatic') 
        self.bwa         = config.get('Tools','bwa')
        self.samtools    = config.get('Tools','samtools')
        self.picard      = config.get('Tools','picard')
        self.gatk        = config.get('Tools','gatk')
        self.freebayes   = config.get('Tools','freebayes')
        self.bcftools    = config.get('Tools','bcftools')
        self.qualimap    = config.get('Tools','qualimap')
    	self.fastqc	     = config.get('Tools','fastqc')


        # annovar settings
        self.convert_to_annovar                 = os.path.join(config.get('Annovar','annovar_home'), 
                                                               config.get('Annovar','convert_to_annovar'))
        self.annovar_annotate                   = os.path.join(config.get('Annovar','annovar_home'),
                                                               config.get('Annovar','annovar_annotate'))
        self.table_annovar                      = os.path.join(config.get('Annovar','annovar_home'), 
                                                               config.get('Annovar','table_annovar'))
        self.annovar_human_db                   = os.path.join(config.get('Annovar','annovar_home'),
                                                               config.get('Annovar','annovar_human_db'))
        self.annovar_1000genomes_eur            = config.get('Annovar','annovar_1000genomes_eur')
        self.annovar_1000genomes_eur_maf_cutoff = config.get('Annovar','annovar_1000genomes_eur_maf_cutoff')
        self.annovar_inhouse_dbs                = config.get('Annovar','annovar_inhouse_dbs')
        self.omim_gene_phenotype_map_file       = config.get('Annovar','omim_gene_phenotype_map_file')
        
    #def load_settings_from_JSON(self, settings):
        #""" Update setting with JSON dictionary """
        
        #d = json.loads(settings, parse_int=int)
        #if d.has_key('target_tasks'):
            #self.target_tasks = d['target_tasks']
        #if d.has_key('num_jobs'):
            #self.num_jobs = int(d['num_jobs'])
    

    def is_runnable(self):
        """ returns true if all required settings are set """
        
        if len(target_tasks) < 1: 
            return False
        # check task names?
        
        if len(self.input_fastqs) < 2 or \
            not os.path.exists(self.runs_scratch_dir):
            return False
            
        return True


    def write_config(self, config_file):
        """ Save the config in a file """
        
        # write root paths
        
        # write reference data
        
        # write tool paths
        
        pass
        
        
