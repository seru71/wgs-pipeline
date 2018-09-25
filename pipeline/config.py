#
#
#
# Pipeline configurator
#

import os, json
import StringIO
import sys
import logging
import logging.handlers



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



class PipelineConfig:


    # Attributes:
    #
    #  - reference_root
    #  - scratch_root
    #  - results_archive
    #  - fastq_archive
    #  - tmp_dir
    


    # Tool paths
    # - docker_bin
    

    def __init__(self):
        
        self.drmaa_session = None
        
        self.dockerize = True
        
        self.logger = None
        
        self.bcl2fastq   = None
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
        self.annovar_1000genomes_eur_nad_cutoff = None
        self.annovar_inhouse_dbs                = None
        self.omim_gene_phenotype_map            = None



        # run settings
        self.reference_root = None
        self.scratch_root = None
        self.run_folder   = None
#        self.input_fastqs = None
        self.run_id       = None
        self.runs_scratch_dir = None
        self.tmp_dir      = None

        self.adapters  = None
        self.reference = None

        self.target_tasks    = []
        self.log_file        = None
        self.num_jobs        = 1
        self.verbosity_level = 0
        self.dry_run         = False
        self.rebuild_mode    = False
        self.run_on_bcl_tile = None


    def set_logger(self, logger):
        self.logger = logger
   
    def set_runfolder(self, runfolder):
        if runfolder != None and \
            os.path.exists(runfolder) and \
            os.path.exists(os.path.join(runfolder,'SampleSheet.csv')):

            self.run_folder = runfolder
            
        else:
            raise Exception("Incorrect runfolder\'s path [%s] or missing SampleSheet file." % runfolder)

    def set_num_jobs(self,jobs):
        if jobs is not None:
            self.num_jobs = jobs


    def set_input_fastqs(self, fastqs):
        self.input_fastqs = fastqs
    
    
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

        import ConfigParser
        config = ConfigParser.ConfigParser()
        config.read(cfg_file)
        
        
        self.reference_root = config.get('Paths','reference-root')
        
        self.scratch_root = os.getcwd()
        try:
            self.scratch_root   = config.get('Paths','scratch-root')
        except ConfigParser.NoOptionError:
            self.logger.info('Scratch-root setting is missing. Using current directory: %s' % self.scratch_root)


        if (self.run_folder != None):
            self.run_id = os.path.basename(self.run_folder)
        else:
            raise Exception('Set runfolder with PipelineConfig.set_runfolder() before loading settings')
                  
        
        #
        # TODO
        # needs to be updated on update of settings
        #
        self.runs_scratch_dir = os.path.join(self.scratch_root, self.run_id) if self.run_folder != None else self.scratch_root
        self.logger.info('Run\'s scratch directory: %s' % self.runs_scratch_dir)
          
        # optional results and fastq archive dirs  
        self.results_archive = None
        try:
            self.results_archive = config.get('Paths','results-archive')
        except ConfigParser.NoOptionError:
            self.logger.info('No results-archive provided. Results will not be archived outside of the run\'s scratch directory.')
        
        self.fastq_archive = None
        try:
            self.fastq_archive = config.get('Paths','fastq-archive')
        except ConfigParser.NoOptionError:
            self.logger.info('No fastq-archive provided. Fastq files will not be archived outside of the run\'s scratch directory.')
    
        
        # optional /tmp dir
        self.tmp_dir = '/tmp'
        try:
            self.tmp_dir = config.get('Paths','tmp-dir')
        except ConfigParser.NoOptionError:
            self.logger.info('No tmp-dir provided. /tmp will be used.')
  
  
                  
            
        # reference files
        self.reference = os.path.join(self.reference_root, config.get('Resources','reference-genome'))
        self.capture = os.path.join(self.reference_root, config.get('Resources','capture-regions-bed'))
        self.capture_qualimap = os.path.join(self.reference_root, config.get('Resources','capture-regions-bed-for-qualimap'))
        self.capture_plus = os.path.join(self.reference_root, config.get('Resources', 'capture-plus-regions-bed'))
        self.gene_coordinates = os.path.join(self.reference_root, config.get('Resources', 'gene-coordinates'))
        
        self.adapters = os.path.join(self.reference_root, config.get('Resources', 'adapters-fasta'))
        
        # tools
        self.bcl2fastq   = config.get('Tools','bcl2fastq')
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
        self.convert_to_annovar                 = config.get('Annovar','convert_to_annovar')
        self.annovar_annotate                   = config.get('Annovar','annovar_annotate')
        self.table_annovar                      = config.get('Annovar','table_annovar')
        self.annovar_human_db                   = config.get('Annovar','annovar_human_db')
        self.annovar_1000genomes_eur            = config.get('Annovar','annovar_1000genomes_eur')
        self.annovar_1000genomes_eur_nad_cutoff = config.get('Annovar','annovar_1000genomes_eur_nad_cutoff')
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
        
        if self.run_folder is None or \
            not os.path.exists(self.run_folder) or \
            not os.path.exists(os.path.join(self.run_folder, self.run_id, 'SampleSheet.csv')):
            return False
            
        return True


    def write_config(self, config_file):
        """ Save the config in a file """
        
        # write root paths
        
        # write reference data
        
        # write tool paths
        
        pass
        
        
