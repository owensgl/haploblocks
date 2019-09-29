"""
  Define the objects and parameters involved in the SNP-calling experiment
"""
import io
import json
import os, os.path

from collections import OrderedDict
from snakemake.logging import logger

from helpers import Contig

# singleton ExperimentParams
Experiment = None

class ExperimentParams(object):
    """Generates unique identifiers based on the parameters of the experiment"""

    VC_SINGULARITY_IMAGE = "docker://broadinstitute/gatk:4.0.1.2"
    VCF_SINGULARITY_IMAGE = "docker://broadinstitute/gatk:4.0.1.2"
    VARIANT_MODEL_SINGULARITY_IMAGE = "docker://broadinstitute/gatk:4.0.6.0"
    BEAGLE_SINGULARITY_IMAGE = "docker://rieseberglab/analytics:2"
    PLINK_SINGULARITY_IMAGE = "docker://rieseberglab/analytics:2"

    GENOTYPE_GVCFS_ARGS = " ".join([
            "--only-output-calls-starting-in-intervals",
            "--use-new-qual-calculator"
            ])

    def __init__(self, config):
        self._vcf_call_uid = None
        self._vcf_params = None
        self._goldset_uid = None
        self._gold_params = None
        self._vcf_filter_params = None
        self._vcf_filter_uid = None
        self.config = config

    @classmethod
    def from_config(cls, config):
        global Experiment
        exp = cls(config)
        Experiment = exp
        return exp

    @property
    def sample_slice_uid(self):
        return self.config['sample_slice_uid']

    @property
    def contig_list_uid(self):
        return self.config['contig_list_uid']

    @property
    def vcf_call_uid(self):
        """return a UID that uniquely describes a call for SNPs.
           the UID should cover all inputs and parameters used in the
           invocation.
        """

        if self._vcf_call_uid:
            return self._vcf_call_uid

        sample_uid = self.sample_slice_uid

        # FIXME fold in SHA1 for all input files

        self._vcf_params = {
            "sample_slice_uid": sample_uid,
            "singularity_img": self.VCF_SINGULARITY_IMAGE,
            "params": { 'gvcf_args': Experiment.GENOTYPE_GVCFS_ARGS }
        }
        with io.BytesIO() as datafile:
            datafile.write(json.dumps(self._vcf_params, sort_keys=True, separators=(",",  ":")).encode('utf-8'))
            datafile.seek(0)
            self._vcf_call_uid = self.config['paramsdb'].put_stream(datafile)
        logger.run_info("Experiment.vcf_call_uid: %s  (%s)" % (self._vcf_call_uid, json.dumps(self._vcf_params)))
        return self._vcf_call_uid

    @property
    def goldset_uid(self):
        """
        return a UID that uniquely describes the current call for a gold set.
        """
        if self._goldset_uid:
            return self._goldset_uid

        # remove empty lines and comments
        filter_data = []
        with open(self.config['goldset_filter'], "r") as filterfd:
            # remove comments
            for line in filterfd:
                line = line.strip()
                if not line or line.startswith("#"): continue
                filter_data.append(line)
        filter_text = "\n".join(filter_data)

        logger.run_info("Filter text to apply:\n    %s" % (filter_text.replace("\n", "\n    ")),)

        self._gold_params = {
            "vcf_call_uid": self.vcf_call_uid,
            "sample_slice_uid": self.sample_slice_uid,
            "contig_list_uid": self.contig_list_uid,
            "goldset_filter_text": filter_text
        }

        with io.BytesIO() as datafile:
            datafile.write(json.dumps(self._gold_params, sort_keys=True, separators=(",",  ":")).encode('utf-8'))
            datafile.seek(0)
            self._goldset_uid = self.config['paramsdb'].put_stream(datafile)
        return self._goldset_uid

    @property
    def tranchevalues(self):
        tranche_values = [ float(x) for x in self.config['recalibration_tranches'] ]
        if [ x for x in tranche_values if x < 0 ]:
            raise Exception("Tranche values should be positive.")
        as_string = [ "%3.2f" % (x,) for x in sorted(tranche_values) ]
        return "-".join(as_string)

    def vcf_sample_slice(self, vcf_call_uid):
        my_call_uid = self.vcf_call_uid
        if vcf_call_uid != my_call_uid:
            raise Exception("the sample slice for vcf_call_uid {} is unknown".format(vcf_call_uid))

        sample_slice_uid = self._vcf_params['sample_slice_uid']
        return get_sample_slice_by_uid(sample_slice_uid)

    @property
    def vcf_filter_uid(self):
        if self.config.get('vcf_filter_uid'):
            return self.config.get('vcf_filter_uid')

        # construct filter from known settings

        if "vcf_filter" not in self.config:
            raise Exception("Configuration key `vcf_filter` must indicate the file path to the filter to apply")

        exclude_chrom_prefixes = self.config.get('vcf_filter_exclude_chrom_prefixes', [])
        if not isinstance(exclude_chrom_prefixes, list):
            raise Exception("Parameter exclude_chrom_prefix must be a list of strings")

        params = {
            'vcf_call_uid': self.vcf_call_uid,
            'contig_list_uid': self.contig_list_uid,
            'tranchevalues': self.tranchevalues,
            'exclude_chrom_prefixes': exclude_chrom_prefixes
        }

        filter_obj = VcfFilter(**params)
        filter_obj.load_gatk_args(self.config['vcf_filter'])
        self.config['vcf_filter_uid'] = filter_obj.filter_uid

        # circular -- will hit cache
        return self.vcf_filter_uid

    @property
    def beagle_call_uid(self):
        if self.config.get('beagle_call_uid'):
            return self.config.get('beagle_call_uid')

        if not self.config.get('recombination_map'):
            raise Exception("The recombination_map entry must be set in the config")

        params = {
            'recombination_map': self.config.get('recombination_map'),
            'vcf_filter_uid': self.vcf_filter_uid,
            'beagle_args': self.config.get('beagle_args', [])
        }

        beagle_obj = BeagleCall(**params)
        self.config['beagle_call_uid'] = beagle_obj.beagle_uid

        # circular -- will hit cache
        return self.beagle_call_uid

    @property
    def plink_call_uid(self):
        if self.config.get('plink_call_uid'):
            return self.config.get('plink_call_uid')

        params = {
            'beagle_call_uid': self.beagle_call_uid,
            'plink_args': self.config.get('plink_args', [])
        }
        plink_obj = PlinkCall(**params)
        self.config['plink_call_uid'] = plink_uid

        # circular -- will hit cache
        return self.plink_call_uid

class BeagleCall(object):
    """represents the result of a beagle call on a filtered set"""
    cache = {}
    def __init__(self, **kwargs):
        self.params = {
            'singularity_img': Experiment.BEAGLE_SINGULARITY_IMAGE,
            'beagle_url': "https://faculty.washington.edu/browning/beagle/beagle.10Jun18.811.jar",
            'beagle_sha1': "0a3234a77f455c31c5084b38ef7e0b5b8d208723",
            'recombination_map': kwargs['recombination_map'],
            'vcf_filter_uid': kwargs['vcf_filter_uid'],
            'beagle_args': [x for x in kwargs.get('beagle_args', [])]
        }
        self._beagle_uid = None

    @classmethod
    def from_dict(cls, obj):
        return cls(**obj)

    @classmethod
    def from_beagle_uid(cls, beagle_uid):
        """loads the filter information from the given uid"""
        cached = cls.cache.get(beagle_uid)
        if cached:
            return cached

        db = Experiment.config['paramsdb']
        with db.get(beagle_uid) as fil:
            obj = cls.from_dict(json.load(fil))

        cls.cache[beagle_uid] = obj
        return obj

    def to_params(self):
        return json.dumps(self.params, sort_keys=True, separators=(",", ":"))

    @property
    def beagle_uid(self):
        if self._beagle_uid:
            return self._beagle_uid

        with io.BytesIO() as datafile:
            datafile.write(json.dumps(self.params, sort_keys=True, separators=(",",  ":")).encode('utf-8'))
            datafile.seek(0)
            self._beagle_uid = Experiment.config['paramsdb'].put_stream(datafile)

        BeagleCall.cache[self._beagle_uid] = self
        return self._beagle_uid

class PlinkCall(object):
    """convert beagle-phased vcf to tped format for gwas"""
    cache = {}
    def __init__(self, **kwargs):
        self.params = {
            'singularity_img': Experiment.PLINK_SINGULARITY_IMAGE,
            'plink_args': [x for x in kwargs.get('plink_args', [])],
            'beagle_call_uid': kwargs['beagle_call_uid']
        }
        self._plink_uid = None

    @classmethod
    def from_dict(cls, obj):
        return cls(**obj)

    @classmethod
    def from_plink_uid(cls, plink_uid):
        """loads the filter information from the given uid"""
        cached = cls.cache.get(plink_uid)
        if cached:
            return cached

        db = Experiment.config['paramsdb']
        with db.get(plink_uid) as fil:
            obj = cls.from_dict(json.load(fil))

        cls.cache[plink_uid] = obj
        return obj

    @property
    def plink_uid(self):
        if self._plink_uid:
            return self._plink_uid

        with io.BytesIO() as datafile:
            datafile.write(json.dumps(self.params, sort_keys=True, separators=(",",  ":")).encode('utf-8'))
            datafile.seek(0)
            self._plink_uid = Experiment.config['paramsdb'].put_stream(datafile)

        PlinkCall.cache[self._plink_uid] = self
        return self._plink_uid

class VcfFilter(object):
    """represents the operation of applying a filter on a VCF call set"""

    cache = {}
    def __init__(self, **kwargs):
        # loaded/saved from file
        self.params = {
            # missing a goldsetuid
            'vcf_call_uid': kwargs['vcf_call_uid'],
            'contig_list_uid': kwargs['contig_list_uid'],
            'tranchevalues': kwargs['tranchevalues'],
            'filter_name': kwargs.get('filter_name'),
            'gatk_args': kwargs.get('gatk_args', []),
            'exclude_chrom_prefixes': kwargs.get('exclude_chrom_prefixes', [])
        }
        self._sample_slice = None
        self._filter_uid = None

    @property
    def sample_slice(self):
        if not self._sample_slice:
            self._sample_slice = Experiment.vcf_sample_slice(self.params['vcf_call_uid'])
        return self._sample_slice

    def load_gatk_args(self, filename):
        logger.run_info("Loading filter %s" % (filename,))
        filter_name = os.path.basename(filename)
        if filter_name.endswith(".py"):
            filter_name = filter_name[:-3]

        with open(filename, "r") as filterfd:
            source = filterfd.read()
        binary = compile(source, filename, "exec")
        result = {}
        exec(binary, result)
        filter_func = result.get('get_filter', None)

        if not filter_func:
            raise Exception("filter %s is missing a get_filter() function" % (filename,))

        filter_params = {
            'MAX_AN': 2*len(self.sample_slice),
        }

        self.params['filter_name'] = filter_name
        self.params['gatk_args'] = list(filter_func(filter_params))

    @property
    def filter_uid(self):
        if self._filter_uid:
            return self._filter_uid

        if not self.params['gatk_args']:
            raise Exception("Filter not loaded.")

        with io.BytesIO() as datafile:
            datafile.write(json.dumps(self.params, sort_keys=True, separators=(",",  ":")).encode('utf-8'))
            datafile.seek(0)
            self._filter_uid = Experiment.config['paramsdb'].put_stream(datafile)

        VcfFilter.cache[self._filter_uid] = self
        return self._filter_uid

    @classmethod
    def from_dict(cls, obj):
        return cls(**obj)

    @classmethod
    def from_filter_uid(cls, filter_uid):
        """loads the filter information from the given uid"""
        cached = cls.cache.get(filter_uid)
        if cached:
            return cached

        db = Experiment.config['paramsdb']
        with db.get(filter_uid) as fil:
            obj = cls.from_dict(json.load(fil))

        cls.cache[filter_uid] = obj
        return obj


def get_gold_params_by_uid(goldset_uid):
    cached = get_gold_params_by_uid.cache.get(goldset_uid)
    if cached:
        return cached

    db = Experiment.config['paramsdb']
    sample_slice = OrderedDict()

    with db.get(goldset_uid) as inputfd:
        obj = json.load(inputfd)

    get_gold_params_by_uid.cache[goldset_uid] = obj
    return obj
get_gold_params_by_uid.cache = {}

def get_sample_slice_by_uid(sample_slice_uid):
    """gets an ordereddict of sample definitions matching the sample slice uid specification given"""

    cached = get_sample_slice_by_uid.cache.get(sample_slice_uid)
    if cached:
        return cached

    db = Experiment.config['paramsdb']
    sample_slice = OrderedDict()

    with db.get(sample_slice_uid) as sample_list:
        for line in sample_list:
            line = line.strip()
            if not line: continue
            record = json.loads(line)
            sample_slice[record['name']] = record['val']

    get_sample_slice_by_uid.cache[sample_slice_uid] = sample_slice
    return sample_slice

get_sample_slice_by_uid.cache = {}

def get_sample_slice_by_name(all_samples, names):
    """
    gets an ordereddict of sample definitions matching names listed in `names`

    all_samples is a list of all samples. each entry of the list is a list of source files for that sample (one per lane).
    """
    sample_slice = OrderedDict()
    selection = OrderedDict()

    for name in names:
        selection[name] = True

    for sample_lanes in all_samples:
        # we use [0] but all lanes have the same sample name.
        for first_key in sample_lanes.keys():
            break
        sample_name = sample_lanes[first_key]['name']
        if sample_name in selection:
            sample_slice[sample_name] = sample_lanes
    return sample_slice

def get_sample_slice_by_number(all_samples, samplerange):
    """
    gets an ordereddict of sample definitions matching a range specification in
    number syntax.

    all_samples is a list of all samples. each entry of the list is a list of source files for that sample (one per lane).
    """
    sample_slice = OrderedDict()

    selected_indices = {}

    for subrange in samplerange.split(","):
        start, end = subrange.split("-", 1) if "-" in subrange else (subrange, subrange)
        if start == "": start = "1"
        if end == "": end = str(len(all_samples))
        start, end = int(start, 10), int(end, 10)
        for sample_i in range(max(0, start - 1), min(end, len(all_samples))):
            selected_indices[sample_i] = True

    for sample_i in sorted(selected_indices.keys()):
        sample_lanes = all_samples[sample_i]

        # we use [0] but all lanes have the same sample name.
        for first_key in sample_lanes.keys():
            break
        sample_name = sample_lanes[first_key]['name']
        sample_slice[sample_name] = sample_lanes

    return sample_slice

def _contigs_from_file(contigfd):
    """load a list of contigs from plain text file,
       oone contig by line

       file syntax is:
          CHROM POS END LEN\n
    """
    def _valid_contig(line):
        line = line.strip()
        if not line: return False
        if line.startswith("#") or line.startswith("/"): return False
        return True

    def _parse_window(line):
        toks = line.strip().split()
        return Contig(toks[0], int(toks[1], 10), int(toks[2], 10), int(toks[3], 10))

    return [ _parse_window(line) for line in contigfd if _valid_contig(line) ]


def get_contig_slice_by_uid(contig_list_uid):
    """gets a list of Contigs based on the contig_list identifier"""
    cached = get_contig_slice_by_uid.cache.get(contig_list_uid)
    if cached:
        return cached

    db = Experiment.config['paramsdb']
    sample_slice = OrderedDict()
    with db.get(Experiment.config['contig_list_uid']) as contigfd:
        contig_slice = _contigs_from_file(contigfd)

    get_contig_slice_by_uid.cache[contig_list_uid] = contig_slice
    return contig_slice

get_contig_slice_by_uid.cache = {}

def get_contig_slice_by_number(all_contigs, inputrange):
    """gets an ordereddict of sample definitions matching a range specification"""

    selected_indices = {}

    # we only select each value once.
    # e.g. 1,2,1-3,2-4 == 1-4
    for subrange in inputrange.split(","):
        start, end = subrange.split("-", 1) if "-" in subrange else (subrange, subrange)
        if start == "": start = "1"
        if end == "": end = str(len(all_contigs))
        start, end = int(start, 10), int(end, 10)
        for sample_i in range(start - 1, end):
            selected_indices[sample_i] = True

    contig_slice = [ all_contigs[i] for i in sorted(selected_indices.keys())
                     if i >= 0 and i < len(all_contigs) ]
    return contig_slice


