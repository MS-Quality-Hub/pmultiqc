import json
import os
from datetime import datetime
from pathlib import Path
from typing import Mapping, Sequence, Union

from mzqc import MZQCFile as qc

from multiqc import config
from multiqc.plots.bargraph import InputDatasetT

# define common metric units
metric_unit_count = {"accession": "UO:0000189",
                     "name": "count unit"}

class MzQCExporterModule():

    def __init__(self):
        from pmultiqc.modules.common.logging import get_logger
        self.log = get_logger(self.__class__.__module__)

        self.run_quality_base_metadata = []   # metadata, which should be applied to each runquality (e.g. inputfile from FASTA filenames)
        self.run_quality_metrics = {}   # intermediate store for lists of run_qualities per sample
        self.run_quality_metadata = {}  # intermediate store for metadata per sample: each has a dict 
                                        # with the keys 'input_files' and 'analysis_software', values are lists of these

 
    def create_export(self):
        """
        Create the mzQC object(s) and write the mzQC file
        """

        self.log.info("Starting mzQC export ...")

        cv_ms = qc.ControlledVocabulary(
            name="Proteomics Standards Initiative Mass Spectrometry Ontology",
            version="4.1.212",
            uri="https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo")

        mzqc_data = {
            "run_qualities" : [],
            "set_qualities": [],
        }

        # add all base metadata to each RunQualities' metadata
        for label in self.run_quality_metrics.keys():
            for metadata in self.run_quality_base_metadata:
                self.add_metadata_for_run_quality(label, metadata)

        # create quality metrics and metadata per RunQuality
        for label, qualitymetrics in self.run_quality_metrics.items():
            meta = qc.MetaDataParameters(label = label,
                                         inputFiles = self.run_quality_metadata[label]['input_files'],
                                         analysisSoftware = self.run_quality_metadata[label]['analysis_software'])
            rq = qc.RunQuality(metadata = meta, qualityMetrics = self.run_quality_metrics[label])
            mzqc_data['run_qualities'].append(rq)

        # create the mzQC object storing all the data
        mzqc = qc.MzQcFile(version = "1.0.0",
                           creationDate = datetime.now().isoformat(),
                           runQualities = mzqc_data["run_qualities"],
                           setQualities = mzqc_data["set_qualities"],
                           controlledVocabularies = [cv_ms])
        
        # write the mzQC file
        output_dir = Path(config.output_dir) if config.output_dir is not None else Path("./")  ## use current dir
        mzqc_filename = os.path.join(output_dir, "pmultiqc.mzqc")
        with open(mzqc_filename, "w") as mzqc_file:
            mzqc_file.write(json.dumps(json.loads(qc.JsonSerialisable.to_json(mzqc)), indent = 2))
        self.log.info(f"Done exporting mzQC to {mzqc_filename}")


    def _create_empty_run(self, label: str):
        # initialize metadata for RunQuality
        self.run_quality_metrics[label] = []

        # initialize metadata store for RunQuality
        self.run_quality_metadata[label] = {
            'input_files': [],
            'analysis_software': [],
        }
    
    
    def add_metric_to_run_quality(self, label: str, qm: qc.QualityMetric):
        if label not in self.run_quality_metrics.keys():
            self._create_empty_run(label)
        self.run_quality_metrics[label].append(qm)


    def add_base_metadata_to_run_quality(self, metadata: qc.InputFile | qc.AnalysisSoftware):
        """
        This function adds information which should be used for all runs in the mzQC export, like
        the used AnalysisSoftware (if equal), same FASTAs etc.
        """
        # the infromation is only collected at this step, and added later during creation of the mzQC object
        self.run_quality_base_metadata.append(metadata)
    

    def add_metadata_for_run_quality(self, label: str, metadata: qc.InputFile | qc.AnalysisSoftware):
        """ 
        Insert into the runquality with key 'label':
         - InputFile into input_files and
         - AnalysisSoftware into analysis_software
        """
        if label not in self.run_quality_metrics.keys():
            self._create_empty_run(label)
        
        if isinstance(metadata, qc.InputFile):
            self.run_quality_metadata[label]['input_files'].append(metadata)
        elif isinstance(metadata, qc.AnalysisSoftware):
            self.run_quality_metadata[label]['analysis_software'].append(metadata)
    

    def _add_count_metric_for_run_qualities(self, accession: str, name: str, data: Union[InputDatasetT, Sequence[InputDatasetT]]):
        # data should be a mapping from "file name" to "categories -> values"
        for label, sample_data in data.items():
            metric_count = None
            if isinstance(sample_data, Mapping):
                # add up the counts for all categories
                metric_count = sum(x for x in sample_data.values())
            else:
                # try to cast to int
                metric_count = int(sample_data)

            # create a QualityMetric for this sample
            qm = qc.QualityMetric(accession=accession,
                              name=name,
                              value=metric_count,     
                              unit=metric_unit_count)
            
            self.add_metric_to_run_quality(label=label, qm=qm)

    def add_peptide_id_count(self, data: Union[InputDatasetT, Sequence[InputDatasetT]]):
        """
        Add peptide identification count metrics to mzQC export.
        
        Args:
            data: Peptide count data mapping from sample names to counts
        """
        self._add_count_metric_for_run_qualities("MS:1003250", "count of identified peptidoforms", data)

    def add_protein_id_count(self, data: Union[InputDatasetT, Sequence[InputDatasetT]]):
        """
        Add protein identification count metrics to mzQC export.
        
        Args:
            data: Protein count data mapping from sample names to counts
        """
        self._add_count_metric_for_run_qualities("MS:1002404", "count of identified proteins", data)
