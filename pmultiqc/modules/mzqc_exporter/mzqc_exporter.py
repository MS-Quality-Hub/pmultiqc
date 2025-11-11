import json
import os
from datetime import datetime
from pathlib import Path
from typing import Mapping, Sequence, Union

from mzqc import MZQCFile as qc

from multiqc import config
from multiqc.plots.bargraph import InputDatasetT

# define common metric units
metric_unit_count = {"unit_accession": "UO:0000189",
                     "unit_name": "count unit"}

class MzQCExporterModule():

    def __init__(self):
        from pmultiqc.modules.common.logging import get_logger
        self.log = get_logger(self.__class__.__module__)

        self.base_metadata = []         # metadata, which should be applied to each run
        self.run_quality_metrics = {}   # intermediate store for lists of run_qualities per run
        self.run_metadata = {}          # intermediate store for metadata per run, each run has a dict 
                                        # with the keys 'input_files' and 'analysis_software', values are lists of these

 
    def create_export(self):
        """
        Finally, create the mzQC object(s) and write the file
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

        # add all base metadata to each sample's metadata
        for sample_name in self.run_quality_metrics.keys():
            for metadata in self.base_metadata:
                self.add_metadata_for_run(sample_name, metadata)
        
        self.check_for_input_files()

        # create quality metrics and metadata per sample
        for sample_name, qualitymetrics in self.run_quality_metrics.items():
            meta = qc.MetaDataParameters(label = sample_name,
                                         inputFiles = self.run_metadata[sample_name]['input_files'],
                                         analysisSoftware = self.run_metadata[sample_name]['analysis_software'])
            rq = qc.RunQuality(metadata = meta, qualityMetrics = self.run_quality_metrics[sample_name])
            mzqc_data['run_qualities'].append(rq)

        # create the mzQC object storing all the data
        mzqc = qc.MzQcFile(version = "1.0.0",
                           creationDate = datetime.now().isoformat(),
                           runQualities = mzqc_data["run_qualities"],
                           setQualities = mzqc_data["set_qualities"],
                           controlledVocabularies = [cv_ms])
        
        # write out the mzQC file
        output_dir = Path(config.output_dir) if config.output_dir is not None else None
        mzqc_filename = None
        if output_dir is not None:
            # TODO: set the pmultiqc-output by parameters
            mzqc_filename = os.path.join(output_dir, "pmultiqc.mzqc")
            with open(mzqc_filename, "w") as mzqc_file:
                mzqc_file.write(json.dumps(json.loads(qc.JsonSerialisable.to_json(mzqc)), indent = 2))
        
        self.log.info(f"Done exporting mzQC to {mzqc_filename}")


    def check_for_input_files(self):
        """
        Checks whether there are input files annotated for the sample_names.
        If not, create stubs which need to be fixed later
        """

        for sample_name in self.run_metadata.keys():
            if len(self.run_metadata[sample_name]['input_files']) < 1:
                self.log.warning(f"No input file given for {sample_name}, creating stub")
                
                input_file_stub = qc.InputFile(name=sample_name,
                                               location="UNKNOWN", 
                                               fileFormat=None, 
                                               fileProperties=[])
                
                self.run_metadata[sample_name]['input_files'].append(input_file_stub)



    def create_run(self, run_id: str):
        # initialize metadata for run
        self.run_quality_metrics[run_id] = []

        # initialize metadata store for run
        self.run_metadata[run_id] = {
            'input_files': [],
            'analysis_software': [],
        }
    
    
    def add_metric_to_run(self, sample_name: str, qm: qc.QualityMetric):
        if sample_name not in self.run_quality_metrics.keys():
            self.create_run(sample_name)
        self.run_quality_metrics[sample_name].append(qm)


    def add_base_metadata(self, metadata: qc.InputFile | qc.AnalysisSoftware):
        """
        This function adds information which should be used for all runs in the mzQC export, like
        the used AnalysisSoftware (if equal), same FASTAs etc.
        """
        # the infromation is only collected at this step, and added later during creation of the mzQC object
        self.base_metadata.append(metadata)
    

    def add_metadata_for_run(self, sample_name: str, metadata: qc.InputFile | qc.AnalysisSoftware):
        """ 
        Insert
         - InputFile into input_files and
         - AnalysisSoftware into analysis_software
        """
        if sample_name not in self.run_quality_metrics.keys():
            self.create_run(sample_name)
        
        if isinstance(metadata, qc.InputFile):
            self.run_metadata[sample_name]['input_files'].append(metadata)
        elif isinstance(metadata, qc.AnalysisSoftware):
            self.run_metadata[sample_name]['analysis_software'].append(metadata)
    

    def add_metric(self,
                   data: Union[InputDatasetT, Sequence[InputDatasetT]],
                   accession: str = None,
                   ) -> bool | None:
        if accession == "MS:1002404":
            self.add_count_metric_per_sample("MS:1002404", "count of identified proteins", data)
        elif accession == "MS:1003250":
            self.add_count_metric_per_sample("MS:1003250", "count of identified peptidoforms", data)
        
        return True
    

    def add_count_metric_per_sample(self, accession: str, name: str, data: Union[InputDatasetT, Sequence[InputDatasetT]]) -> bool | None:
        # data should be a mapping from "file name" to "categories -> values"
        for sample_name, sample_data in data.items():
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
            
            self.add_metric_to_run(sample_name=sample_name, qm=qm)
