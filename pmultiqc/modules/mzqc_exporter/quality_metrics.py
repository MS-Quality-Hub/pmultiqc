"""
Quality metrics module for mzQC export functionality.

This module contains specific quality metric implementations that can be added to mzQC files.
"""

from typing import Mapping, Sequence, Union

from mzqc import MZQCFile as qc
from multiqc.plots.bargraph import InputDatasetT

# define common metric units
METRIC_UNIT_COUNT = {
    "unit_accession": "UO:0000189",
    "unit_name": "count unit"
}


class QualityMetrics:
    """
    Class containing quality metric implementations for mzQC export.
    
    This class provides methods to create specific quality metrics that can be
    added to mzQC files. Each method handles the conversion of data to proper
    mzQC QualityMetric objects.
    """
    
    def __init__(self, mzqc_exporter):
        """
        Initialize with reference to the mzQC exporter.
        
        Args:
            mzqc_exporter: MzQCExporterModule instance for adding metrics
        """
        self.mzqc_exporter = mzqc_exporter
    
    def add_peptide_id_count(self, data: Union[InputDatasetT, Sequence[InputDatasetT]]):
        """
        Add peptide identification count metrics to mzQC export.
        
        Args:
            data: Peptide count data mapping from sample names to counts
        """
        self._add_count_metric_for_run_qualities(
            "MS:1003250", 
            "count of identified peptidoforms", 
            data
        )
    
    def add_protein_id_count(self, data: Union[InputDatasetT, Sequence[InputDatasetT]]):
        """
        Add protein identification count metrics to mzQC export.
        
        Args:
            data: Protein count data mapping from sample names to counts
        """
        self._add_count_metric_for_run_qualities(
            "MS:1002404", 
            "count of identified proteins", 
            data
        )
    
    def _add_count_metric_for_run_qualities(self, accession: str, name: str, data: Union[InputDatasetT, Sequence[InputDatasetT]]):
        """
        Helper method to add count-based metrics to run qualities.
        
        Args:
            accession: PSI-MS CV accession for the metric
            name: Human-readable name for the metric
            data: Data mapping from sample names to counts
        """
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
            qm = qc.QualityMetric(
                accession=accession,
                name=name,
                value=metric_count,     
                unit=METRIC_UNIT_COUNT
            )
            
            self.mzqc_exporter.add_metric_to_run_quality(label=label, qm=qm)