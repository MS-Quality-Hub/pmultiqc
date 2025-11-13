import os
from pathlib import Path
from typing import List, Optional
import pandas as pd

from mzqc import MZQCFile as qc

from pmultiqc.modules.mzqc_exporter import MzQCExporterModule
from pmultiqc.modules.mzqc_exporter.file_format_utils import FileFormatUtils


class MzIdentMLAdapter:
    """
    Adapter class to handle mzIdentML-specific mzQC CV entry creation.
    Separates mzQC logic from the main mzIdentML module.
    """
    
    def __init__(self, mzqc_exporter: MzQCExporterModule):
        """
        Initialize the adapter with an mzQC exporter instance.
        
        Args:
            mzqc_exporter: MzQCExporterModule instance for adding metadata
        """
        self.mzqc_exporter = mzqc_exporter
    
    def process_metadata(self, mzml_ms_df: pd.DataFrame, ms_paths: List[str]):
        """
        Process mzIdentML metadata and create appropriate mzQC CV entries.
        
        Args:
            mzml_ms_df: DataFrame containing mzML MS data with filename column
            ms_paths: List of paths to MS files
        """
        self._create_input_file_entries(mzml_ms_df, ms_paths)
    
    def _create_input_file_entries(self, mzml_ms_df: pd.DataFrame, ms_paths: List[str]):
        """
        Create InputFile CV entries for MS files based on sample names from mzML data.
        
        Args:
            mzml_ms_df: DataFrame containing mzML MS data with filename column
            ms_paths: List of paths to MS files
        """
        if mzml_ms_df is None or mzml_ms_df.empty:
            return
            
        # Get unique sample names from the mzML data
        for sample_name in set(mzml_ms_df["filename"].unique()):
            sample_path = self._find_sample_path(sample_name, ms_paths)
            
            # Get the appropriate CV term and name for the file format
            cv_accession, cv_name = FileFormatUtils.filename_to_cv(sample_path or sample_name)
            
            # Create InputFile for MS file with proper format CV term
            input_file = qc.InputFile(
                name=sample_name,
                location=sample_path,
                fileFormat=qc.CvParameter(accession=cv_accession, name=cv_name),
                fileProperties=[]
            )
            
            self.mzqc_exporter.add_metadata_for_run_quality(sample_name, input_file)
    
    def _find_sample_path(self, sample_name: str, ms_paths: List[str]) -> Optional[str]:
        """
        Find the full path for a sample based on its name.
        
        Args:
            sample_name: Name of the sample to find
            ms_paths: List of paths to search through
            
        Returns:
            Full path to the sample file if found, None otherwise
        """
        for file_path in ms_paths:
            if os.path.basename(file_path).startswith(sample_name):
                return file_path
        return None
    