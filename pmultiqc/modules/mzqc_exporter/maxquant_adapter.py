from pathlib import Path
from typing import List

from mzqc import MZQCFile as qc

from pmultiqc.modules.mzqc_exporter import MzQCExporterModule
from pmultiqc.modules.mzqc_exporter.file_format_utils import FileFormatUtils

from pmultiqc.modules.maxquant.maxquant import MQMetaData

class MaxQuantAdapter:
    """
    Adapter class to handle MaxQuant-specific mzQC CV entry creation.
    Separates mzQC logic from the main MaxQuant module.
    """
    
    def __init__(self, mzqc_exporter : MzQCExporterModule):
        """
        Initialize the adapter with an mzQC exporter instance.
        
        Args:
            mzqc_exporter: MzQCExporterModule instance for adding metadata
        """
        self.mzqc_exporter = mzqc_exporter
    
    def process_metadata(self, combined_data : MQMetaData):
        """
        Process MaxQuant metadata and create appropriate mzQC CV entries.
        
        Args:
            combined_data: MQMetaData object containing version, FASTA paths, raw file paths
        """
        self._create_analysis_software_entry(combined_data.version)
        self._create_fasta_file_entries(combined_data.fastafile_paths)
        self._create_raw_file_entries(combined_data.rawfile_paths)
    
    def _create_analysis_software_entry(self, version : str):
        """
        Create AnalysisSoftware CV entry for MaxQuant.
        
        Args:
            combined_data: MQMetaData object containing version information
        """
        maxquant_mzqc = qc.AnalysisSoftware(
            accession="MS:1001583", 
            name="MaxQuant",
            description="MaxQuant is a quantitative proteomics software package designed for analyzing large mass spectrometric data sets. It is specifically aimed at high resolution MS data.",
            version=version, 
            uri="https://www.maxquant.org/"
        )
        self.mzqc_exporter.add_base_metadata_to_run_quality(maxquant_mzqc)
    
    def _create_fasta_file_entries(self, fastafile_paths : List[str]):
        """
        Create InputFile CV entries for FASTA files.
        
        Args:
            combined_data: list of FASTA file paths
        """
        for fastafile in fastafile_paths:
            fasta_input_file = qc.InputFile(
                name = Path(fastafile).name,
                location = fastafile,
                fileFormat = qc.CvParameter(accession="MS:1001348", name="FASTA format")
            )
            self.mzqc_exporter.add_base_metadata_to_run_quality(fasta_input_file)
    
    def _create_raw_file_entries(self, rawfile_paths : List[str]):
        """
        Create InputFile CV entries for raw files, using the filename prefix (=stem) as label.
        
        Args:
            rawfile_paths: list of raw file paths
        """
        for rawfile in rawfile_paths:
            # Get the appropriate CV term and name for the file format
            cv_accession, cv_name = FileFormatUtils.filename_to_cv(rawfile)
            
            # Create InputFile for raw file with proper format CV term
            raw_input_file = qc.InputFile(
                name=Path(rawfile).name,
                location=rawfile,
                fileFormat=qc.CvParameter(accession=cv_accession, name=cv_name)
            )
            self.mzqc_exporter.add_metadata_for_run_quality(Path(rawfile).stem, raw_input_file)
    