from pathlib import Path
from typing import Optional

from mzqc import MZQCFile as qc


class MaxQuantAdapter:
    """
    Adapter class to handle MaxQuant-specific mzQC CV entry creation.
    Separates mzQC logic from the main MaxQuant module.
    """
    
    def __init__(self, mzqc_exporter):
        """
        Initialize the adapter with an mzQC exporter instance.
        
        Args:
            mzqc_exporter: MzQCExporterModule instance for adding metadata
        """
        self.mzqc_exporter = mzqc_exporter
    
    def process_metadata(self, combined_data):
        """
        Process MaxQuant metadata and create appropriate mzQC CV entries.
        
        Args:
            combined_data: MQMetaData object containing version, FASTA paths, raw file paths
        """
        self._create_analysis_software_entry(combined_data)
        self._create_fasta_file_entries(combined_data)
        self._create_raw_file_entries(combined_data)
    
    def _create_analysis_software_entry(self, combined_data):
        """
        Create AnalysisSoftware CV entry for MaxQuant.
        
        Args:
            combined_data: MQMetaData object containing version information
        """
        maxquant_mzqc = qc.AnalysisSoftware(
            accession="MS:1001583", 
            name="MaxQuant",
            description="MaxQuant is a quantitative proteomics software package designed for analyzing large mass spectrometric data sets. It is specifically aimed at high resolution MS data.",
            version=combined_data.version, 
            uri="https://www.maxquant.org/"
        )
        self.mzqc_exporter.add_base_metadata(maxquant_mzqc)
    
    def _create_fasta_file_entries(self, combined_data):
        """
        Create InputFile CV entries for FASTA files.
        
        Args:
            combined_data: MQMetaData object containing FASTA file paths
        """
        for fastafile in combined_data.fastafile_paths:
            fasta_input_file = qc.InputFile(
                name = Path(fastafile).name,
                location = fastafile,
                fileFormat = qc.CvParameter(accession="MS:1001348", name="FASTA format")
            )
            self.mzqc_exporter.add_base_metadata(fasta_input_file)
    
    def _create_raw_file_entries(self, combined_data):
        """
        Create InputFile CV entries for raw files.
        
        Args:
            combined_data: MQMetaData object containing raw file paths
        """
        for rawfile in combined_data.rawfile_paths:
            # Get the appropriate CV term and name for the file format
            cv_accession, cv_name = self._filename_to_cv(rawfile)
            
            # Create InputFile for raw file with proper format CV term
            raw_input_file = qc.InputFile(
                name=Path(rawfile).name,
                location=rawfile,
                fileFormat=qc.CvParameter(accession=cv_accession, name=cv_name)
            )
            self.mzqc_exporter.add_base_metadata(raw_input_file)
    
    def _filename_to_cv(self, filepath):
        """
        For a given filename (e.g. "test.mzML"), check the suffix and translate it to a PSI-MS CV term.
        
        The following mapping is currently known:
        .raw    : MS:1000563 ! Thermo RAW format
        .mzML   : MS:1000584 ! mzML format
        .mzData : MS:1000564 ! PSI mzData format
        .wiff   : MS:1000562 ! ABI WIFF format
        .pkl    : MS:1000565 ! Micromass PKL format
        .mzXML  : MS:1000566 ! ISB mzXML format
        .yep    : MS:1000567 ! Bruker/Agilent YEP format
        .dta    : MS:1000613 ! Sequest DTA format
        .mzMLb  : MS:1002838 ! mzMLb format
        
        Falls back to 'MS:1000560 ! mass spectrometer file format' if no match could be found.
        Upper/lowercase is ignored.
        
        Args:
            filepath: A filename (with optional path)
            
        Returns:
            tuple: (cv_accession, cv_name) e.g. ('MS:1000584', 'mzML format')
        """
        # Convert to lowercase for case-insensitive comparison
        filepath_lower = filepath.lower()
        
        # File extension to CV mapping
        format_mapping = {
            '.raw': ('MS:1000563', 'Thermo RAW format'),
            '.mzml': ('MS:1000584', 'mzML format'),
            '.mzdata': ('MS:1000564', 'PSI mzData format'),
            '.wiff': ('MS:1000562', 'ABI WIFF format'),
            '.pkl': ('MS:1000565', 'Micromass PKL format'),
            '.mzxml': ('MS:1000566', 'ISB mzXML format'),
            '.yep': ('MS:1000567', 'Bruker/Agilent YEP format'),
            '.dta': ('MS:1000613', 'Sequest DTA format'),
            '.mzmlb': ('MS:1002838', 'mzMLb format'),
        }
        
        # Check each known extension
        for extension, (accession, name) in format_mapping.items():
            if filepath_lower.endswith(extension):
                return accession, name
        
        # Fallback for unknown file formats
        return 'MS:1000560', 'mass spectrometer file format'