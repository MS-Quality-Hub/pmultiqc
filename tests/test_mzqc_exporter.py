import json
import tempfile
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest
from mzqc import MZQCFile as qc

from pmultiqc.modules.mzqc_exporter.mzqc_exporter import MzQCExporterModule
from pmultiqc.modules.mzqc_exporter.maxquant_adapter import MaxQuantAdapter
from pmultiqc.modules.mzqc_exporter.mzidentml_adapter import MzIdentMLAdapter
from pmultiqc.modules.mzqc_exporter.file_format_utils import FileFormatUtils
from pmultiqc.modules.maxquant.maxquant import MQMetaData


class TestMzQCExporterModule:
    """
    Test class for MzQCExporterModule functionality.
    
    These tests verify that the mzQC exporter can create valid mzQC files
    with proper metadata and quality metrics.
    """

    @pytest.fixture
    def mzqc_exporter(self):
        """Create a fresh MzQCExporterModule instance for each test."""
        return MzQCExporterModule()

    @pytest.fixture
    def sample_peptide_data(self):
        """Sample peptide count data for testing."""
        return {
            "sample1": {"peptides": 1500, "proteins": 300},
            "sample2": {"peptides": 1200, "proteins": 250},
            "sample3": 1800  # Test simple integer value
        }

    @pytest.fixture
    def sample_protein_data(self):
        """Sample protein count data for testing."""
        return {
            "sample1": 300,
            "sample2": 250,
            "sample3": 350
        }

    def test_mzqc_exporter_initialization(self, mzqc_exporter):
        """Test that MzQCExporterModule initializes correctly."""
        assert mzqc_exporter is not None
        assert hasattr(mzqc_exporter, 'log')
        assert mzqc_exporter.run_quality_base_metadata == []
        assert mzqc_exporter.run_quality_metrics == {}
        assert mzqc_exporter.run_quality_metadata == {}

    def test_create_empty_run(self, mzqc_exporter):
        """Test that _create_empty_run initializes run data correctly."""
        label = "test_sample"
        mzqc_exporter._create_empty_run(label)
        
        assert label in mzqc_exporter.run_quality_metrics
        assert mzqc_exporter.run_quality_metrics[label] == []
        assert label in mzqc_exporter.run_quality_metadata
        assert mzqc_exporter.run_quality_metadata[label]['input_files'] == []
        assert mzqc_exporter.run_quality_metadata[label]['analysis_software'] == []

    def test_add_metric_to_run_quality(self, mzqc_exporter):
        """Test adding quality metrics to run quality."""
        label = "test_sample"
        
        # Create a test quality metric
        qm = qc.QualityMetric(
            accession="MS:1003250",
            name="count of identified peptidoforms",
            value=1500,
            unit={"accession": "UO:0000189", "name": "count unit"}
        )
        
        mzqc_exporter.add_metric_to_run_quality(label, qm)
        
        assert label in mzqc_exporter.run_quality_metrics
        assert len(mzqc_exporter.run_quality_metrics[label]) == 1
        assert mzqc_exporter.run_quality_metrics[label][0] == qm

    def test_add_base_metadata_to_run_quality(self, mzqc_exporter):
        """Test adding base metadata that applies to all runs."""
        # Create test analysis software
        analysis_software = qc.AnalysisSoftware(
            accession="MS:1001583",
            name="MaxQuant",
            version="2.0.3.0",
            uri="https://www.maxquant.org/"
        )
        
        mzqc_exporter.add_base_metadata_to_run_quality(analysis_software)
        
        assert len(mzqc_exporter.run_quality_base_metadata) == 1
        assert mzqc_exporter.run_quality_base_metadata[0] == analysis_software

    def test_add_metadata_for_run_quality_input_file(self, mzqc_exporter):
        """Test adding InputFile metadata to specific run quality."""
        label = "test_sample"
        
        input_file = qc.InputFile(
            name="test.raw",
            location="/path/to/test.raw",
            fileFormat=qc.CvParameter(accession="MS:1000563", name="Thermo RAW format")
        )
        
        mzqc_exporter.add_metadata_for_run_quality(label, input_file)
        
        assert label in mzqc_exporter.run_quality_metadata
        assert len(mzqc_exporter.run_quality_metadata[label]['input_files']) == 1
        assert mzqc_exporter.run_quality_metadata[label]['input_files'][0] == input_file

    def test_add_metadata_for_run_quality_analysis_software(self, mzqc_exporter):
        """Test adding AnalysisSoftware metadata to specific run quality."""
        label = "test_sample"
        
        analysis_software = qc.AnalysisSoftware(
            accession="MS:1001583",
            name="MaxQuant",
            version="2.0.3.0"
        )
        
        mzqc_exporter.add_metadata_for_run_quality(label, analysis_software)
        
        assert label in mzqc_exporter.run_quality_metadata
        assert len(mzqc_exporter.run_quality_metadata[label]['analysis_software']) == 1
        assert mzqc_exporter.run_quality_metadata[label]['analysis_software'][0] == analysis_software

    def test_add_peptide_id_count(self, mzqc_exporter, sample_peptide_data):
        """Test adding peptide identification count metrics."""
        mzqc_exporter.add_peptide_id_count(sample_peptide_data)
        
        # Check that metrics were added for all samples
        assert "sample1" in mzqc_exporter.run_quality_metrics
        assert "sample2" in mzqc_exporter.run_quality_metrics
        assert "sample3" in mzqc_exporter.run_quality_metrics
        
        # Check sample1 (mapping data)
        sample1_metrics = mzqc_exporter.run_quality_metrics["sample1"]
        assert len(sample1_metrics) == 1
        assert sample1_metrics[0].accession == "MS:1003250"
        assert sample1_metrics[0].name == "count of identified peptidoforms"
        assert sample1_metrics[0].value == 1800  # sum of peptides + proteins
        
        # Check sample3 (simple integer)
        sample3_metrics = mzqc_exporter.run_quality_metrics["sample3"]
        assert len(sample3_metrics) == 1
        assert sample3_metrics[0].value == 1800

    def test_add_protein_id_count(self, mzqc_exporter, sample_protein_data):
        """Test adding protein identification count metrics."""
        mzqc_exporter.add_protein_id_count(sample_protein_data)
        
        # Check that metrics were added for all samples
        assert "sample1" in mzqc_exporter.run_quality_metrics
        assert "sample2" in mzqc_exporter.run_quality_metrics
        assert "sample3" in mzqc_exporter.run_quality_metrics
        
        # Check sample1
        sample1_metrics = mzqc_exporter.run_quality_metrics["sample1"]
        assert len(sample1_metrics) == 1
        assert sample1_metrics[0].accession == "MS:1002404"
        assert sample1_metrics[0].name == "count of identified proteins"
        assert sample1_metrics[0].value == 300

    @patch('pmultiqc.modules.mzqc_exporter.mzqc_exporter.config')
    def test_create_export(self, mock_config, mzqc_exporter, sample_peptide_data):
        """Test creating and exporting mzQC file."""
        # Setup mock config
        with tempfile.TemporaryDirectory() as temp_dir:
            mock_config.output_dir = temp_dir
            
            # Add some test data
            mzqc_exporter.add_peptide_id_count(sample_peptide_data)
            
            # Add base metadata
            analysis_software = qc.AnalysisSoftware(
                accession="MS:1001583",
                name="MaxQuant",
                version="2.0.3.0"
            )
            mzqc_exporter.add_base_metadata_to_run_quality(analysis_software)
            
            # Create export
            mzqc_exporter.create_export()
            
            # Check that file was created
            mzqc_file_path = Path(temp_dir) / "pmultiqc.mzqc"
            assert mzqc_file_path.exists()
            
            # Verify file content
            with open(mzqc_file_path, 'r') as f:
                mzqc_data = json.load(f)
            
            assert "mzQC" in mzqc_data
            assert "version" in mzqc_data["mzQC"]
            assert "runQualities" in mzqc_data["mzQC"]
            assert "controlledVocabularies" in mzqc_data["mzQC"]
            # setQualities is optional and may not be present if no set qualities are added
            
            # Check that run qualities were created
            run_qualities = mzqc_data["mzQC"]["runQualities"]
            assert len(run_qualities) == 3  # sample1, sample2, sample3
            
            # Check that each run quality has the expected structure
            for rq in run_qualities:
                assert "metadata" in rq
                assert "qualityMetrics" in rq
                assert len(rq["qualityMetrics"]) == 1  # One peptide count metric


class TestMaxQuantAdapter:
    """
    Test class for MaxQuantAdapter functionality.
    
    These tests verify that the MaxQuant adapter correctly processes
    metadata and creates appropriate mzQC CV entries.
    """

    @pytest.fixture
    def mzqc_exporter(self):
        """Create a fresh MzQCExporterModule instance for each test."""
        return MzQCExporterModule()

    @pytest.fixture
    def maxquant_adapter(self, mzqc_exporter):
        """Create a MaxQuantAdapter instance for each test."""
        return MaxQuantAdapter(mzqc_exporter)

    @pytest.fixture
    def sample_metadata(self):
        """Sample MaxQuant metadata for testing."""
        return MQMetaData(
            version="2.0.3.0",
            fastafile_paths=["/path/to/database.fasta", "/path/to/contaminants.fasta"],
            rawfile_paths=["/path/to/sample1.raw", "/path/to/sample2.mzML", "/path/to/sample3.wiff"]
        )

    def test_maxquant_adapter_initialization(self, maxquant_adapter, mzqc_exporter):
        """Test that MaxQuantAdapter initializes correctly."""
        assert maxquant_adapter is not None
        assert maxquant_adapter.mzqc_exporter == mzqc_exporter

    def test_process_metadata(self, maxquant_adapter, sample_metadata):
        """Test processing of MaxQuant metadata."""
        maxquant_adapter.process_metadata(sample_metadata)
        
        # Check that base metadata was added (analysis software + FASTA files)
        assert len(maxquant_adapter.mzqc_exporter.run_quality_base_metadata) == 3  # 1 software + 2 FASTA files
        
        # Check that run-specific metadata was added for each raw file
        assert "sample1" in maxquant_adapter.mzqc_exporter.run_quality_metadata
        assert "sample2" in maxquant_adapter.mzqc_exporter.run_quality_metadata
        assert "sample3" in maxquant_adapter.mzqc_exporter.run_quality_metadata

    def test_create_analysis_software_entry(self, maxquant_adapter):
        """Test creation of analysis software CV entry."""
        version = "2.0.3.0"
        maxquant_adapter._create_analysis_software_entry(version)
        
        assert len(maxquant_adapter.mzqc_exporter.run_quality_base_metadata) == 1
        
        software = maxquant_adapter.mzqc_exporter.run_quality_base_metadata[0]
        assert isinstance(software, qc.AnalysisSoftware)
        assert software.accession == "MS:1001583"
        assert software.name == "MaxQuant"
        assert software.version == version
        assert software.uri == "https://www.maxquant.org/"

    def test_create_fasta_file_entries(self, maxquant_adapter):
        """Test creation of FASTA file CV entries."""
        fasta_paths = ["/path/to/database.fasta", "/path/to/contaminants.fasta"]
        maxquant_adapter._create_fasta_file_entries(fasta_paths)
        
        assert len(maxquant_adapter.mzqc_exporter.run_quality_base_metadata) == 2
        
        for i, fasta_file in enumerate(maxquant_adapter.mzqc_exporter.run_quality_base_metadata):
            assert isinstance(fasta_file, qc.InputFile)
            assert fasta_file.name == Path(fasta_paths[i]).name
            assert fasta_file.location == fasta_paths[i]
            assert fasta_file.fileFormat.accession == "MS:1001348"
            assert fasta_file.fileFormat.name == "FASTA format"

    def test_create_raw_file_entries(self, maxquant_adapter):
        """Test creation of raw file CV entries."""
        raw_paths = ["/path/to/sample1.raw", "/path/to/sample2.mzML"]
        maxquant_adapter._create_raw_file_entries(raw_paths)
        
        # Check that run-specific metadata was created
        assert "sample1" in maxquant_adapter.mzqc_exporter.run_quality_metadata
        assert "sample2" in maxquant_adapter.mzqc_exporter.run_quality_metadata
        
        # Check sample1 (RAW format)
        sample1_files = maxquant_adapter.mzqc_exporter.run_quality_metadata["sample1"]["input_files"]
        assert len(sample1_files) == 1
        assert sample1_files[0].name == "sample1.raw"
        assert sample1_files[0].location == "/path/to/sample1.raw"
        assert sample1_files[0].fileFormat.accession == "MS:1000563"
        assert sample1_files[0].fileFormat.name == "Thermo RAW format"
        
        # Check sample2 (mzML format)
        sample2_files = maxquant_adapter.mzqc_exporter.run_quality_metadata["sample2"]["input_files"]
        assert len(sample2_files) == 1
        assert sample2_files[0].name == "sample2.mzML"
        assert sample2_files[0].fileFormat.accession == "MS:1000584"
        assert sample2_files[0].fileFormat.name == "mzML format"



class TestMzQCIntegration:
    """
    Integration tests for the complete mzQC export workflow.
    """

    @pytest.fixture
    def complete_setup(self):
        """Set up a complete mzQC export scenario."""
        mzqc_exporter = MzQCExporterModule()
        adapter = MaxQuantAdapter(mzqc_exporter)
        
        # Sample metadata
        metadata = MQMetaData(
            version="2.0.3.0",
            fastafile_paths=["/path/to/database.fasta"],
            rawfile_paths=["/path/to/sample1.raw", "/path/to/sample2.mzML"]
        )
        
        # Sample data
        peptide_data = {"sample1": 1500, "sample2": 1200}
        protein_data = {"sample1": 300, "sample2": 250}
        
        return mzqc_exporter, adapter, metadata, peptide_data, protein_data

    @patch('pmultiqc.modules.mzqc_exporter.mzqc_exporter.config')
    def test_complete_mzqc_workflow(self, mock_config, complete_setup):
        """Test the complete mzQC export workflow."""
        mzqc_exporter, adapter, metadata, peptide_data, protein_data = complete_setup
        
        with tempfile.TemporaryDirectory() as temp_dir:
            mock_config.output_dir = temp_dir
            
            # Process metadata
            adapter.process_metadata(metadata)
            
            # Add quality metrics
            mzqc_exporter.add_peptide_id_count(peptide_data)
            mzqc_exporter.add_protein_id_count(protein_data)
            
            # Create export
            mzqc_exporter.create_export()
            
            # Verify the exported file
            mzqc_file_path = Path(temp_dir) / "pmultiqc.mzqc"
            assert mzqc_file_path.exists()
            
            with open(mzqc_file_path, 'r') as f:
                mzqc_data = json.load(f)
            
            # Verify structure
            assert "mzQC" in mzqc_data
            mzqc_content = mzqc_data["mzQC"]
            
            # Check run qualities
            run_qualities = mzqc_content["runQualities"]
            assert len(run_qualities) == 2  # sample1, sample2
            
            # Check that each run quality has the expected metrics and metadata
            for rq in run_qualities:
                # Should have 2 quality metrics (peptide + protein counts)
                assert len(rq["qualityMetrics"]) == 2
                
                # Should have metadata with input files and analysis software
                metadata = rq["metadata"]
                assert "inputFiles" in metadata
                assert "analysisSoftware" in metadata
                
                # Should have 2 input files (1 raw file + 1 FASTA file)
                assert len(metadata["inputFiles"]) == 2
                
                # Should have 1 analysis software (MaxQuant)
                assert len(metadata["analysisSoftware"]) == 1
                assert metadata["analysisSoftware"][0]["name"] == "MaxQuant"
            
            # Check controlled vocabularies
            assert len(mzqc_content["controlledVocabularies"]) == 1
            cv = mzqc_content["controlledVocabularies"][0]
            assert cv["name"] == "Proteomics Standards Initiative Mass Spectrometry Ontology"


class TestMzIdentMLAdapter:
    """
    Test class for MzIdentMLAdapter functionality.
    
    These tests verify that the mzIdentML adapter correctly processes
    metadata and creates appropriate mzQC CV entries.
    """

    @pytest.fixture
    def mzqc_exporter(self):
        """Create a fresh MzQCExporterModule instance for each test."""
        return MzQCExporterModule()

    @pytest.fixture
    def mzidentml_adapter(self, mzqc_exporter):
        """Create a MzIdentMLAdapter instance for each test."""
        return MzIdentMLAdapter(mzqc_exporter)

    @pytest.fixture
    def sample_mzml_df(self):
        """Sample mzML DataFrame for testing."""
        return pd.DataFrame({
            'filename': ['sample1', 'sample1', 'sample2', 'sample2', 'sample3'],
            'scan_number': [1, 2, 1, 2, 1],
            'retention_time': [10.5, 11.2, 9.8, 10.1, 12.3]
        })

    @pytest.fixture
    def sample_ms_paths(self):
        """Sample MS file paths for testing."""
        return [
            "/path/to/sample1.mzML",
            "/path/to/sample2.raw", 
            "/path/to/sample3.wiff"
        ]

    def test_mzidentml_adapter_initialization(self, mzidentml_adapter, mzqc_exporter):
        """Test that MzIdentMLAdapter initializes correctly."""
        assert mzidentml_adapter is not None
        assert mzidentml_adapter.mzqc_exporter == mzqc_exporter

    def test_process_metadata(self, mzidentml_adapter, sample_mzml_df, sample_ms_paths):
        """Test processing of mzIdentML metadata."""
        mzidentml_adapter.process_metadata(sample_mzml_df, sample_ms_paths)
        
        # Check that run-specific metadata was added for each unique sample
        unique_samples = set(sample_mzml_df["filename"].unique())
        assert len(unique_samples) == 3  # sample1, sample2, sample3
        
        for sample in unique_samples:
            assert sample in mzidentml_adapter.mzqc_exporter.run_quality_metadata
            input_files = mzidentml_adapter.mzqc_exporter.run_quality_metadata[sample]["input_files"]
            assert len(input_files) == 1
            assert input_files[0].name == sample

    def test_process_metadata_empty_dataframe(self, mzidentml_adapter):
        """Test processing with empty DataFrame."""
        empty_df = pd.DataFrame()
        ms_paths = ["/path/to/sample1.mzML"]
        
        # Should not raise an error
        mzidentml_adapter.process_metadata(empty_df, ms_paths)
        
        # Should not create any metadata
        assert len(mzidentml_adapter.mzqc_exporter.run_quality_metadata) == 0

    def test_process_metadata_none_dataframe(self, mzidentml_adapter):
        """Test processing with None DataFrame."""
        ms_paths = ["/path/to/sample1.mzML"]
        
        # Should not raise an error
        mzidentml_adapter.process_metadata(None, ms_paths)
        
        # Should not create any metadata
        assert len(mzidentml_adapter.mzqc_exporter.run_quality_metadata) == 0

    def test_create_input_file_entries(self, mzidentml_adapter, sample_mzml_df, sample_ms_paths):
        """Test creation of input file CV entries."""
        mzidentml_adapter._create_input_file_entries(sample_mzml_df, sample_ms_paths)
        
        # Check sample1 (mzML format)
        sample1_files = mzidentml_adapter.mzqc_exporter.run_quality_metadata["sample1"]["input_files"]
        assert len(sample1_files) == 1
        assert sample1_files[0].name == "sample1"
        assert sample1_files[0].location == "/path/to/sample1.mzML"
        assert sample1_files[0].fileFormat.accession == "MS:1000584"
        assert sample1_files[0].fileFormat.name == "mzML format"
        
        # Check sample2 (RAW format)
        sample2_files = mzidentml_adapter.mzqc_exporter.run_quality_metadata["sample2"]["input_files"]
        assert len(sample2_files) == 1
        assert sample2_files[0].name == "sample2"
        assert sample2_files[0].location == "/path/to/sample2.raw"
        assert sample2_files[0].fileFormat.accession == "MS:1000563"
        assert sample2_files[0].fileFormat.name == "Thermo RAW format"
        
        # Check sample3 (WIFF format)
        sample3_files = mzidentml_adapter.mzqc_exporter.run_quality_metadata["sample3"]["input_files"]
        assert len(sample3_files) == 1
        assert sample3_files[0].name == "sample3"
        assert sample3_files[0].location == "/path/to/sample3.wiff"
        assert sample3_files[0].fileFormat.accession == "MS:1000562"
        assert sample3_files[0].fileFormat.name == "ABI WIFF format"

    def test_find_sample_path(self, mzidentml_adapter):
        """Test finding sample paths."""
        ms_paths = [
            "/path/to/sample1.mzML",
            "/path/to/sample2_extra.raw",
            "/path/to/different_sample3.wiff"
        ]
        
        # Test exact match
        path = mzidentml_adapter._find_sample_path("sample1", ms_paths)
        assert path == "/path/to/sample1.mzML"
        
        # Test prefix match
        path = mzidentml_adapter._find_sample_path("sample2", ms_paths)
        assert path == "/path/to/sample2_extra.raw"
        
        # Test no match
        path = mzidentml_adapter._find_sample_path("nonexistent", ms_paths)
        assert path is None



class TestFileFormatUtils:
    """
    Test class for FileFormatUtils functionality.
    
    These tests verify that the file format utility correctly maps
    file extensions to PSI-MS CV terms.
    """

    def test_filename_to_cv_known_formats(self):
        """Test filename to CV mapping for known file formats."""
        test_cases = [
            ("test.raw", "MS:1000563", "Thermo RAW format"),
            ("test.mzML", "MS:1000584", "mzML format"),
            ("test.mzData", "MS:1000564", "PSI mzData format"),
            ("test.wiff", "MS:1000562", "ABI WIFF format"),
            ("test.pkl", "MS:1000565", "Micromass PKL format"),
            ("test.mzXML", "MS:1000566", "ISB mzXML format"),
            ("test.yep", "MS:1000567", "Bruker/Agilent YEP format"),
            ("test.dta", "MS:1000613", "Sequest DTA format"),
            ("test.mzMLb", "MS:1002838", "mzMLb format"),
        ]
        
        for filename, expected_accession, expected_name in test_cases:
            accession, name = FileFormatUtils.filename_to_cv(filename)
            assert accession == expected_accession, f"Failed for {filename}"
            assert name == expected_name, f"Failed for {filename}"

    def test_filename_to_cv_case_insensitive(self):
        """Test that filename to CV mapping is case insensitive."""
        test_cases = [
            ("TEST.RAW", "MS:1000563", "Thermo RAW format"),
            ("Test.MzML", "MS:1000584", "mzML format"),
            ("test.WIFF", "MS:1000562", "ABI WIFF format"),
        ]
        
        for filename, expected_accession, expected_name in test_cases:
            accession, name = FileFormatUtils.filename_to_cv(filename)
            assert accession == expected_accession, f"Failed for {filename}"
            assert name == expected_name, f"Failed for {filename}"

    def test_filename_to_cv_unknown_format(self):
        """Test filename to CV mapping for unknown file formats."""
        unknown_files = ["test.xyz", "test.unknown", "test"]
        
        for filename in unknown_files:
            accession, name = FileFormatUtils.filename_to_cv(filename)
            assert accession == "MS:1000560", f"Failed for {filename}"
            assert name == "mass spectrometer file format", f"Failed for {filename}"

    def test_filename_to_cv_empty_or_none(self):
        """Test filename to CV mapping with empty or None input."""
        test_cases = ["", None]
        
        for filename in test_cases:
            accession, name = FileFormatUtils.filename_to_cv(filename)
            assert accession == "MS:1000560"
            assert name == "mass spectrometer file format"

    def test_filename_to_cv_with_path(self):
        """Test filename to CV mapping with full file paths."""
        test_cases = [
            ("/path/to/data/sample.raw", "MS:1000563", "Thermo RAW format"),
            ("/path/to/data/sample.mzML", "MS:1000584", "mzML format"),
            ("C:\\data\\sample.wiff", "MS:1000562", "ABI WIFF format"),
        ]
        
        for filepath, expected_accession, expected_name in test_cases:
            accession, name = FileFormatUtils.filename_to_cv(filepath)
            assert accession == expected_accession, f"Failed for {filepath}"
            assert name == expected_name, f"Failed for {filepath}"