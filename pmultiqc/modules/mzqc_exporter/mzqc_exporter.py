import json
import os
from datetime import datetime
from pathlib import Path

from mzqc import MZQCFile as qc

from multiqc import config

class MzQCExporterModule():

    def __init__(self):
        from pmultiqc.modules.common.logging import get_logger
        self.log = get_logger(self.__class__.__module__)


    def export_mzqc(self, mzqc_data) -> bool | None:
        """
        Exports MZQC data to a MZQC file format.

        This function takes structured MZQC data and writes it to a file in the
        MZQC (Mass Spectrometry QC) format. The MZQC format is designed to
        capture quality control metrics from mass spectrometry experiments in a
        standardized way that can be shared and processed by various tools.

        Args:
            mzqc_data (dict): A dictionary containing the mzQC data structure.
                              This must include 'run_qualities' and 'set_qualities',
                              which are lists of QualityMetric objects

        Returns:
            bool | None: Returns True if the export was successful, False if it failed,
                         or None if no export was performed due to invalid input.

        Example:
            >>> exporter = MzQCExporterModule()
            >>> qm = qc.QualityMetric( ... )
            >>> mzqc_data = {
            ...     "run_qualities": [qm],
            ...     "set_qualities": [],
            ... }
            >>> exporter.export_mzqc(mzqc_data)

        Note:
            This function requires the 'mzqc' library to be installed and properly
            configured. The input data structure must conform to the MZQC schema
            specification.
        
        """
        self.log.info("Starting mzQC export ...")

        cv_ms = qc.ControlledVocabulary(
            name="Proteomics Standards Initiative Mass Spectrometry Ontology",
            version="4.1.212",
            uri="https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo")

        mzqc = qc.MzQcFile(version="1.0.0",
                           creationDate=datetime.now().isoformat(),
                           runQualities=mzqc_data["run_qualities"],
                           setQualities=mzqc_data["set_qualities"],
                           controlledVocabularies=[cv_ms])

        output_dir = Path(config.output_dir) if config.output_dir is not None else None
        if output_dir is not None:
            # TODO: set the pmultiqc-output by parameters
            mzqc_filename = os.path.join(output_dir, "pmultiqc.mzqc")
            with open(mzqc_filename, "w") as mzqc_file:
                mzqc_file.write(json.dumps(json.loads(qc.JsonSerialisable.to_json(mzqc)), indent=2))
    
        self.log.info(f"Done exporting mzQC to {mzqc_filename}")
        return True