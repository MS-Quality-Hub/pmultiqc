import json

from datetime import datetime
from mzqc import MZQCFile as qc

class MzQCExporterModule():

    def __init__(self):
        from pmultiqc.modules.common.logging import get_logger
        self.log = get_logger(self.__class__.__module__)


    def export_mzqc(self, mzqc_data) -> bool | None:
        self.log.info("Starting mzQC export ...")

        cv_ms = qc.ControlledVocabulary(
            name="Proteomics Standards Initiative Mass Spectrometry Ontology",
            version="4.1.197",
            uri="https://github.com/HUPO-PSI/psi-ms-CV/blob/master/psi-ms.obo")


        self.log.info(f"data: {mzqc_data}" )

        mzqc = qc.MzQcFile(version="1.0.0",
                           creationDate=datetime.now().isoformat(),
                           runQualities=mzqc_data["run_qualities"],
                           setQualities=mzqc_data["set_qualities"],
                           controlledVocabularies=[cv_ms])

        print(f"mzQC: \n{json.dumps(json.loads(qc.JsonSerialisable.to_json(mzqc)), indent=2)}")
        # with open(mzqc_filename, "w") as mzqc_file:
        #     mzqc_file.write(json.dumps(json.loads(qc.JsonSerialisable.to_json(mzqc)), indent=2))

        self.log.info("Done exporting mzQC")
        return True