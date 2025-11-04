import os

class MzQCExporterModule():

    def __init__(self):
        from pmultiqc.modules.common.logging import get_logger
        self.log = get_logger(self.__class__.__module__)


    def export_mzqc(self, mzqc_data) -> bool | None:
        self.log.info("Starting mzQC export ...")


        self.log.info(f"data: {mzqc_data}" )


        self.log.info("Done exporting mzQC")
        return True