import os

from ..common.file_utils import file_prefix


# MaxQuant File Paths
def maxquant_file_path(find_log_files):

    required_files = [
        "parameters.txt",
        "summary.txt",
        "proteinGroups.txt",
        "evidence.txt",
        "msms.txt",
        "msScans.txt",
        "msmsScans.txt",
    ]

    maxquant_paths = {}

    txt_root = None

    # MaxQuant Data
    for maxquant_file in find_log_files("pmultiqc/maxquant_result", filecontents=False):
        if maxquant_file["fn"] in required_files:
            f_path = os.path.join(maxquant_file["root"], maxquant_file["fn"])
            maxquant_paths[file_prefix(f_path)] = f_path
            ## grab the path to the txt's, so we can search for the mqpar.xml two folders up
            if txt_root is None:
                txt_root = maxquant_file["root"]

    ## find ../../mqpar.xml (using manual glob, since multiqc will not look outside the search path)
    import glob
    maxquant_paths["mqpar"] = (glob.glob(os.path.join(txt_root, "../../mqpar.xml")) or [None])[0]

    # SDRF
    # "*sdrf.tsv"
    for sdrf_file in find_log_files("pmultiqc/sdrf", filecontents=False):
        maxquant_paths["sdrf"] = os.path.join(sdrf_file["root"], sdrf_file["fn"])

    return maxquant_paths