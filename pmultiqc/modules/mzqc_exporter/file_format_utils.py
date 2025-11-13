"""
Utility functions for file format detection and CV term mapping.
"""


class FileFormatUtils:
    """
    Utility class for mapping file formats to PSI-MS CV terms.
    """
    
    # File extension to CV mapping
    _FORMAT_MAPPING = {
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
    
    # Fallback CV term for unknown formats
    _FALLBACK_CV = ('MS:1000560', 'mass spectrometer file format')
    
    @classmethod
    def filename_to_cv(cls, filepath: str) -> tuple[str, str]:
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
        if not filepath:
            return cls._FALLBACK_CV
            
        # Convert to lowercase for case-insensitive comparison
        filepath_lower = filepath.lower()
        
        # Check each known extension
        for extension, (accession, name) in cls._FORMAT_MAPPING.items():
            if filepath_lower.endswith(extension):
                return accession, name
        
        # Fallback for unknown file formats
        return cls._FALLBACK_CV