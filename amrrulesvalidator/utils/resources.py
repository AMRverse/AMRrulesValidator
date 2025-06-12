"""Resource management for CARD and AMRFinderPlus data files."""

import io
import urllib.request
from pathlib import Path
from typing import Optional

import obonet


def parse_obo_file(obo_file):
    # Load the OBO file
    aro_obo = obonet.read_obo(obo_file)

    # Extract all ARO terms
    aro_terms = [node for node in aro_obo.nodes if node.startswith('ARO:')]

    return aro_terms


def download_amrfp_files(url):

    # Use urllib to download the file
    with urllib.request.urlopen(url) as response:
        content = response.read().decode('utf-8')  # Decode the response as UTF-8
    
    # Parse the downloaded content as a TSV file
    content_io = io.StringIO(content)

    return content_io


class ResourceManager:
    """Manages cached resource files for validation."""
    
    def __init__(self):
        """Initialize resource manager with default resource directory."""
        self.dir = Path(__file__).parent.parent / "resources"
        self._aro_terms_cache: Optional[list] = None
    
    def aro_terms(self) -> list:
        """Get ARO terms from cached OBO file, loading if necessary."""
        if self._aro_terms_cache is None:
            aro_obo_file = self.dir / "aro.obo"
            if aro_obo_file.exists():
                self._aro_terms_cache = parse_obo_file(str(aro_obo_file))
            else:
                # Return empty list if file doesn't exist yet
                self._aro_terms_cache = []
        return self._aro_terms_cache
    
    def drug_names(self) -> list:
        """Get list of valid drug names from CARD data."""
        # TODO: Implement parsing of aro_categories.tsv or equivalent
        return []
    
    def drug_classes(self) -> list:
        """Get list of valid drug classes from CARD data."""
        # TODO: Implement parsing of aro_categories.tsv or equivalent
        return []
