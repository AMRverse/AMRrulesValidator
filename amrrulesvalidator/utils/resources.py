"""Resource management for CARD and AMRFinderPlus data files."""

import io
import tarfile
import tempfile
import urllib.request
from pathlib import Path
from typing import Optional

import obonet


def extract_all_aro_terms(obo_file):
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
        # Create the resources directory if it doesn't exist
        self.dir.mkdir(parents=True, exist_ok=True)
        self._aro_terms_cache: Optional[list] = None
    
    def aro_terms(self) -> list:
        """Get ARO terms from cached OBO file, loading if necessary."""
        if self._aro_terms_cache is None:
            aro_obo_file = self.dir / "aro.obo"
            if aro_obo_file.exists():
                self._aro_terms_cache = extract_all_aro_terms(str(aro_obo_file))
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
    
    def download_card_archives(self):
        """
        Download and extract CARD ontology and data files into the resource directory.
        """
        # URLs for the CARD archives
        card_ontology_url = "https://card.mcmaster.ca/download/5/ontology-v4.0.1.tar.bz2"
        card_data_url = "https://card.mcmaster.ca/download/0/broadstreet-v4.0.1.tar.bz2"

        # Files to extract from each archive
        card_ontology_files = ["aro.obo", "ncbi_taxonomy.tsv"]
        card_data_files = ["aro_categories.tsv"]

        # Ensure the resource directory exists
        if not self.dir.exists():
            print(f"Creating resource directory: {self.dir}")
            self.dir.mkdir(parents=True, exist_ok=True)
            
        print(f"Resource directory: {self.dir}")

        # Helper function to download and extract files
        def download_and_extract(url, files_to_extract):
            with tempfile.NamedTemporaryFile(delete=False) as temp_file:
                # Download the archive
                print(f"Downloading {url}...")
                try:
                    with urllib.request.urlopen(url) as response:
                        temp_file.write(response.read())
                        temp_file.flush()
                except Exception as e:
                    print(f"Error downloading {url}: {e}")
                    return False

                temp_path = temp_file.name
            
            try:
                # Extract specific files
                with tarfile.open(temp_path, "r:bz2") as tar:
                    # List all members for debugging
                    all_members = tar.getmembers()
                    print(f"Archive contains {len(all_members)} files")
                    
                    # For each file we want to extract
                    for target_file in files_to_extract:
                        found = False
                        # Look for exact match or file within subdirectory
                        for member in all_members:
                            basename = Path(member.name).name
                            if basename == target_file:
                                print(f"Extracting {member.name} as {target_file}...")
                                # Extract but rename to the target filename
                                member_obj = tar.extractfile(member)
                                if member_obj:
                                    with open(self.dir / target_file, 'wb') as f:
                                        f.write(member_obj.read())
                                    found = True
                                    break
                        
                        if not found:
                            print(f"Warning: Could not find {target_file} in the archive")
                            
                return True
            except Exception as e:
                print(f"Error extracting files: {e}")
                return False
            finally:
                # Clean up the temp file
                try:
                    Path(temp_path).unlink(missing_ok=True)
                except:
                    pass

        # Download and extract files from both archives
        ontology_success = download_and_extract(card_ontology_url, card_ontology_files)
        data_success = download_and_extract(card_data_url, card_data_files)
        
        if ontology_success and data_success:
            print("CARD archives downloaded and extracted successfully.")
        else:
            print("Warning: Some files may not have been extracted correctly.")
            
        # Verify that the files exist
        missing_files = []
        for file_name in card_ontology_files + card_data_files:
            file_path = self.dir / file_name
            if not file_path.exists():
                missing_files.append(file_name)
                
        if missing_files:
            print(f"Warning: The following files are missing: {', '.join(missing_files)}")
            return False
        else:
            print("All required files are present in the resources directory.")
            return True
