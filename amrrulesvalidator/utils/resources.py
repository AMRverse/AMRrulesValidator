"""Resource management for CARD and AMRFinderPlus data files."""

import csv
import io
import tarfile
import tempfile
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

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
        self._drug_names_cache: Optional[list] = None
        self._drug_classes_cache: Optional[list] = None
        self._amrfp_db_version: Optional[str] = None
    
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
    
    def _load_drug_data(self):
        """Load drug data from card_drug_names.tsv if it exists."""
        drug_names_file = self.dir / "card_drug_names.tsv"
        if drug_names_file.exists():
            drug_names = []
            drug_classes = []
            with open(drug_names_file, 'r', newline='') as file:
                reader = csv.DictReader(file, delimiter='\t')
                for row in reader:
                    if 'Drug Name' in row and row['Drug Name'].strip():
                        drug_names.append(row['Drug Name'].strip())
                    if 'Drug Class' in row and row['Drug Class'].strip():
                        drug_classes.append(row['Drug Class'].strip())
            
            # Remove duplicates
            self._drug_names_cache = list(set(drug_names))
            self._drug_classes_cache = list(set(drug_classes))
        else:
            self._drug_names_cache = []
            self._drug_classes_cache = []
    
    def drug_names(self) -> list:
        """Get list of valid drug names from CARD data."""
        if self._drug_names_cache is None:
            self._load_drug_data()
        return self._drug_names_cache
    
    def drug_classes(self) -> list:
        """Get list of valid drug classes from CARD data."""
        if self._drug_classes_cache is None:
            self._load_drug_data()
        return self._drug_classes_cache
    
    def _parse_obo_for_descendants(self, obo_file_path: str, term_id: str) -> List[str]:
        """Find all descendants of a given term in the OBO file."""
        children = []
        with open(obo_file_path, 'r') as file:
            current_term = None
            for line in file:
                line = line.strip()
                if line.startswith("[Term]"):
                    current_term = None
                elif line.startswith("id: "):
                    current_term = line.split("id: ")[1]
                elif line.startswith("is_a: ") and current_term:
                    parent = line.split("is_a: ")[1].split(" ! ")[0]
                    if parent == term_id:
                        children.append(current_term)
        return children
    
    def _extract_card_drugs(self, obo_file_path: str, categories_file_path: str) -> List[Tuple[str, str, str]]:
        """Extract drug names and classes from CARD ontology files."""
        # Load the OBO file
        card_ontology = obonet.read_obo(obo_file_path)
        
        # Load the categories TSV
        with open(categories_file_path, 'r', newline='') as file:
            card_categories = csv.DictReader(file, delimiter='\t')
            
            # Get drug classes from categories
            drug_classes: Dict[str, str] = {}
            for row in card_categories:
                if row['ARO Category'] == 'Drug Class':
                    drug_classes[row['ARO Name']] = row['ARO Accession']
        
        # Get term names
        id_to_name = {id_: data.get("name") for id_, data in card_ontology.nodes(data=True)}
        
        # Collect results
        output_rows = []
        
        # Extract drugs for each drug class
        for drug_class, aro_accession in drug_classes.items():
            try:
                children = self._parse_obo_for_descendants(obo_file_path, aro_accession)
                for child in children:
                    if child in id_to_name:
                        child_name = id_to_name[child]
                        output_rows.append((child, child_name, drug_class))
            except Exception as e:
                print(f"Error processing {aro_accession} for {drug_class}: {e}")
                continue
        
        # Additional specific betalactam and other drug classes
        betalac_aros = [
            'ARO:3009105', 'ARO:3009106', 'ARO:3009107', 'ARO:3009108', 
            'ARO:3009109', 'ARO:3009123', 'ARO:3009124', 'ARO:3009125',
            'ARO:3000035', 'ARO:3007783', 'ARO:0000022', 'ARO:3007629',
            'ARO:3000707'
        ]
        
        for aro_accession in betalac_aros:
            if aro_accession in id_to_name:
                drug_class = id_to_name[aro_accession]
                children = self._parse_obo_for_descendants(obo_file_path, aro_accession)
                for child in children:
                    if child in id_to_name:
                        child_name = id_to_name[child]
                        output_rows.append((child, child_name, drug_class))
        
        return output_rows
    
    def _generate_card_drug_names_file(self):
        """Generate the card_drug_names.tsv file from CARD data files."""
        obo_file = self.dir / "aro.obo"
        categories_file = self.dir / "aro_categories.tsv"
        output_file = self.dir / "card_drug_names.tsv"
        
        if not obo_file.exists() or not categories_file.exists():
            print("Warning: Required CARD files are missing. Cannot generate drug names file.")
            return False
        
        try:
            # Extract drug data
            drug_data = self._extract_card_drugs(str(obo_file), str(categories_file))
            
            # Write to TSV file
            with open(output_file, 'w', newline='') as file:
                writer = csv.writer(file, delimiter='\t')
                writer.writerow(['ARO Accession', 'Drug Name', 'Drug Class'])
                for row in drug_data:
                    writer.writerow(row)
            
            print(f"Generated card_drug_names.tsv with {len(drug_data)} entries")
            
            # Reset cache to force reload
            self._drug_names_cache = None
            self._drug_classes_cache = None
            
            return True
        except Exception as e:
            print(f"Error generating drug names file: {e}")
            return False
    
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
            # Now generate the drug names file
            drugs_success = self._generate_card_drug_names_file()
            if drugs_success:
                print("Drug names file generated successfully.")
            else:
                print("Warning: Failed to generate drug names file.")
            return drugs_success and ontology_success and data_success
    
    def get_amrfp_db_version(self) -> str:
        """
        Get the AMRFinderPlus database version from the downloaded version.txt file.
        
        Returns:
            str: The version string or "Unknown" if the file is not available
        """
        version_file = self.dir / "version.txt"
        if version_file.exists():
            try:
                with open(version_file, 'r') as f:
                    self._amrfp_db_version = f.read().strip()
                return self._amrfp_db_version
            except Exception as e:
                print(f"Error reading AMRFinderPlus version file: {e}")
                return "Unknown"
        else:
            print("AMRFinderPlus version file not found.")
            return "Unknown"
    
    def download_amrfp_resources(self):
        """
        Download AMRFinderPlus reference files into the resource directory.
        """
        # URLs for the AMRFinderPlus resources
        refseq_url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneCatalog.txt'
        amrfp_nodes_url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/ReferenceGeneHierarchy.txt'
        amrfp_version_url = 'https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/version.txt'
        
        # Files to download
        file_urls = {
            'ReferenceGeneCatalog.txt': refseq_url,
            'ReferenceGeneHierarchy.txt': amrfp_nodes_url,
            'version.txt': amrfp_version_url
        }
        
        # TODO: HMM file is missing in this implementation
        # Need to determine where to get the HMM file from as it's not available via direct URL
        
        success = True
        for filename, url in file_urls.items():
            target_path = self.dir / filename
            print(f"Downloading {filename} from {url}...")
            
            try:
                with urllib.request.urlopen(url) as response:
                    content = response.read()
                    
                with open(target_path, 'wb') as f:
                    f.write(content)
                
                print(f"Successfully downloaded {filename}")
                
                # Get the database version
                if filename == 'version.txt':
                    self._amrfp_db_version = content.decode('utf-8').strip()
                    print(f"AMRFinderPlus database version: {self._amrfp_db_version}")
                    
            except Exception as e:
                print(f"Error downloading {filename}: {e}")
                success = False
        
        if success:
            print("AMRFinderPlus resources downloaded successfully.")
            # Make sure to update the cached version
            self.get_amrfp_db_version()
        else:
            print("Warning: Some AMRFinderPlus resources could not be downloaded.")
        
        return success
    
    def setup_all_resources(self):
        """
        Download and set up all required resources (CARD and AMRFinderPlus).
        """
        card_success = self.download_card_archives()
        amrfp_success = self.download_amrfp_resources()
        
        if card_success and amrfp_success:
            print("All resources have been successfully set up.")
            return True
        else:
            print("Warning: Some resources could not be set up properly.")
            return False
