"""
Module for analyzing phosphosites on-the-fly when structural data is not available in the database.
This module processes protein structure and sequence data to provide metrics for phosphosites.
"""

import io
import logging
from typing import Dict, List, Optional, Union
import numpy as np
from Bio.PDB import PDBParser, Selection, NeighborSearch

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def process_phosphosites_OTF(structure: str, sequence: str, phosphosites: List[Dict]) -> List[Dict]:
    """
    Process phosphorylation sites that don't have structural data in the database.
    Calculates pLDDT, nearby residues, surface accessibility, and other metrics on-the-fly.
    
    Args:
        structure: AlphaFold PDB structure string
        sequence: Protein sequence string
        phosphosites: List of phosphosite dictionaries
    
    Returns:
        Enhanced list of phosphosite dictionaries with calculated metrics
    """
    try:
        # Parse the structure
        parser = PDBParser(QUIET=True)
        structure_obj = parser.get_structure("protein", io.StringIO(structure))
        
        # Create neighbor search object for efficient proximity calculations
        all_atoms = list(structure_obj.get_atoms())
        ns = NeighborSearch(all_atoms)
        
        # Create a map of residue IDs to residue objects for quick lookup
        residue_map = {}
        for model in structure_obj:
            for chain in model:
                for residue in chain:
                    # Store only standard residues (not hetero-residues)
                    if residue.get_id()[0] == " ":
                        residue_map[residue.get_id()[1]] = residue
        
        # Process each phosphosite
        for site in phosphosites:
            # Skip sites that already have structural data
            if site.get('StructuralSimAvailable', False):
                continue
                
            if 'mean_plddt' in site and site['mean_plddt'] not in [None, 0, '0', 'N/A']:
                continue
                
            resno = site.get('resno')
            if not resno:
                continue
                
            # Process this site
            logger.info(f"Processing phosphosite {site.get('site')} on the fly")
            
            # Extract motif if not already present
            if 'motif' not in site or not site['motif']:
                site['motif'] = extract_motif_OTF(sequence, resno-1, 7)  # 0-indexed for sequence
                
            # Only proceed if the residue exists in the structure
            if resno in residue_map:
                target_residue = residue_map[resno]
                
                # Calculate site pLDDT
                site_plddt = calculate_site_plddt_OTF(target_residue)
                site['site_plddt'] = round(site_plddt, 2) if site_plddt is not None else 0
                
                # Calculate motif pLDDT
                motif_plddt = calculate_motif_plddt_OTF(structure_obj, sequence, resno, residue_map)
                site['mean_plddt'] = round(motif_plddt, 2) if motif_plddt is not None else 0
                
                # Get center atom for the residue
                center_atom = None
                if 'CA' in target_residue:
                    center_atom = target_residue['CA']
                else:
                    # Use the first atom if CA not available
                    try:
                        center_atom = next(target_residue.get_atoms())
                    except StopIteration:
                        # No atoms in residue
                        center_atom = None
                
                if center_atom:
                    # Calculate nearby residues
                    nearby_atoms = ns.search(center_atom.get_coord(), 10)  # 10Å radius
                    
                    # Count unique residues (excluding the target residue)
                    nearby_residues = set()
                    for atom in nearby_atoms:
                        parent = atom.get_parent()
                        if parent != target_residue:
                            nearby_residues.add(parent)
                    
                    site['nearby_count'] = len(nearby_residues)
                    
                    # Calculate surface accessibility
                    surface_accessibility = calculate_surface_accessibility_OTF(len(nearby_residues))
                    site['surface_accessibility'] = round(surface_accessibility, 2)
                else:
                    site['nearby_count'] = 0
                    site['surface_accessibility'] = 0.0
            else:
                # Residue not in structure
                site['site_plddt'] = 0
                site['mean_plddt'] = 0
                site['nearby_count'] = 0
                site['surface_accessibility'] = 0.0
                
            # Calculate AA frequencies if motif is available
            if site.get('motif'):
                calculate_aa_frequencies_OTF(site)
        
        return phosphosites
    except Exception as e:
        logger.error(f"Error processing phosphosites on the fly: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return phosphosites

def extract_motif_OTF(sequence: str, center_pos: int, window_size: int = 7) -> str:
    """
    Extract a motif from a sequence centered at a position with a specific window size.
    
    Args:
        sequence: The protein sequence
        center_pos: The center position (0-indexed)
        window_size: The number of residues to include on each side
    
    Returns:
        The motif sequence
    """
    if not sequence:
        return ""
        
    start = max(0, center_pos - window_size)
    end = min(len(sequence), center_pos + window_size + 1)
    
    # Pad with X if needed
    prefix = 'X' * max(0, window_size - center_pos)
    suffix = 'X' * max(0, window_size - (len(sequence) - center_pos - 1))
    
    return prefix + sequence[start:end] + suffix

def calculate_site_plddt_OTF(residue) -> Optional[float]:
    """
    Calculate the pLDDT score for a residue based on the B-factor.
    
    Args:
        residue: Biopython residue object
        
    Returns:
        pLDDT score or None if not available
    """
    try:
        # pLDDT is stored in the B-factor field in AlphaFold structures
        b_factors = []
        for atom in residue:
            b_factors.append(atom.get_bfactor())
            
        if b_factors:
            # Convert from bfactor to pLDDT (AlphaFold specific)
            return np.mean(b_factors)
        else:
            return None
    except Exception as e:
        logger.error(f"Error calculating site pLDDT: {e}")
        return None

def calculate_motif_plddt_OTF(structure, sequence, resno, residue_map) -> Optional[float]:
    """
    Calculate the average pLDDT for the motif around a phosphosite.
    
    Args:
        structure: Biopython structure object
        sequence: Protein sequence
        resno: Residue number of the phosphosite
        residue_map: Map of residue numbers to residue objects
        
    Returns:
        Average pLDDT for the motif or None if not available
    """
    try:
        # Define the window for the motif (-7 to +7)
        window_size = 7
        start_resno = max(1, resno - window_size)
        end_resno = min(len(sequence), resno + window_size)
        
        # Collect pLDDT values for each residue in the window
        plddt_values = []
        for i in range(start_resno, end_resno + 1):
            if i in residue_map:
                residue = residue_map[i]
                site_plddt = calculate_site_plddt_OTF(residue)
                if site_plddt is not None:
                    plddt_values.append(site_plddt)
        
        # Calculate average if we have values
        if plddt_values:
            return np.mean(plddt_values)
        else:
            return None
    except Exception as e:
        logger.error(f"Error calculating motif pLDDT: {e}")
        return None

def calculate_surface_accessibility_OTF(nearby_count: int) -> float:
    """
    Calculate a simple estimate of surface accessibility based on nearby residue count.
    Lower nearby counts typically indicate higher surface accessibility.
    
    Args:
        nearby_count: Number of residues within 10Å
        
    Returns:
        Estimated surface accessibility percentage
    """
    # Typically, buried residues have ~30-40 neighbors, surface residues ~10-15
    # Convert to a 0-100% scale where 0 = buried, 100 = fully exposed
    if nearby_count <= 0:
        return 100.0  # Likely an error, but treat as fully exposed
    
    # Values tuned based on empirical observations
    max_nearby = 40  # Completely buried
    min_nearby = 5   # Completely exposed
    
    # Constrain the nearby count
    constrained_count = max(min(nearby_count, max_nearby), min_nearby)
    
    # Linear mapping from nearby_count to accessibility
    accessibility = 100 * (max_nearby - constrained_count) / (max_nearby - min_nearby)
    return max(0, min(100, accessibility))  # Ensure result is between 0-100

def calculate_aa_frequencies_OTF(site: Dict) -> None:
    """
    Calculate amino acid type frequencies in the motif and add them to the site data.
    
    Args:
        site: Site dictionary with a motif field
        
    Returns:
        None (updates site dictionary in place)
    """
    motif = site.get('motif', '')
    if not motif or len(motif) < 3:
        return
    
    # Define amino acid groups
    aa_groups = {
        'polar': 'STYCNQ',
        'nonpolar': 'AVILMFWPG',
        'acidic': 'DE',
        'basic': 'KRH'
    }
    
    # Find the center position
    center_pos = len(motif) // 2
    
    # Count amino acids in each group (excluding the phosphosite itself)
    counts = {group: 0 for group in aa_groups}
    total = 0
    
    for i, aa in enumerate(motif):
        # Skip the phosphosite and any non-standard amino acids
        if i == center_pos or aa == 'X':
            continue
            
        total += 1
        for group, aas in aa_groups.items():
            if aa in aas:
                counts[group] += 1
                break
    
    # Calculate percentages
    if total > 0:
        for group, count in counts.items():
            site[f'{group}_aa_percent'] = round((count / total) * 100, 2)

import random
import re
from typing import List, Dict, Optional, Union

def enhance_phosphosites_with_default_flags(phosphosites: List[Dict]) -> List[Dict]:
    """
    Enhance a list of phosphosite dictionaries with default flag values.
    
    Args:
        phosphosites: List of phosphosite dictionaries
        
    Returns:
        Enhanced list of phosphosite dictionaries
    """
    if not phosphosites:
        return []
    
    # For each phosphosite, determine known and structural data flags
    for site in phosphosites:
        # Check if site is already properly flagged
        if 'is_known' not in site:
            # Set is_known flag based on is_known_phosphosite if available
            if 'is_known_phosphosite' in site:
                try:
                    is_known_val = int(site['is_known_phosphosite'])
                    site['is_known'] = bool(is_known_val)
                except (ValueError, TypeError):
                    # If we can't convert to int, use the bool() conversion
                    site['is_known'] = bool(site['is_known_phosphosite'])
            else:
                # Default to False if no information is available
                site['is_known'] = False
        
        # Set StructuralSimAvailable flag if not present
        if 'StructuralSimAvailable' not in site:
            # Set based on structural_sim_available if available
            if 'structural_sim_available' in site:
                try:
                    struct_sim_val = int(site['structural_sim_available'])
                    site['StructuralSimAvailable'] = bool(struct_sim_val)
                except (ValueError, TypeError):
                    # If we can't convert to int, use the bool() conversion
                    site['StructuralSimAvailable'] = bool(site['structural_sim_available'])
            else:
                # Default to False if no information is available
                site['StructuralSimAvailable'] = False
    
    return phosphosites

def mark_known_sites_by_pattern(phosphosites: List[Dict], site_patterns: List[str]) -> List[Dict]:
    """
    Mark phosphosites as known based on pattern matching.
    
    Args:
        phosphosites: List of phosphosite dictionaries
        site_patterns: List of site patterns to match (e.g., "S15", "T.*", "Y[1-9]")
        
    Returns:
        Updated list of phosphosite dictionaries
    """
    if not phosphosites or not site_patterns:
        return phosphosites
    
    # Compile patterns
    compiled_patterns = [re.compile(pattern) for pattern in site_patterns]
    
    # Mark sites as known if they match any pattern
    for site in phosphosites:
        site_str = site.get('site', '')
        if not site_str:
            continue
            
        for pattern in compiled_patterns:
            if pattern.match(site_str):
                site['is_known'] = True
                # Also set is_known_phosphosite for consistency
                site['is_known_phosphosite'] = 1
                break
    
    return phosphosites

def mark_structural_data_availability(phosphosites: List[Dict], available_percentage: float = 30) -> List[Dict]:
    """
    Mark a percentage of phosphosites as having structural data available.
    
    Args:
        phosphosites: List of phosphosite dictionaries
        available_percentage: Percentage of sites to mark as having data (0-100)
        
    Returns:
        Updated list of phosphosite dictionaries
    """
    if not phosphosites:
        return phosphosites
    
    # Determine number of sites to mark
    num_sites = len(phosphosites)
    num_to_mark = int(num_sites * available_percentage / 100)
    
    # Get random sample of indices
    random.seed(42)  # For reproducibility
    indices_to_mark = random.sample(range(num_sites), min(num_to_mark, num_sites))
    
    # Mark selected sites
    for i in indices_to_mark:
        phosphosites[i]['StructuralSimAvailable'] = True
        # Also set the original field for consistency
        phosphosites[i]['structural_sim_available'] = 1
    
    return phosphosites

def generate_demo_phosphosite_data(uniprot_id: str, sequence: str) -> List[Dict]:
    """
    Generate demo phosphosite data for a protein sequence.
    
    Args:
        uniprot_id: UniProt ID of the protein
        sequence: Protein sequence
        
    Returns:
        List of phosphosite dictionaries
    """
    if not sequence:
        return []
    
    phosphosites = []
    
    # Find all S, T, Y residues in the sequence
    for i, aa in enumerate(sequence):
        if aa in ['S', 'T', 'Y']:
            resno = i + 1  # 1-based indexing
            
            # Extract motif (-7 to +7)
            start = max(0, i - 7)
            end = min(len(sequence), i + 8)
            motif = sequence[start:end]
            
            # Generate random metrics
            mean_plddt = random.uniform(60, 90)
            site_plddt = random.uniform(60, 90)
            nearby_count = random.randint(5, 30)
            surface_accessibility = random.uniform(20, 80)
            
            # Create site dictionary
            site = {
                'site': f"{aa}{resno}",
                'resno': resno,
                'siteType': aa,
                'motif': motif,
                'mean_plddt': mean_plddt,
                'site_plddt': site_plddt,
                'nearby_count': nearby_count,
                'surface_accessibility': surface_accessibility,
                'is_known': False,
                'StructuralSimAvailable': False
            }
            
            phosphosites.append(site)
    
    # Mark some sites as known (20%)
    mark_known_sites = int(len(phosphosites) * 0.2)
    known_indices = random.sample(range(len(phosphosites)), mark_known_sites)
    for i in known_indices:
        phosphosites[i]['is_known'] = True
        phosphosites[i]['is_known_phosphosite'] = 1
    
    # Mark some sites as having structural data (30%)
    struct_indices = random.sample(range(len(phosphosites)), int(len(phosphosites) * 0.3))
    for i in struct_indices:
        phosphosites[i]['StructuralSimAvailable'] = True
        phosphosites[i]['structural_sim_available'] = 1
    
    return phosphosites

def find_all_sty_sites(sequence: str, existing_phosphosites: List[Dict] = None) -> List[Dict]:
    """
    Find all S, T, Y sites in a sequence and create basic phosphosite entries.
    Merges with existing phosphosites if provided.
    
    Args:
        sequence: Protein sequence
        existing_phosphosites: List of existing phosphosite dictionaries to merge with
        
    Returns:
        Complete list of phosphosite dictionaries with all S, T, Y sites
    """
    if not sequence:
        return existing_phosphosites or []
    
    # Create a lookup of existing sites by resno
    existing_lookup = {}
    if existing_phosphosites:
        for site in existing_phosphosites:
            if 'resno' in site:
                existing_lookup[site['resno']] = site
    
    # Create a complete list including all S, T, Y sites
    all_phosphosites = []
    
    # Find all S, T, Y residues in the sequence
    for i, aa in enumerate(sequence):
        if aa in ['S', 'T', 'Y']:
            resno = i + 1  # 1-based indexing
            
            # Check if we have existing data for this site
            if resno in existing_lookup:
                all_phosphosites.append(existing_lookup[resno])
                continue
            
            # Extract motif (-7 to +7)
            start = max(0, i - 7)
            end = min(len(sequence), i + 8)
            motif = sequence[start:end]
            
            # Create a basic site entry
            site = {
                'site': f"{aa}{resno}",
                'resno': resno,
                'siteType': aa,
                'motif': motif,
                'mean_plddt': 0,
                'site_plddt': 0,
                'nearby_count': 0,
                'surface_accessibility': 0,
                'is_known': False,
                'StructuralSimAvailable': False
            }
            
            all_phosphosites.append(site)
    
    return all_phosphosites