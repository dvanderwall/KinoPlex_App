"""
Flask web application for the Protein Explorer.
"""

from flask import Flask, render_template, request, jsonify, redirect, url_for
import os
import sys
import logging
import networkx as nx
import requests
import numpy as np
from Bio.PDB import PDBParser, Selection, NeighborSearch
import io
import re
import json
from typing import Dict, List, Optional
import shutil

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Add the parent directory to the path to import protein_explorer
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import protein_explorer as pe

# Import essential modules from protein_explorer
from protein_explorer.data.scaffold import (
    get_uniprot_id_from_gene,
    get_protein_by_id,
    get_alphafold_structure
)

# Database configuration - add before other imports to ensure proper initialization
USE_DATABASE = os.environ.get('USE_DATABASE', 'True').lower() in ('true', '1', 't', 'yes')

# Import database integration module
from protein_explorer.db.db_integration import (
    use_database,
    check_database_health,
    find_structural_matches,
    get_phosphosite_data,
    enhance_phosphosite,
    enhance_structural_matches,
    find_sequence_matches,
    get_site_kinase_scores,
    predict_kinases,
    get_similar_sites,
    get_network_kinase_scores,
    get_heatmap_data,
    get_all_phosphosites,
    get_kinase_comparison_data,
    get_network_kinase_comparison,
    get_known_kinase_info,
    categorize_kinases_by_family,
    analyze_motif_conservation
)

# Enable database mode if configured
if USE_DATABASE:
    use_database(True)
    logger.info("Database mode enabled")
else:
    logger.info("File-based mode enabled")

# Import other modules after database configuration
from protein_explorer.analysis.phospho import analyze_phosphosites
from protein_explorer.analysis.phospho_analyzer import (
    preload_structural_data, 
    get_phosphosites,
    analyze_phosphosite_context,
    enhance_site_visualization, 
    create_comparative_motif_visualization,
    analyze_residue_distributions
)

# Import visualization functions
from protein_explorer.visualization.network import create_phosphosite_network

from protein_explorer.analysis.sequence_analyzer import (
    preload_sequence_data,
    create_sequence_network_data,
    get_motif_enrichment,
    create_sequence_motif_visualization
)

from protein_explorer.analysis.kinase_predictor import (
    load_kinase_scores
)

# Import the network-based functions
from protein_explorer.analysis.network_kinase_predictor import (
    predict_kinases_network,
    get_network_heatmap_data, 
    get_kinase_family_distribution_network
)

# Perform a database health check on startup
try:
    if USE_DATABASE:
        logger.info("Performing database health check...")
        health_status = check_database_health()
        print("HEALTH STATUS:", health_status)
        if health_status['status'] == 'healthy':
            logger.info(f"Database connection healthy: {health_status.get('latency_ms', 'N/A')}ms latency")
        else:
            logger.warning(f"Database health check failed: {health_status.get('error', 'Unknown error')}")
            logger.warning("Falling back to file-based mode")
            use_database(False)
except Exception as e:
    logger.warning(f"Error checking database health: {e}")
    logger.warning("Falling back to file-based mode")
    use_database(False)

# Conditionally preload data files only if not using database
if not USE_DATABASE:
    try:
        print("Database mode is disabled. Preloading structural similarity data from files...")
        
        # Check for feather file first, fall back to parquet
        current_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(current_dir)
        feather_file = os.path.join(parent_dir, 'Combined_Kinome_10A_Master_Filtered_2.feather')
        if os.path.exists(feather_file):
            preload_structural_data(feather_file)
            print(f"Structural similarity data preloaded from: {feather_file}")
        else:
            parquet_file = os.path.join(parent_dir, 'Combined_Kinome_10A_Master_Filtered_2.parquet')
            if os.path.exists(parquet_file):
                preload_structural_data(parquet_file)
                print(f"Structural similarity data preloaded from: {parquet_file}")
            else:
                print("Warning: Could not find structural data files. Data will be loaded on first request (may cause delay).")
        
        print("Preloading sequence similarity data...")
        seq_file = os.path.join(parent_dir, 'Sequence_Similarity_Edges.parquet')
        if os.path.exists(seq_file):
            preload_sequence_data(seq_file)
            print(f"Sequence similarity data preloaded from: {seq_file}")
        else:
            print("Warning: Could not find sequence similarity data file.")
    except Exception as e:
        print(f"Warning: Failed to preload data files: {e}")
        print("Data will be loaded on first request (may cause delay)")
else:
    print("Database mode is enabled. Skipping file preloading.")

# Ensure cache directory exists
def ensure_cache_directory():
    """Ensure the cache directory exists at application startup."""
    # Paths for cache
    cache_dir = os.path.expanduser("~/.protein_explorer/cache")
    
    # Create cache directory if it doesn't exist
    if not os.path.exists(cache_dir):
        try:
            print(f"Creating cache directory: {cache_dir}")
            os.makedirs(cache_dir, exist_ok=True)
            print(f"Cache directory created successfully: {cache_dir}")
        except Exception as e:
            # If we can't create in home directory, use a temporary directory
            import tempfile
            alt_cache_dir = os.path.join(tempfile.gettempdir(), "protein_explorer_cache")
            print(f"Failed to create cache in home directory: {e}")
            print(f"Using alternative cache directory: {alt_cache_dir}")
            
            try:
                os.makedirs(alt_cache_dir, exist_ok=True)
            except Exception as e3:
                print(f"Failed to create alternative cache directory: {e3}")
                print("Application may have issues with caching")

# Clear and recreate the cache directory
cache_dir = os.path.expanduser("~/.protein_explorer/cache")
if os.path.exists(cache_dir):
    shutil.rmtree(cache_dir)
    os.makedirs(cache_dir, exist_ok=True)
    print("Cache cleared.")
else:
    print("Cache directory does not exist.")
ensure_cache_directory()

# Initialize Flask app
app = Flask(__name__)

@app.route('/')
def index():
    """Render the home page."""
    return render_template('index.html')

@app.route('/search', methods=['GET', 'POST'])
def search():
    """Search for a protein by UniProt ID or gene symbol."""
    if request.method == 'POST':
        # Get search parameters
        identifier = request.form.get('identifier', '')
        id_type = request.form.get('id_type', 'uniprot')
        
        if not identifier:
            return render_template('search.html', error="Please enter an identifier")
        
        try:
            # Redirect to protein page
            return redirect(url_for('protein', identifier=identifier, id_type=id_type))
        except Exception as e:
            logger.error(f"Error in search: {e}")
            return render_template('search.html', error=str(e))
    
    # GET request
    return render_template('search.html')

@app.route('/protein/<identifier>')
def protein(identifier):
    """Display protein information and structure with full supplementary structural data."""
    try:
        print(f"DEBUG: Processing protein view for identifier: {identifier}")
        id_type = request.args.get('id_type', 'uniprot')
        
        # Get protein data
        if id_type.lower() == 'uniprot':
            print(f"DEBUG: Getting protein by UniProt ID")
            protein_data = get_protein_by_id(uniprot_id=identifier)
        else:
            print(f"DEBUG: Getting protein by gene symbol")
            pre_filt_id = get_uniprot_id_from_gene(identifier)
            protein_data = get_protein_by_id(pre_filt_id)
        
        print(f"DEBUG: Protein data retrieved, has_structure: {protein_data.get('has_structure', False)}")
        
        # TEMPORARY FIX: For specific proteins, manually set has_structure=True
        if identifier in ['P04637', 'P53_HUMAN']:
            protein_data['has_structure'] = True
            print(f"DEBUG: Manually set has_structure=True for {identifier}")
        
        # Check if protein has a structure
        structure_html = None
        structure = None
        phosphosite_html = None
        phosphosites = None
        
        if protein_data.get('has_structure', False):
            print(f"DEBUG: Protein has structure, retrieving...")
            
            # Get structure 
            try:
                # Try direct download first for known proteins
                if identifier in ['P04637', 'P53_HUMAN']:
                    import requests
                    url = f"https://alphafold.ebi.ac.uk/files/AF-{protein_data['uniprot_id']}-F1-model_v4.pdb"
                    response = requests.get(url)
                    if response.status_code == 200:
                        structure = response.text
                        print(f"DEBUG: Got structure directly, length: {len(structure)}")
                    else:
                        # Fall back to normal method
                        structure = get_alphafold_structure(protein_data['uniprot_id'])
                else:
                    # Use normal method
                    structure = get_alphafold_structure(protein_data['uniprot_id'])
                
                # Get the protein sequence from metadata if available
                sequence = None
                if protein_data.get('metadata', {}).get('sequence', {}).get('value'):
                    sequence = protein_data['metadata']['sequence']['value']
                
                # Create structure visualization
                if structure:
                    print(f"DEBUG: Creating structure visualization")
                    structure_html = pe.visualization.visualize_structure(
                        structure,
                        sequence=sequence
                    )
                    
                    # Process phosphorylation sites - NEW APPROACH
                    if sequence:
                        print(f"DEBUG: Processing all potential phosphorylation sites from sequence")
                        
                        # Helper function to extract motif
                        def extract_motif(sequence, center_pos, window_size):
                            """Extract a motif of size 2*window_size+1 around center_pos."""
                            start = max(0, center_pos - window_size)
                            end = min(len(sequence), center_pos + window_size + 1)
                            
                            # Pad with X if needed
                            prefix = 'X' * max(0, window_size - center_pos)
                            suffix = 'X' * max(0, window_size - (len(sequence) - center_pos - 1))
                            
                            return prefix + sequence[start:end] + suffix
                        
                        # Find all potential STY sites in the sequence
                        all_phosphosites = []
                        for i, aa in enumerate(sequence):
                            if aa in ['S', 'T', 'Y']:
                                resno = i + 1  # 1-based residue numbering
                                all_phosphosites.append({
                                    'site': f"{aa}{resno}",
                                    'resno': resno,
                                    'siteType': aa,
                                    'motif': extract_motif(sequence, i, 7),  # Extract motif around the site
                                    'mean_plddt': 0,
                                    'site_plddt': 0,
                                    'nearby_count': 0,
                                    'surface_accessibility': 0,
                                    'is_known': False,  # Default to not known
                                    'StructuralSimAvailable': False  # Default to no structural data
                                })
                        
                        print(f"DEBUG: Found {len(all_phosphosites)} potential STY sites in the sequence")
                        
                        # Get known phosphosites from the database
                        try:
                            known_phosphosites = get_all_phosphosites(protein_data['uniprot_id'])
                            print(f"DEBUG: Retrieved {len(known_phosphosites)} known phosphosites from database")
                            
                            # Create a lookup dictionary for known sites
                            known_site_lookup = {}
                            for site in known_phosphosites:
                                if 'site' in site:
                                    site_key = site['site']
                                    known_site_lookup[site_key] = site

                                    # Print the first few raw site entries for debugging
                                    if len(known_site_lookup) <= 5:
                                        print(f"DEBUG: Raw DB entry for {site_key}:")
                                        for key, value in site.items():
                                            if key != 'full_data':  # Skip full_data to avoid huge output
                                                print(f"    {key}: {value}")

                            # Merge the information
                            for site in all_phosphosites:
                                site_key = site['site']
                                if site_key in known_site_lookup:
                                    # Update with known site data
                                    known_site = known_site_lookup[site_key]
                                    
                                    # Copy relevant fields
                                    for field in ['motif', 'mean_plddt', 'site_plddt', 'nearby_count', 
                                                'surface_accessibility']:
                                        if field in known_site and known_site[field] is not None:
                                            # Round numeric values if needed
                                            if field in ['mean_plddt', 'site_plddt', 'surface_accessibility']:
                                                try:
                                                    site[field] = round(float(known_site[field]), 2)
                                                except (ValueError, TypeError):
                                                    site[field] = known_site[field]
                                            else:
                                                site[field] = known_site[field]
                                    
                                    # Explicitly set is_known based on is_known_phosphosite field
                                    # First, try is_known_phosphosite
                                    if 'is_known_phosphosite' in known_site:
                                        try:
                                            # First try to convert to int and then to boolean
                                            is_known_val = int(float(known_site['is_known_phosphosite']))
                                            site['is_known'] = bool(is_known_val)
                                            print(f"DEBUG: Setting is_known={site['is_known']} for {site_key}, raw value was {known_site['is_known_phosphosite']}")
                                        except (ValueError, TypeError):
                                            # If that fails, try direct boolean conversion
                                            try:
                                                site['is_known'] = bool(known_site['is_known_phosphosite'])
                                                print(f"DEBUG: Setting is_known={site['is_known']} for {site_key}, using bool conversion on {known_site['is_known_phosphosite']}")
                                            except (ValueError, TypeError):
                                                # If even that fails, fallback to True because it's in the database
                                                site['is_known'] = True
                                                print(f"DEBUG: Setting is_known=True for {site_key}, fallback because the site is in the database")
                                    elif 'is_known' in known_site:
                                        # Try the is_known field next
                                        try:
                                            site['is_known'] = bool(int(known_site['is_known']))
                                        except (ValueError, TypeError):
                                            try:
                                                site['is_known'] = bool(known_site['is_known'])
                                            except (ValueError, TypeError):
                                                site['is_known'] = True
                                    else:
                                        # If site is in known_site_lookup but no is_known_phosphosite field,
                                        # we'll assume it's known (since it's in the phosphosite table)
                                        site['is_known'] = True
                                        print(f"DEBUG: Setting is_known=True for {site_key} (no is_known_phosphosite field)")
                                            
                                    # Explicitly handle StructuralSimAvailable field
                                    if 'StructuralSimAvailable' in known_site:
                                        # Convert to integer and then to boolean
                                        try:
                                            struct_sim_val = int(float(known_site['StructuralSimAvailable']))
                                            site['StructuralSimAvailable'] = bool(struct_sim_val)
                                            print(f"DEBUG: Setting StructuralSimAvailable={site['StructuralSimAvailable']} for {site_key}, raw value was {known_site['StructuralSimAvailable']}")
                                        except (ValueError, TypeError):
                                            # If it can't be converted to int, try direct boolean conversion
                                            try:
                                                site['StructuralSimAvailable'] = bool(known_site['StructuralSimAvailable'])
                                                print(f"DEBUG: Setting StructuralSimAvailable={site['StructuralSimAvailable']} for {site_key}, using bool conversion on {known_site['StructuralSimAvailable']}")
                                            except (ValueError, TypeError):
                                                # If even that fails, default to False
                                                site['StructuralSimAvailable'] = False
                                                print(f"DEBUG: Setting StructuralSimAvailable=False for {site_key}, conversion failed")
                                    else:
                                        # Default to False if field is missing
                                        site['StructuralSimAvailable'] = False
                                        print(f"DEBUG: Setting StructuralSimAvailable=False for {site_key} (field missing)")

                            # Print the first few sites to verify the merge
                            for i, site in enumerate(all_phosphosites[:10]):
                                print(f"DEBUG: Site {i+1}: {site['site']}, is_known={site['is_known']}, StructuralSimAvailable={site.get('StructuralSimAvailable', False)}")
                                
                        except Exception as e:
                            print(f"DEBUG: Error getting known phosphosites: {e}")
                            import traceback
                            print(traceback.format_exc())
                            
                            # If database query fails, use basic analysis
                            try:
                                # As a fallback, analyze structure to get pLDDT and nearby counts
                                analyzed_sites = analyze_phosphosites(sequence, structure, uniprot_id=protein_data.get('uniprot_id'))
                                
                                # Create a lookup for the analyzed sites
                                analyzed_lookup = {site['site']: site for site in analyzed_sites if 'site' in site}
                                
                                # Merge with all_phosphosites
                                for site in all_phosphosites:
                                    site_key = site['site']
                                    if site_key in analyzed_lookup:
                                        analyzed_site = analyzed_lookup[site_key]
                                        # Copy relevant fields
                                        for field in ['mean_plddt', 'nearby_count', 'motif', 'is_known']:
                                            if field in analyzed_site and analyzed_site[field] is not None:
                                                # Round numeric values if needed
                                                if field in ['mean_plddt', 'site_plddt', 'surface_accessibility']:
                                                    try:
                                                        site[field] = round(float(analyzed_site[field]), 2)
                                                    except (ValueError, TypeError):
                                                        site[field] = analyzed_site[field]
                                                else:
                                                    site[field] = analyzed_site[field]
                            except Exception as fallback_e:
                                print(f"DEBUG: Error in fallback analysis: {fallback_e}")
                                
                        # Process sites without structural data on-the-fly
                        try:
                            # Import the on-the-fly analyzer
                            from protein_explorer.analysis.OnTheFly import process_phosphosites_OTF
                            
                            # Process sites with missing structural data
                            all_phosphosites = process_phosphosites_OTF(structure, sequence, all_phosphosites)
                            print(f"DEBUG: Processed phosphosites on-the-fly")
                        except ImportError:
                            print(f"DEBUG: OnTheFly module not found, skipping on-the-fly processing")
                        except Exception as e:
                            print(f"DEBUG: Error processing phosphosites on-the-fly: {e}")
                            import traceback
                            print(traceback.format_exc())
                        
                        # Create enhanced table visualization using the merged data
                        try:
                            from protein_explorer.analysis.enhanced_table import enhance_phosphosite_table
                            phosphosite_html = enhance_phosphosite_table(all_phosphosites, protein_data['uniprot_id'])
                            phosphosites = all_phosphosites  # Store for template access
                        except Exception as e2:
                            print(f"DEBUG: Error using enhanced table: {e2}")
                            # Fall back to original HTML generation if enhanced_table fails
                            
                            phosphosite_html = f"""
                            <div class="card mt-4">
                                <div class="card-header">
                                    <h5 class="mb-0">Phosphorylation Site Analysis</h5>
                                </div>
                                <div class="card-body p-0">
                                    <div class="table-responsive">
                                        <table class="table table-striped table-hover phosphosite-table">
                                            <thead class="thead-light">
                                                <tr>
                                                    <th>Site</th>
                                                    <th>Motif (-7 to +7)</th>
                                                    <th>Mean pLDDT</th>
                                                    <th>Site pLDDT</th>
                                                    <th>Nearby Residues (10Å)</th>
                                                    <th>Surface Access.</th>
                                                    <th>Known</th>
                                                    <th>Actions</th>
                                                </tr>
                                            </thead>
                                            <tbody id="phosphosite-table">
                            """
                            
                            for site in all_phosphosites:
                                # Determine row style based on structural data availability
                                row_style = ""
                                if site.get('StructuralSimAvailable', False):
                                    row_style = 'style="background-color: #c8e6c9;"'  # Green
                                else:
                                    row_style = 'style="background-color: #ffcc80;"'  # Orange
                                
                                # Format numeric values for display
                                try:
                                    mean_plddt = float(site.get('mean_plddt', 0))
                                    mean_plddt_text = f"{mean_plddt:.2f}" if mean_plddt > 0 else "N/A"
                                except (ValueError, TypeError):
                                    mean_plddt_text = site.get('mean_plddt', 'N/A')
                                
                                try:
                                    site_plddt = float(site.get('site_plddt', 0))
                                    site_plddt_text = f"{site_plddt:.2f}" if site_plddt > 0 else "N/A"
                                except (ValueError, TypeError):
                                    site_plddt_text = site.get('site_plddt', 'N/A')
                                
                                try:
                                    surface_access = float(site.get('surface_accessibility', 0))
                                    surface_access_text = f"{surface_access:.2f}%" if surface_access > 0 else "N/A"
                                except (ValueError, TypeError):
                                    surface_access_text = site.get('surface_accessibility', 'N/A')
                                
                                # Add data attributes to the row for better visualization
                                data_attrs = f"""
                                    data-site="{site['site']}" 
                                    data-resno="{site['resno']}" 
                                    data-type="{site['siteType']}"
                                    data-plddt="{site.get('mean_plddt', 0)}" 
                                    data-nearby="{site.get('nearby_count', 0)}"
                                    data-surface="{site.get('surface_accessibility', 0)}"
                                    data-known="{str(site.get('is_known', False)).lower()}"
                                """
                                
                                phosphosite_html += f"""
                                <tr {row_style} {data_attrs}>
                                    <td><a href="/site/{protein_data['uniprot_id']}/{site['site']}" class="site-link" data-resno="{site['resno']}">{site['site']}</a></td>
                                    <td><code class="motif-sequence">{site.get('motif', '')}</code></td>
                                    <td>{mean_plddt_text}</td>
                                    <td>{site_plddt_text}</td>
                                    <td>{site.get('nearby_count', 'N/A')}</td>
                                    <td>{surface_access_text}</td>
                                    <td>{"Yes" if site.get('is_known', False) else "No"}</td>
                                    <td>
                                        <a href="/site/{protein_data['uniprot_id']}/{site['site']}" class="btn btn-sm btn-outline-primary">
                                            Details
                                        </a>
                                    </td>
                                </tr>
                                """
                            
                            phosphosite_html += """
                                            </tbody>
                                        </table>
                                    </div>
                                </div>
                            </div>
                            
                            <script>
                                document.addEventListener('DOMContentLoaded', function() {
                                    const siteLinks = document.querySelectorAll('.site-link');
                                    siteLinks.forEach(link => {
                                        link.addEventListener('click', function(e) {
                                            e.preventDefault();
                                            const resno = this.getAttribute('data-resno');
                                            
                                            // Find the span in the sequence viewer
                                            const sequenceSpans = document.querySelectorAll('.sequence-viewer span');
                                            if (sequenceSpans.length > 0) {
                                                // Find and click the span for this residue
                                                const index = parseInt(resno) - 1;
                                                if (index >= 0 && index < sequenceSpans.length) {
                                                    sequenceSpans[index].click();
                                                }
                                            }
                                        });
                                    });
                                });
                            </script>
                            """
                    else:
                        print(f"DEBUG: No sequence available for phosphosite analysis")
                        phosphosite_html = "<div class='alert alert-warning'>No protein sequence available for phosphosite analysis.</div>"
            except Exception as e:
                print(f"DEBUG: Error getting structure: {e}")
                structure_html = f'<div class="alert alert-danger">Error loading structure: {str(e)}</div>'
        
        print(f"DEBUG: Rendering template with phosphosite_html: {'Present' if phosphosite_html else 'None'}")
        
        # Make sure we include necessary scripts for visualization
        extra_scripts = """
        <script src="https://cdn.jsdelivr.net/npm/chart.js@3.9.1/dist/chart.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-annotation@2.0.1/dist/chartjs-plugin-annotation.min.js"></script>
        <script>
        // Check if Chart.js loaded successfully
        document.addEventListener('DOMContentLoaded', function() {
            console.log('DOM loaded, Chart.js status:', typeof Chart !== 'undefined' ? 'Loaded' : 'Not loaded');
            
            // If not loaded, try to load it again
            if (typeof Chart === 'undefined') {
                console.warn('Chart.js not loaded, attempting to load it again');
                const script = document.createElement('script');
                script.src = 'https://cdn.jsdelivr.net/npm/chart.js@3.9.1/dist/chart.min.js';
                script.onload = function() {
                    console.log('Chart.js loaded successfully on retry');
                    // Load visualization script
                    const vizScript = document.createElement('script');
                    vizScript.src = "{{ url_for('static', filename='js/phosphosite-visualization.js') }}";
                    document.head.appendChild(vizScript);
                };
                document.head.appendChild(script);
            }
        });
        </script>
        """
        
        return render_template(
            'protein.html',
            protein=protein_data,
            structure_html=structure_html,
            phosphosites=phosphosites,  # Pass the phosphosites data
            phosphosite_html=phosphosite_html,  # Pass the HTML for rendering
            extra_scripts=extra_scripts  # Add extra scripts for visualization
        )
    except Exception as e:
        print(f"DEBUG: Exception in protein view: {e}")
        logger.error(f"Error in protein view: {e}")
        return render_template('error.html', error=str(e))


@app.route('/analyze', methods=['GET', 'POST'])
def analyze():
    """Analyze multiple proteins."""
    if request.method == 'POST':
        # Get proteins from form
        proteins_text = request.form.get('proteins', '')
        analysis_type = request.form.get('analysis_type', 'network')
        
        # Parse protein list
        protein_list = [p.strip() for p in proteins_text.split(',') if p.strip()]
        
        if not protein_list:
            return render_template('analyze.html', error="Please enter at least one protein")
        
        try:
            # Determine ID type (assume all are the same type)
            id_type = 'gene' if '_' not in protein_list[0] else 'uniprot'
            
            # Convert to UniProt IDs if needed
            if id_type == 'gene':
                uniprot_ids = []
                for gene in protein_list:
                    try:
                        protein = pe.data.get_protein_by_id(gene_symbol=gene)
                        uniprot_ids.append(protein['uniprot_id'])
                    except Exception as e:
                        logger.warning(f"Could not convert gene {gene} to UniProt ID: {e}")
                protein_list = uniprot_ids
            
            # Perform the analysis
            if analysis_type == 'network':
                # Build network
                network = pe.navigation.build_interaction_network(protein_list, max_depth=1)
                
                # Create network visualization
                network_html = pe.visualization.visualize_network(
                    network,
                    highlight_nodes=protein_list,
                    title=f"Protein Interaction Network"
                )
                
                # Find common interactors
                common_interactors = pe.navigation.find_common_interactors(network, protein_list)
                
                return render_template(
                    'analyze_results.html',
                    analysis_type=analysis_type,
                    proteins=protein_list,
                    network_html=network_html,
                    common_interactors=common_interactors,
                    node_count=network.number_of_nodes(),
                    edge_count=network.number_of_edges()
                )
            elif analysis_type == 'structure':
                # Get structures and compare
                structures = {}
                for uniprot_id in protein_list:
                    try:
                        structure = get_alphafold_structure(uniprot_id)
                        if structure:
                            structures[uniprot_id] = structure
                    except Exception as e:
                        logger.warning(f"Could not get structure for {uniprot_id}: {e}")
                
                if not structures:
                    return render_template('analyze.html', 
                                          error="Could not find structures for any of the proteins")
                
                # Create structure visualization
                structure_html = None
                if len(structures) == 1:
                    # Single structure
                    uniprot_id, structure = next(iter(structures.items()))
                    structure_html = pe.visualization.visualize_structure(structure)
                else:
                    # Compare structures
                    structure_html = pe.visualization.compare_structures(list(structures.values()))
                
                return render_template(
                    'analyze_results.html',
                    analysis_type=analysis_type,
                    proteins=protein_list,
                    structure_html=structure_html,
                    structures_found=list(structures.keys())
                )
            else:
                return render_template('analyze.html', error=f"Unknown analysis type: {analysis_type}")
                
        except Exception as e:
            logger.error(f"Error in analysis: {e}")
            return render_template('analyze.html', error=str(e))
    
    # GET request
    return render_template('analyze.html')

@app.route('/api/phosphosites/<uniprot_id>', methods=['GET'])
def api_phosphosites(uniprot_id):
    """API endpoint for phosphorylation site analysis."""
    try:
        # Get all phosphosites for the protein using integration function
        phosphosites = get_all_phosphosites(uniprot_id)
        
        # If no phosphosites found through DB, try original analysis method
        if not phosphosites:
            # Get protein data
            protein_data = get_protein_by_id(uniprot_id=uniprot_id)
            
            # Get sequence
            sequence = protein_data.get('metadata', {}).get('sequence', {}).get('value')
            if not sequence:
                return jsonify({'error': 'Protein sequence not found'}), 404
                
            # Get structure
            structure = get_alphafold_structure(uniprot_id)
            if not structure:
                return jsonify({'error': 'Protein structure not found'}), 404
                
            # Analyze phosphosites
            phosphosites = analyze_phosphosites(sequence, structure, uniprot_id)
            
            # Enhance with supplementary data
            for site in phosphosites:
                if 'resno' in site:
                    site_id = f"{uniprot_id}_{site['resno']}"
                    supp_data = get_phosphosite_data(site_id)
                    
                    if supp_data:
                        for key, value in supp_data.items():
                            if key not in site and value is not None:
                                site[key] = value
        
        return jsonify(phosphosites)
    except Exception as e:
        logger.error(f"Error in API endpoint: {e}")
        return jsonify({'error': str(e)}), 400

@app.route('/api/structure/<uniprot_id>', methods=['GET'])
def api_structure(uniprot_id):
    """API endpoint for protein structure."""
    try:
        structure = get_alphafold_structure(uniprot_id)
        if not structure:
            return jsonify({'error': 'Structure not found'}), 404
        
        # Return basic info about the structure
        structure_info = pe.data.parse_pdb_structure(structure)
        return jsonify(structure_info)
    except Exception as e:
        return jsonify({'error': str(e)}), 400

@app.route('/api/network/<uniprot_id>', methods=['GET'])
def api_network(uniprot_id):
    """API endpoint for protein network."""
    try:
        depth = int(request.args.get('depth', 1))
        network = pe.navigation.build_interaction_network([uniprot_id], max_depth=depth)
        
        # Convert network to dict for JSON serialization
        network_data = {
            'nodes': list(network.nodes()),
            'edges': list(network.edges()),
            'node_count': network.number_of_nodes(),
            'edge_count': network.number_of_edges()
        }
        
        return jsonify(network_data)
    except Exception as e:
        return jsonify({'error': str(e)}), 400

    
@app.route('/phosphosite', methods=['GET', 'POST'])
def phosphosite_analysis():
    """Phosphosite structural analysis page with improved error handling and supplementary data."""
    import os
    from flask import request, render_template
    
    # Initialize variables
    results = None
    error = None
    
    if request.method == 'POST':
        # Get identifier from form
        identifier = request.form.get('identifier', '')
        id_type = request.form.get('id_type', 'uniprot')
        
        if not identifier:
            return render_template('phosphosite.html', error="Please enter an identifier")
        
        try:
            # Get protein data
            if id_type.lower() == 'uniprot':
                print(f"DEBUG: Getting protein by UniProt ID: {identifier}")
                protein_info = {"uniprot_id": identifier, "gene_symbol": "Unknown", "name": "Unknown Protein"}
                try:
                    protein_data = get_protein_by_id(uniprot_id=identifier)
                    protein_info = {
                        "uniprot_id": identifier,
                        "gene_symbol": protein_data.get("gene_symbol", "Unknown"),
                        "name": protein_data.get("metadata", {}).get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "Unknown Protein")
                    }
                except Exception as e:
                    logger.warning(f"Error getting detailed protein info: {e}")
            else:
                # Get UniProt ID from gene symbol
                try:
                    print(f"DEBUG: Getting protein by gene symbol: {identifier}")
                    pre_filt_id = get_uniprot_id_from_gene(identifier)
                    protein_data = get_protein_by_id(pre_filt_id)
                    protein_info = {
                        "uniprot_id": pre_filt_id,
                        "gene_symbol": identifier,
                        "name": protein_data.get("metadata", {}).get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "Unknown Protein")
                    }
                except Exception as e:
                    logger.error(f"Error getting protein from gene symbol: {e}")
                    return render_template('phosphosite.html', error=f"Could not find protein with gene symbol: {identifier}")
            
            # Get sequence data
            sequence = None
            try:
                full_protein_data = get_protein_by_id(uniprot_id=protein_info["uniprot_id"])
                sequence = full_protein_data.get('metadata', {}).get('sequence', {}).get('value')
                print(f"DEBUG: Got protein sequence of length: {len(sequence) if sequence else 'None'}")
            except Exception as e:
                logger.warning(f"Error getting protein sequence: {e}")
            
            # Get structure for analysis if needed
            structure = None
            try:
                # Special handling for TP53/P04637
                if protein_info["uniprot_id"] in ['P04637', 'P53_HUMAN']:
                    import requests
                    url = f"https://alphafold.ebi.ac.uk/files/AF-{protein_info['uniprot_id']}-F1-model_v4.pdb"
                    response = requests.get(url)
                    if response.status_code == 200:
                        structure = response.text
                        print(f"DEBUG: Got structure directly, length: {len(structure)}")
                    else:
                        # Fall back to normal method
                        structure = get_alphafold_structure(protein_info["uniprot_id"])
                else:
                    structure = get_alphafold_structure(protein_info["uniprot_id"])
                print(f"DEBUG: Got AlphaFold structure")
            except Exception as e:
                logger.warning(f"Error getting AlphaFold structure: {e}")
            
            # Get phosphosites from database
            phosphosites = get_all_phosphosites(protein_info["uniprot_id"])
            print(f"DEBUG: Got {len(phosphosites) if phosphosites else 0} phosphosites from database")
            
            # If we have sequence but no sites, or sites array is empty, analyze all STY sites
            if sequence and (not phosphosites or len(phosphosites) == 0):
                print(f"DEBUG: No phosphosites from database, analyzing all STY sites from sequence")
                # Find all STY sites in the sequence
                from protein_explorer.analysis.OnTheFly import find_all_sty_sites
                phosphosites = find_all_sty_sites(sequence)
                print(f"DEBUG: Found {len(phosphosites)} STY sites in sequence")
                
                # If we have structure, do full analysis
                if structure:
                    try:
                        print(f"DEBUG: Analyzing phosphosites from structure")
                        analyzed_sites = analyze_phosphosites(sequence, structure, protein_info["uniprot_id"])
                        
                        # Create lookup dictionaries
                        analyzed_lookup = {site['site']: site for site in analyzed_sites if 'site' in site}
                        
                        # Update phosphosites with analysis results
                        for site in phosphosites:
                            site_key = site['site']
                            if site_key in analyzed_lookup:
                                analyzed_site = analyzed_lookup[site_key]
                                for key, value in analyzed_site.items():
                                    if value is not None and (key not in site or site[key] is None):
                                        # Round numeric values if needed
                                        if key in ['mean_plddt', 'site_plddt', 'surface_accessibility']:
                                            try:
                                                site[key] = round(float(value), 2)
                                            except (ValueError, TypeError):
                                                site[key] = value
                                        else:
                                            site[key] = value
                    except Exception as e:
                        logger.error(f"Error analyzing phosphosites from structure: {e}")
            elif sequence and phosphosites:
                # We have both sequence and sites from database
                # Make sure we have all STY sites by merging database sites with sequence sites
                from protein_explorer.analysis.OnTheFly import find_all_sty_sites
                phosphosites = find_all_sty_sites(sequence, phosphosites)
                print(f"DEBUG: After merging, have {len(phosphosites)} total STY sites")
            
            # Process phosphosites for correct boolean values
            for site in phosphosites:
                # Handle is_known properly
                if 'is_known_phosphosite' in site:
                    # Convert to integer and then to boolean (0 becomes False, any other number becomes True)
                    try:
                        # First try to convert to float, then int (handles strings like "1.0")
                        is_known_val = int(float(site['is_known_phosphosite']))
                        site['is_known'] = bool(is_known_val)
                        print(f"DEBUG: Setting is_known={site['is_known']} for {site.get('site', 'unknown')}, raw value was {site['is_known_phosphosite']}")
                    except (ValueError, TypeError):
                        # If it can't be converted to int, try direct boolean conversion
                        try:
                            site['is_known'] = bool(site['is_known_phosphosite'])
                            print(f"DEBUG: Setting is_known={site['is_known']} for {site.get('site', 'unknown')}, using bool conversion on {site['is_known_phosphosite']}")
                        except (ValueError, TypeError):
                            # Last resort, default to False
                            site['is_known'] = False
                elif 'is_known' not in site:
                    # Default to False if field is missing
                    site['is_known'] = False
                
                # Handle StructuralSimAvailable properly
                if 'StructuralSimAvailable' in site:
                    # Convert to integer and then to boolean
                    try:
                        # First try to convert to float, then int (handles strings like "1.0")
                        struct_sim_val = int(float(site['StructuralSimAvailable']))
                        site['StructuralSimAvailable'] = bool(struct_sim_val)
                        print(f"DEBUG: Setting StructuralSimAvailable={site['StructuralSimAvailable']} for {site.get('site', 'unknown')}, raw value was {site['StructuralSimAvailable']}")
                    except (ValueError, TypeError):
                        # If it can't be converted to int, try direct boolean conversion
                        try:
                            site['StructuralSimAvailable'] = bool(site['StructuralSimAvailable'])
                            print(f"DEBUG: Setting StructuralSimAvailable={site['StructuralSimAvailable']} for {site.get('site', 'unknown')}, using bool conversion on {site['StructuralSimAvailable']}")
                        except (ValueError, TypeError):
                            # Last resort, default to False
                            site['StructuralSimAvailable'] = False
                else:
                    # Default to False if field is missing
                    site['StructuralSimAvailable'] = False
                    
                # Round numeric fields to 2 decimal places
                for field in ['mean_plddt', 'site_plddt', 'surface_accessibility']:
                    if field in site and site[field] not in [None, 'N/A', '']:
                        try:
                            site[field] = round(float(site[field]), 2)
                        except (ValueError, TypeError):
                            # Skip if not convertible to float
                            pass
            
            # Process sites without structural data on-the-fly
            if structure and sequence:
                try:
                    # Import and use the on-the-fly processor
                    from protein_explorer.analysis.OnTheFly import process_phosphosites_OTF
                    phosphosites = process_phosphosites_OTF(structure, sequence, phosphosites)
                    print(f"DEBUG: Processed phosphosites on-the-fly to fill in missing structural data")
                    
                    # Round all numeric values after on-the-fly processing
                    for site in phosphosites:
                        for field in ['mean_plddt', 'site_plddt', 'surface_accessibility']:
                            if field in site and site[field] not in [None, 'N/A', '']:
                                try:
                                    site[field] = round(float(site[field]), 2)
                                except (ValueError, TypeError):
                                    # Skip if not convertible to float
                                    pass
                except ImportError:
                    print(f"DEBUG: OnTheFly.process_phosphosites_OTF not found, skipping on-the-fly processing")
                except Exception as e:
                    print(f"DEBUG: Error in on-the-fly processing: {e}")
                    import traceback
                    print(traceback.format_exc())
            
            # Debug the processed phosphosite data
            from protein_explorer.utils.debug_helper import debug_phosphosite_data
            debug_phosphosite_data(phosphosites)
            
            # Get structural matches for phosphosites
            structural_matches = {}
            if phosphosites:
                try:
                    structural_matches = find_structural_matches(
                        protein_info["uniprot_id"], 
                        phosphosites
                    )
                    
                    # Enhance with supplementary data
                    for site, matches in structural_matches.items():
                        if matches:
                            structural_matches[site] = enhance_structural_matches(matches, site)
                except Exception as e:
                    logger.error(f"Error finding structural matches: {e}")
            
            # Create enhanced HTML for the phosphosites
            phosphosites_html = None
            try:
                from protein_explorer.analysis.enhanced_table import enhance_phosphosite_table
                phosphosites_html = enhance_phosphosite_table(
                    phosphosites, 
                    protein_info["uniprot_id"]
                )
            except Exception as e:
                logger.error(f"Error creating enhanced table: {e}")
                # Fallback to basic HTML table if enhanced table fails
                phosphosites_html = create_fallback_table_html(phosphosites, protein_info["uniprot_id"])
                
            comparison_container_html = '''
            <div id="phosphosite-comparison-container" class="comparison-container mb-4">
                <!-- Phosphosite comparison visualizations will be loaded here -->
                <div class="text-center py-4">
                    <div class="spinner-border text-primary" role="status">
                        <span class="visually-hidden">Loading...</span>
                    </div>
                    <p class="mt-2">Loading comparison visualizations...</p>
                </div>
            </div>
            '''
            
            # Prepare data for comparison visualizations
            known_unknown_stats = {
                'known_count': sum(1 for site in phosphosites if site.get('is_known', False)),
                'unknown_count': sum(1 for site in phosphosites if not site.get('is_known', False)),
                'total_count': len(phosphosites)
            }
            # Log phosphosite stats for debugging
            print(f"DEBUG: Phosphosite stats - Known: {known_unknown_stats['known_count']}, Unknown: {known_unknown_stats['unknown_count']}, Total: {known_unknown_stats['total_count']}")
            
            # Add data attributes to table rows for phosphosite-comparison.js script
            # This is handled by the enhance_phosphosite_table function or fallback_table_html function
            
            # Additional scripts to load comparison script and initialize
            additional_scripts = '''
            <!-- Ensure Chart.js is loaded before initializing comparison visualization -->
            <script>
                document.addEventListener('DOMContentLoaded', function() {
                    // Log the loading status of scripts
                    console.log('DOM loaded, initializing phosphosite comparison...');
                    console.log('Chart.js loaded:', typeof Chart !== 'undefined');
                    console.log('initPhosphositeComparison function exists:', typeof window.initPhosphositeComparison !== 'undefined');
                    
                    // If Chart.js is not loaded, we need to make sure it loads before our scripts
                    if (typeof Chart === 'undefined') {
                        console.warn('Chart.js not loaded yet, will try to load it again');
                        const script = document.createElement('script');
                        script.src = 'https://cdn.jsdelivr.net/npm/chart.js@3.9.1/dist/chart.min.js';
                        script.onload = function() {
                            console.log('Chart.js loaded successfully on retry');
                            
                            // Load the annotation plugin
                            const annotationPlugin = document.createElement('script');
                            annotationPlugin.src = 'https://cdn.jsdelivr.net/npm/chartjs-plugin-annotation@2.0.1/dist/chartjs-plugin-annotation.min.js';
                            annotationPlugin.onload = function() {
                                // Once Chart.js and plugins are loaded, load visualization scripts
                                loadVisualizationScripts();
                            };
                            document.head.appendChild(annotationPlugin);
                        };
                        document.head.appendChild(script);
                    } else {
                        // Chart.js is already loaded, just check for our scripts
                        if (typeof window.initPhosphositeComparison !== 'function') {
                            console.warn('initPhosphositeComparison function not found, loading script again');
                            loadVisualizationScripts();
                        } else {
                            // Everything is loaded, initialize the comparison view
                            setTimeout(function() {
                                console.log('Initializing phosphosite comparison with delay');
                                initPhosphositeComparison('phosphosite-comparison-container');
                            }, 500); // Small delay to ensure DOM is fully ready
                        }
                    }
                    
                    function loadVisualizationScripts() {
                        // Load phosphosite-visualization.js first
                        const vizScript = document.createElement('script');
                        vizScript.src = "/static/js/phosphosite-visualization.js";
                        vizScript.onload = function() {
                            console.log('phosphosite-visualization.js loaded');
                            
                            // Then load phosphosite-comparison.js
                            const comparisonScript = document.createElement('script');
                            comparisonScript.src = "/static/js/phosphosite-comparison.js";
                            comparisonScript.onload = function() {
                                console.log('phosphosite-comparison.js loaded');
                                console.log('initPhosphositeComparison function now exists:', typeof window.initPhosphositeComparison === 'function');
                                // Initialize the comparison visualization
                                if (typeof window.initPhosphositeComparison === 'function') {
                                    setTimeout(function() {
                                        initPhosphositeComparison('phosphosite-comparison-container');
                                    }, 500); // Small delay to ensure DOM is fully ready
                                }
                            };
                            document.head.appendChild(comparisonScript);
                        };
                        document.head.appendChild(vizScript);
                    }
                });
            </script>
            '''
            
            # Compile results
            results = {
                'protein_info': protein_info,
                'phosphosites': phosphosites,
                'phosphosites_html': phosphosites_html,
                'structural_matches': structural_matches,
                'comparison_container_html': comparison_container_html,
                'known_unknown_stats': known_unknown_stats,
                'additional_scripts': additional_scripts,
                'error': None
            }
            
            # Return results
            return render_template('phosphosite.html', 
                                  protein_info=results['protein_info'],
                                  phosphosites=results['phosphosites'],
                                  phosphosites_html=results.get('phosphosites_html'),
                                  structural_matches=results['structural_matches'],
                                  comparison_container_html=results.get('comparison_container_html'),
                                  known_unknown_stats=results.get('known_unknown_stats'),
                                  additional_scripts=results.get('additional_scripts'),
                                  error=results.get('error'))
                
        except Exception as e:
            logger.error(f"Error in phosphosite analysis: {e}")
            import traceback
            logger.error(traceback.format_exc())
            error = str(e)
            return render_template('phosphosite.html', error=error)
    
    # GET request - show empty form
    return render_template('phosphosite.html')

# Helper function to create a fallback HTML table if enhanced_table fails
def create_fallback_table_html(phosphosites, protein_uniprot_id):
    """Create a basic HTML table for phosphosites if the enhanced table function fails."""
    if not phosphosites:
        return "<div class='alert alert-warning'>No phosphosite data available.</div>"
    
    html = """
    <div class="card mt-4">
        <div class="card-header">
            <h5 class="mb-0">Phosphorylation Site Analysis</h5>
            <small class="ms-4 text-muted" style="font-size: 0.85rem;">
                <!-- Green legend box -->
                <span style="background-color: #c8e6c9; display: inline-block; width: 15px; height: 15px; margin-right: 5px; border: 1px solid #bbb;"></span>
                Has Structural Similarity Data
                &nbsp;&nbsp;
                <!-- Orange legend box -->
                <span style="background-color: #ffcc80; display: inline-block; width: 15px; height: 15px; margin-right: 5px; border: 1px solid #bbb;"></span>
                No Structural Similarity Data
            </small>
        </div>
        <div class="card-body p-0">
            <div class="table-responsive">
                <table class="table table-striped table-hover phosphosite-table">
                    <thead class="thead-light">
                        <tr>
                            <th>Site</th>
                            <th>Motif (-7 to +7)</th>
                            <th>Mean pLDDT</th>
                            <th>Site pLDDT</th>
                            <th>Nearby Residues (10Å)</th>
                            <th>Surface Access.</th>
                            <th>Known</th>
                            <th>Actions</th>
                        </tr>
                    </thead>
                    <tbody id="phosphosite-table">
    """
    
    for site in phosphosites:
        # Determine row style based on structural data availability
        row_style = ""
        if site.get('StructuralSimAvailable', False):
            row_style = 'style="background-color: #c8e6c9;"'  # Green
        else:
            row_style = 'style="background-color: #ffcc80;"'  # Orange
        
        # Format numeric values for display
        try:
            mean_plddt = float(site.get('mean_plddt', 0))
            mean_plddt_text = f"{mean_plddt:.2f}" if mean_plddt > 0 else "N/A"
        except (ValueError, TypeError):
            mean_plddt_text = site.get('mean_plddt', 'N/A')
        
        try:
            site_plddt = float(site.get('site_plddt', 0))
            site_plddt_text = f"{site_plddt:.2f}" if site_plddt > 0 else "N/A"
        except (ValueError, TypeError):
            site_plddt_text = site.get('site_plddt', 'N/A')
        
        try:
            surface_access = float(site.get('surface_accessibility', 0))
            surface_access_text = f"{surface_access:.2f}%" if surface_access > 0 else "N/A"
        except (ValueError, TypeError):
            surface_access_text = site.get('surface_accessibility', 'N/A')
        
        # Get values for data attributes
        resno = site.get('resno', 0)
        site_type = site.get('siteType', '')
        nearby_count = site.get('nearby_count', 0)
        mean_plddt_val = site.get('mean_plddt', 0)
        surface_access_val = site.get('surface_accessibility', 0)
        is_known = str(site.get('is_known', False)).lower()
        
        # Add data attributes to the row for better visualization
        data_attrs = f"""
            data-site="{site.get('site', '')}" 
            data-resno="{resno}" 
            data-type="{site_type}"
            data-plddt="{mean_plddt_val}" 
            data-nearby="{nearby_count}"
            data-surface="{surface_access_val}"
            data-known="{is_known}"
        """
        
        html += f"""
        <tr {row_style} {data_attrs}>
            <td><a href="/site/{protein_uniprot_id}/{site.get('site', '')}" class="site-link" data-resno="{resno}">{site.get('site', '')}</a></td>
            <td><code class="motif-sequence">{site.get('motif', '')}</code></td>
            <td>{mean_plddt_text}</td>
            <td>{site_plddt_text}</td>
            <td>{nearby_count}</td>
            <td>{surface_access_text}</td>
            <td>{"Yes" if site.get('is_known', False) else "No"}</td>
            <td>
                <a href="/site/{protein_uniprot_id}/{site.get('site', '')}" class="btn btn-sm btn-outline-primary">
                    Details
                </a>
            </td>
        </tr>
        """
    
    html += """
                    </tbody>
                </table>
            </div>
        </div>
    </div>
    
    <script>
        document.addEventListener('DOMContentLoaded', function() {
            // Add click handlers to site links
            const siteLinks = document.querySelectorAll('.site-link');
            siteLinks.forEach(link => {
                link.addEventListener('click', function(e) {
                    e.preventDefault();
                    const resno = this.getAttribute('data-resno');
                    
                    // Find the span in the sequence viewer
                    const sequenceSpans = document.querySelectorAll('.sequence-viewer span');
                    if (sequenceSpans.length > 0) {
                        // Find and click the span for this residue
                        const index = parseInt(resno) - 1;
                        if (index >= 0 && index < sequenceSpans.length) {
                            sequenceSpans[index].click();
                        }
                    }
                });
            });
        });
    </script>
    """
    
    return html

# Add this to app.py
@app.route('/api/sequence_matches/<site_id>', methods=['GET'])
def api_sequence_matches(site_id):
    """API endpoint for sequence similarity matches."""
    try:
        # Get query parameters
        top_n = request.args.get('top_n', default=100, type=int)
        min_similarity = request.args.get('min_similarity', default=0.4, type=float)
        
        # Validate site_id format (e.g., UniProtID_ResidueNumber)
        if '_' not in site_id:
            return jsonify({'error': 'Invalid site ID format. Expected: UniProtID_ResidueNumber'}), 400
            
        # Get sequence matches
        matches = find_sequence_matches(site_id, top_n=top_n, min_similarity=min_similarity)
        
        # If no matches, return empty list with a message
        if not matches:
            return jsonify({
                'matches': [],
                'count': 0,
                'message': f'No sequence similarity matches found for {site_id}'
            })
        
        # Create a complete response with metadata
        response = {
            'site_id': site_id,
            'matches': matches,
            'count': len(matches),
            'params': {
                'top_n': top_n,
                'min_similarity': min_similarity
            }
        }
        
        return jsonify(response)
    except Exception as e:
        logger.error(f"Error in sequence matches API: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/sequence_conservation/<site_id>', methods=['GET'])
def api_sequence_conservation(site_id):
    """API endpoint for sequence conservation analysis."""
    try:
        # Get query parameters
        top_n = request.args.get('top_n', default=200, type=int)
        min_similarity = request.args.get('min_similarity', default=0.4, type=float)
        
        # Get the query motif if possible
        query_motif = None
        
        # Parse site_id to extract uniprot_id and site_number
        parts = site_id.split('_')
        if len(parts) >= 2:
            uniprot_id = parts[0]
            site_number = int(parts[1])
            
            # Try to get site data with motif
            try:
                # Get phosphosite data directly using DB integration
                site_data = get_phosphosite_data(site_id)
                if site_data and 'SITE_+/-7_AA' in site_data:
                    query_motif = site_data['SITE_+/-7_AA']
                elif site_data and 'motif' in site_data:
                    query_motif = site_data['motif']
                
                # If not found in DB, try to get from all sites
                if not query_motif:
                    # Get all sites for the protein
                    all_sites = get_all_phosphosites(uniprot_id)
                    
                    # Find the specific site
                    for site in all_sites:
                        if site.get('resno') == site_number:
                            query_motif = site.get('motif')
                            break
            except Exception as e:
                logger.warning(f"Error retrieving motif for {site_id}: {e}")
        
        # Get sequence matches
        matches = find_sequence_matches(site_id, top_n=top_n, min_similarity=min_similarity)
        
        # Analyze conservation
        conservation = analyze_motif_conservation(matches, query_motif=query_motif)
        
        # Add metadata to response
        response = {
            'site_id': site_id,
            'analysis': conservation,
            'match_count': len(matches),
            'query_motif': query_motif
        }
        
        return jsonify(response)
    except Exception as e:
        logger.error(f"Error in sequence conservation API: {e}")
        return jsonify({'error': str(e)}), 500


@app.route('/site/<uniprot_id>/<site>')
def site_detail(uniprot_id, site):
    """Display detailed information about a specific phosphorylation site with enhanced supplementary data."""
    try:
        # Get protein data
        protein_data = get_protein_by_id(uniprot_id=uniprot_id)
        
        # Get the sequence from metadata
        sequence = protein_data.get('metadata', {}).get('sequence', {}).get('value')
        if not sequence:
            return render_template('error.html', error=f"Protein sequence not found for {uniprot_id}")
            
        # Get structure
        structure = get_alphafold_structure(uniprot_id)
        if not structure:
            return render_template('error.html', error=f"Protein structure not found for {uniprot_id}")
        
        # Parse the site string to get type and residue number
        site_match = re.match(r'([A-Z])(\d+)', site)
        if not site_match:
            return render_template('error.html', error=f"Invalid site format: {site}")
            
        site_type = site_match.group(1)
        site_number = int(site_match.group(2))
        
        # Generate site_id in the format UniProtID_ResidueNumber
        site_id = f"{uniprot_id}_{site_number}"
        
        # Get supplementary data for this site
        site_data = get_phosphosite_data(site_id)
        
        # If site_data is None, create a minimal site data dictionary
        if not site_data:
            site_data = {
                'site_id': site_id,
                'site': site,
                'resno': site_number,
                'siteType': site_type
            }
        
        # Get all phosphosites to find the specific site
        all_phosphosites = get_all_phosphosites(uniprot_id)
        
        # If no phosphosites from DB, analyze from structure
        if not all_phosphosites:
            all_phosphosites = analyze_phosphosites(sequence, structure, uniprot_id)
        
        # Find the specific site
        site_data_from_analysis = next((s for s in all_phosphosites if s.get('site') == site), None)
        if site_data_from_analysis:
            # Merge with existing site_data
            for key, value in site_data_from_analysis.items():
                if key not in site_data or site_data[key] is None:
                    site_data[key] = value
        
        # Find structural matches
        structural_matches = []
        try:
            # Create basic site data for structural match search if needed
            basic_site_data = {
                'site': site,
                'resno': site_number,
                'siteType': site_type
            }
            
            # Get structural matches with database integration
            all_matches = find_structural_matches(uniprot_id, [basic_site_data], top_n=None)
            raw_matches = all_matches.get(site, [])
            
            # Filter out self-matches
            raw_matches = [match for match in raw_matches if match.get('rmsd', 0) > 0.01]
            
            # Enhance matches with supplementary data
            structural_matches = enhance_structural_matches(raw_matches, site)
            
            # Sort matches by RMSD
            structural_matches = sorted(structural_matches, key=lambda x: x.get('rmsd', float('inf')))
            
        except Exception as e:
            logger.error(f"Error finding structural matches: {e}")
            structural_matches = []
        
        # Create basic structure visualization
        structure_html = pe.visualization.visualize_structure(
            structure,
            sequence=sequence
        )
        
        # Create specialized visualizations
        # Network visualization
        try:
            network_html = create_phosphosite_network(site, structural_matches, site_data)
        except Exception as e:
            logger.error(f"Error creating network visualization: {e}")
            network_html = f"<div class='alert alert-warning'>Error creating network visualization: {str(e)}</div>"
            
        # Comparative motif visualization
        try:
            if structural_matches and 'motif' in site_data:
                motif_html = create_comparative_motif_visualization(site_data, structural_matches)
            else:
                motif_html = "<div class='alert alert-info'>No motif data available for comparison.</div>"
        except Exception as e:
            logger.error(f"Error creating motif visualization: {e}")
            motif_html = f"<div class='alert alert-warning'>Error creating motif visualization: {str(e)}</div>"
            
        # Residue distribution analysis
        try:
            if structural_matches:
                distribution_data = analyze_residue_distributions(structural_matches)
            else:
                distribution_data = None
        except Exception as e:
            logger.error(f"Error analyzing residue distributions: {e}")
            distribution_data = None
            
        # Enhanced 3D visualization
        try:
            enhanced_3d_html = enhance_site_visualization(uniprot_id, site, site_data)
        except Exception as e:
            logger.error(f"Error creating enhanced 3D visualization: {e}")
            enhanced_3d_html = f"<div class='alert alert-warning'>Error creating enhanced 3D visualization: {str(e)}</div>"
            
        # Structural context analysis
        try:
            context_data = analyze_phosphosite_context(structure, site_number, site_type)
        except Exception as e:
            logger.error(f"Error analyzing structural context: {e}")
            context_data = None
        
        # Find sequence similarity matches
        sequence_matches = []
        sequence_network_data = None
        sequence_conservation = None
        sequence_motif_html = None
        
        try:
            # Get sequence similarity matches using DB integration
            logger.info(f"Finding sequence matches for {site_id}")
            
            sequence_matches = find_sequence_matches(site_id, top_n=200, min_similarity=0.4)
            logger.info(f"Found {len(sequence_matches)} sequence matches")
            
            # Log a sample of the matches to see what they contain
            if sequence_matches:
                logger.info(f"Sample match: {sequence_matches[0]}")
                
                # Create sequence network data
                sequence_network_data = create_sequence_network_data(
                    site_id, 
                    sequence_matches,
                    query_motif=site_data.get('motif')
                )
                logger.info(f"Created network data with {len(sequence_network_data['nodes'])} nodes")
                
                # Create sequence conservation analysis
                sequence_conservation = analyze_motif_conservation(
                    sequence_matches,
                    query_motif=site_data.get('motif')
                )
                
                if sequence_conservation and 'motif_count' in sequence_conservation:
                    logger.info(f"Created conservation analysis with {sequence_conservation['motif_count']} motifs")
                
                # Create sequence motif visualization
                sequence_motif_html = create_sequence_motif_visualization(
                    site_id,
                    site_data.get('motif', ''),
                    sequence_matches
                )
                logger.info(f"Created motif visualization HTML of length {len(sequence_motif_html) if sequence_motif_html else 0}")
            else:
                logger.warning(f"No sequence matches found for {site_id}")
        except Exception as e:
            logger.error(f"Error analyzing sequence similarity: {e}")
            import traceback
            logger.error(traceback.format_exc())
            sequence_matches = []
        
        # Standard kinase prediction data
        structure_kinase_data = {}
        sequence_kinase_data = {}
        
        # Network-based kinase prediction data
        structure_network_kinase_data = {}
        sequence_network_kinase_data = {}
        
        try:
            # Standard structural kinase prediction
            struct_scores = get_site_kinase_scores(site_id, 'structure')
            if struct_scores and 'scores' in struct_scores:
                # Get top kinases
                top_kinases = predict_kinases(site_id, top_n=10, score_type='structure')
                
                # Get kinase families
                kinase_families = categorize_kinases_by_family(top_kinases)
                
                # Get known kinase info
                known_kinase = get_known_kinase_info(site_id, 'structure')
                
                # Prepare data for heatmap visualization
                struct_match_ids = [f"{match['target_uniprot']}_{match['target_site'].replace('S', '').replace('T', '').replace('Y', '')}" 
                                   for match in structural_matches[:20] if match.get('rmsd', 10) < 5.0]
                struct_match_ids.insert(0, site_id)  # Add the query site first
                
                heatmap_data = get_heatmap_data(struct_match_ids, top_n=10, score_type='structure')
                
                # Bundle all data
                structure_kinase_data = {
                    'site_id': site_id,
                    'top_kinases': top_kinases,
                    'kinase_families': kinase_families,
                    'known_kinase': known_kinase,
                    'heatmap': heatmap_data
                }
            
            # Standard sequence kinase prediction
            seq_scores = get_site_kinase_scores(site_id, 'sequence')
            if seq_scores and 'scores' in seq_scores:
                # Get top kinases
                top_kinases = predict_kinases(site_id, top_n=10, score_type='sequence')
                
                # Get kinase families
                kinase_families = categorize_kinases_by_family(top_kinases)
                
                # Get known kinase info
                known_kinase = get_known_kinase_info(site_id, 'sequence')
                
                # Prepare data for heatmap visualization
                seq_match_ids = []
                if sequence_matches:
                    seq_match_ids = [match['target_id'] for match in sequence_matches[:20] 
                                    if match.get('similarity', 0) > 0.6]
                    seq_match_ids.insert(0, site_id)  # Add the query site first
                
                heatmap_data = get_heatmap_data(seq_match_ids, top_n=10, score_type='sequence')
                
                # Bundle all data
                sequence_kinase_data = {
                    'site_id': site_id,
                    'top_kinases': top_kinases,
                    'kinase_families': kinase_families,
                    'known_kinase': known_kinase,
                    'heatmap': heatmap_data
                }
            
            # Network-based kinase predictions
            try:
                # 1. Structure network prediction
                structure_network_kinases = predict_kinases_network(
                    site_id, top_n=10, score_type='structure',
                    similarity_threshold=0.6, rmsd_threshold=3.0
                )
                
                if structure_network_kinases:
                    # Get similar sites
                    struct_similar_sites = get_similar_sites(
                        site_id,
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Get heatmap data with enhanced function
                    struct_heatmap_data = get_network_heatmap_data(
                        site_id, top_n=10, score_type='structure',
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Log heatmap data for debugging
                    logger.info(f"Structure heatmap data: {len(struct_heatmap_data['sites'])} sites, " 
                                f"{len(struct_heatmap_data['kinases'])} kinases, "
                                f"{len(struct_heatmap_data['scores'])} scores")
                    
                    # Get family distribution
                    struct_family_dist = get_kinase_family_distribution_network(
                        site_id, score_type='structure',
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Combine into data structure
                    structure_network_kinase_data = {
                        'site_id': site_id,
                        'top_kinases': structure_network_kinases,
                        'heatmap': struct_heatmap_data,
                        'kinase_families': struct_family_dist,
                        'site_count': len(struct_similar_sites),
                        'rmsd_threshold': 3.0,
                        'similarity_threshold': 0.6
                    }
                
                # 2. Sequence network prediction
                sequence_network_kinases = predict_kinases_network(
                    site_id, top_n=10, score_type='sequence',
                    similarity_threshold=0.6, rmsd_threshold=3.0
                )
                
                if sequence_network_kinases:
                    # Get similar sites
                    seq_similar_sites = get_similar_sites(
                        site_id,
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Get heatmap data with enhanced function
                    seq_heatmap_data = get_network_heatmap_data(
                        site_id, top_n=10, score_type='sequence',
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Log heatmap data for debugging
                    logger.info(f"Sequence heatmap data: {len(seq_heatmap_data['sites'])} sites, " 
                                f"{len(seq_heatmap_data['kinases'])} kinases, "
                                f"{len(seq_heatmap_data['scores'])} scores")
                    
                    # Get family distribution
                    seq_family_dist = get_kinase_family_distribution_network(
                        site_id, score_type='sequence',
                        similarity_threshold=0.6, rmsd_threshold=3.0
                    )
                    
                    # Combine into data structure
                    sequence_network_kinase_data = {
                        'site_id': site_id,
                        'top_kinases': sequence_network_kinases,
                        'heatmap': seq_heatmap_data,
                        'kinase_families': seq_family_dist,
                        'site_count': len(seq_similar_sites),
                        'rmsd_threshold': 3.0,
                        'similarity_threshold': 0.6
                    }
            except Exception as network_e:
                logger.error(f"Error in network kinase prediction: {network_e}")
                import traceback
                logger.error(traceback.format_exc())
                
                # Provide minimal fallback data to avoid template errors
                structure_network_kinase_data = {
                    'site_id': site_id,
                    'top_kinases': [],
                    'heatmap': {'sites': [], 'kinases': [], 'scores': []},
                    'kinase_families': {},
                    'site_count': 0,
                    'rmsd_threshold': 3.0,
                    'similarity_threshold': 0.6,
                    'error': str(network_e)
                }
                
                sequence_network_kinase_data = {
                    'site_id': site_id,
                    'top_kinases': [],
                    'heatmap': {'sites': [], 'kinases': [], 'scores': []},
                    'kinase_families': {},
                    'site_count': 0,
                    'rmsd_threshold': 3.0,
                    'similarity_threshold': 0.6,
                    'error': str(network_e)
                }

        except Exception as e:
            logger.error(f"Error processing kinase data: {e}")
            import traceback
            logger.error(traceback.format_exc())
        
        # Render the template with all the data
        return render_template(
            'site.html',
            protein=protein_data,
            site=site,
            site_data=site_data,
            site_id=site_id,  # Make sure site_id is passed to the template
            structure_html=structure_html,
            structural_matches=structural_matches,
            supplementary_data=site_data,
            network_html=network_html,
            motif_html=motif_html,
            distribution_data=distribution_data,
            enhanced_3d_html=enhanced_3d_html,
            context_data=context_data,
            # Sequence similarity data:
            sequence_matches=sequence_matches,
            sequence_network_data=sequence_network_data,
            sequence_conservation=sequence_conservation,
            sequence_motif_html=sequence_motif_html,
            # Standard kinase prediction data:
            structure_kinase_data=structure_kinase_data,
            sequence_kinase_data=sequence_kinase_data,
            # Network kinase prediction data:
            structure_network_kinase_data=structure_network_kinase_data,
            sequence_network_kinase_data=sequence_network_kinase_data
        )
    except Exception as e:
        logger.error(f"Error in site detail view: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return render_template('error.html', error=str(e))
    

def analyze_motif(motif: str, site_type: str, site_number: int) -> Dict:
    """
    Analyze a phosphosite motif sequence for additional insights.
    
    Args:
        motif: The motif sequence string
        site_type: The site type (S, T, or Y)
        site_number: The residue number
        
    Returns:
        Dictionary with motif analysis
    """
    # Find the position of the phosphosite in the motif
    motif_length = len(motif)
    
    # Identify the center residue position in the motif
    center_pos = None
    for i, aa in enumerate(motif):
        # Look for site_type (S, T, or Y) in appropriate position
        if aa == site_type:
            # Calculate expected position from start of motif
            expected_pos = site_number - 8  # If site is the center of a -7 to +7 motif
            if i == 7:  # Expected position for -7 to +7 motif (0-indexed)
                center_pos = i
                break
    
    if center_pos is None:
        # If we can't find it exactly, assume it's in the middle
        center_pos = motif_length // 2
    
    # Count amino acid types
    aa_groups = {
        'polar': 'STYCNQ',
        'nonpolar': 'AVILMFWPG',
        'acidic': 'DE',
        'basic': 'KRH',
        'other': 'X'
    }
    
    counts = {group: 0 for group in aa_groups}
    upstream_counts = {group: 0 for group in aa_groups}
    downstream_counts = {group: 0 for group in aa_groups}
    
    for i, aa in enumerate(motif):
        # Skip if it's the phosphosite
        if i == center_pos:
            continue
            
        # Determine group
        group = 'other'
        for g, aas in aa_groups.items():
            if aa in aas:
                group = g
                break
                
        # Increment total count
        counts[group] += 1
        
        # Increment upstream/downstream count
        if i < center_pos:
            upstream_counts[group] += 1
        else:
            downstream_counts[group] += 1
    
    # Calculate percentages
    total_aas = motif_length - 1  # Exclude phosphosite
    upstream_total = center_pos
    downstream_total = motif_length - center_pos - 1
    
    percentages = {group: count/total_aas*100 for group, count in counts.items()}
    upstream_percentages = {f"upstream_{group}": count/upstream_total*100 if upstream_total > 0 else 0 
                          for group, count in upstream_counts.items()}
    downstream_percentages = {f"downstream_{group}": count/downstream_total*100 if downstream_total > 0 else 0 
                            for group, count in downstream_counts.items()}
    
    # Extract -3 to +3 submotif (important region for kinase recognition)
    narrow_motif = motif[max(0, center_pos-3):min(motif_length, center_pos+4)]
    narrow_motif_highlighted = f"{narrow_motif[:3]}<strong>{narrow_motif[3:4]}</strong>{narrow_motif[4:]}" if len(narrow_motif) >= 7 else narrow_motif
    
    # Create motif visualization with color coding
    motif_vis = []
    for i, aa in enumerate(motif):
        position = i - center_pos
        # Determine AA type for coloring
        aa_class = "highlighted" if i == center_pos else "other"
        for group, aas in aa_groups.items():
            if aa in aas:
                aa_class = group
                break
                
        motif_vis.append({
            "aa": aa,
            "position": position,
            "class": aa_class
        })
    
    # Check for known kinase recognition motifs
    motif_patterns = [
        {"name": "CDK", "pattern": r"[ST]P.?[KR]", "kinases": ["CDK1", "CDK2", "CDK5"]},
        {"name": "MAPK/ERK", "pattern": r"P.[ST]P", "kinases": ["ERK1", "ERK2", "p38", "JNK"]},
        {"name": "PKA/PKG", "pattern": r"[RK][RK].?[ST]", "kinases": ["PKA", "PKG"]},
        {"name": "CK2", "pattern": r"[ST]..[DE]", "kinases": ["CK2"]},
        {"name": "GSK3", "pattern": r"[ST]..[ST]", "kinases": ["GSK3α", "GSK3β"]},
        {"name": "CK1", "pattern": r"[DE]..[ST].[ST]", "kinases": ["CK1"]},
        {"name": "PLK", "pattern": r"[DE]..[ST]", "kinases": ["PLK1"]},
        {"name": "Akt/PKB", "pattern": r"R.R..[ST]", "kinases": ["Akt1", "Akt2", "Akt3"]},
        {"name": "AGC", "pattern": r"R.[ST]", "kinases": ["PKA", "PKC", "PKG"]}
    ]
    
    potential_kinases = []
    for pattern_info in motif_patterns:
        pattern = pattern_info["pattern"]
        match = re.search(pattern, motif)
        if match:
            potential_kinases.append({
                "pattern_name": pattern_info["name"],
                "pattern": pattern,
                "kinases": pattern_info["kinases"],
                "match": match.group(0)
            })
    
    return {
        "motif_vis": motif_vis,
        "center_position": center_pos,
        "narrow_motif": narrow_motif,
        "narrow_motif_highlighted": narrow_motif_highlighted,
        "aa_counts": counts,
        "aa_percentages": percentages,
        "upstream_percentages": upstream_percentages,
        "downstream_percentages": downstream_percentages,
        "potential_kinases": potential_kinases
    }

@app.route('/site-search', methods=['GET', 'POST'])
def site_search():
    """Search for a specific phosphorylation site."""
    if request.method == 'POST':
        # Get form data
        uniprot_id = request.form.get('uniprot_id', '')
        site = request.form.get('site', '')
        
        if not uniprot_id:
            return render_template('site-search.html', error="Please enter a UniProt ID")
        
        if not site:
            return render_template('site-search.html', error="Please enter a site identifier")
            
        # Validate site format (S, T, or Y followed by a number)
        import re
        if not re.match(r'^[STY]\d+$', site):
            return render_template('site-search.html', error="Invalid site format. Expected format: S123, T45, Y678, etc.")
            
        # Redirect to site detail page
        return redirect(url_for('site_detail', uniprot_id=uniprot_id, site=site))
    
    # GET request - show search form
    return render_template('site-search.html')

# Add a new API endpoint for kinase data
@app.route('/api/kinases/<site_id>', methods=['GET'])
def api_kinases(site_id):
    """API endpoint for kinase prediction scores."""
    try:
        score_type = request.args.get('type', 'structure')
        top_n = int(request.args.get('top_n', 10))
        
        # Get kinase scores using integration function
        scores = get_site_kinase_scores(site_id, score_type)
        if not scores or 'scores' not in scores:
            return jsonify({'error': f'No {score_type} kinase scores found for {site_id}'}), 404
        
        # Get top predicted kinases
        top_kinases = predict_kinases(site_id, top_n, score_type)
        
        # Get known kinase info
        known_kinase = get_known_kinase_info(site_id, score_type)
        
        # Get comparison with the other score type
        other_type = 'sequence' if score_type == 'structure' else 'structure'
        comparison = get_kinase_comparison_data(site_id, [score_type, other_type], top_n)
        
        # Return all data
        return jsonify({
            'site_id': site_id,
            'score_type': score_type,
            'known_kinase': known_kinase,
            'top_kinases': top_kinases,
            'comparison': comparison
        })
    except Exception as e:
        logger.error(f"Error in kinase API: {e}")
        return jsonify({'error': str(e)}), 500

# Add another API endpoint for comparing kinases between sites
@app.route('/api/kinases/compare', methods=['POST'])
def api_compare_kinases():
    """API endpoint for comparing kinase predictions across multiple sites."""
    try:
        data = request.get_json()
        if not data or 'site_ids' not in data:
            return jsonify({'error': 'Missing site_ids parameter'}), 400
        
        site_ids = data['site_ids']
        score_type = data.get('type', 'structure')
        top_n = int(data.get('top_n', 10))
        
        # Get heatmap data
        heatmap_data = get_heatmap_data(site_ids, top_n, score_type)
        
        # Return the data
        return jsonify({
            'heatmap': heatmap_data
        })
    except Exception as e:
        logger.error(f"Error in kinase comparison API: {e}")
        return jsonify({'error': str(e)}), 500
    

# Add a new API endpoint for network-based kinase predictions
@app.route('/api/network-kinases/<site_id>', methods=['GET'])
def api_network_kinases(site_id):
    """API endpoint for network-based kinase prediction scores."""
    try:
        score_type = request.args.get('type', 'structure')
        top_n = int(request.args.get('top_n', 10))
        similarity_threshold = float(request.args.get('similarity_threshold', 0.6))
        rmsd_threshold = float(request.args.get('rmsd_threshold', 3.0))
        
        # Get network-based kinase predictions
        network_kinases = predict_kinases_network(
            site_id, top_n=top_n, score_type=score_type,
            similarity_threshold=similarity_threshold, 
            rmsd_threshold=rmsd_threshold
        )
        
        if not network_kinases:
            return jsonify({'error': f'No {score_type} network kinase scores found for {site_id}'}), 404
        
        # Get similar sites
        similar_sites = get_similar_sites(
            site_id,
            similarity_threshold=similarity_threshold, 
            rmsd_threshold=rmsd_threshold
        )
        
        # Get heatmap data
        heatmap_data = get_network_heatmap_data(
            site_id, top_n=top_n, score_type=score_type,
            similarity_threshold=similarity_threshold, 
            rmsd_threshold=rmsd_threshold
        )
        
        # Get kinase families
        kinase_families = get_kinase_family_distribution_network(
            site_id, score_type=score_type,
            similarity_threshold=similarity_threshold, 
            rmsd_threshold=rmsd_threshold
        )
        
        # Return all data
        return jsonify({
            'site_id': site_id,
            'score_type': score_type,
            'top_kinases': network_kinases,
            'heatmap': heatmap_data,
            'kinase_families': kinase_families,
            'site_count': len(similar_sites),
            'rmsd_threshold': rmsd_threshold,
            'similarity_threshold': similarity_threshold
        })
    except Exception as e:
        logger.error(f"Error in network kinase API: {e}")
        return jsonify({'error': str(e)}), 500

# Add an endpoint to update network predictions with new thresholds
@app.route('/api/update-network-kinases/<site_id>', methods=['POST'])
def update_network_kinases(site_id):
    """Update network-based kinase predictions with new thresholds."""
    try:
        data = request.get_json()
        if not data:
            return jsonify({'error': 'No data provided'}), 400
        
        score_type = data.get('score_type', 'structure')
        top_n = int(data.get('top_n', 10))
        
        # Get appropriate threshold based on score_type
        if score_type == 'structure':
            threshold_param = 'rmsd_threshold'
            threshold_default = 3.0
        else:  # sequence
            threshold_param = 'similarity_threshold'
            threshold_default = 0.6
            
        threshold = float(data.get(threshold_param, threshold_default))
        
        # Set common parameters 
        kwargs = {
            'site_id': site_id,
            'top_n': top_n,
            'score_type': score_type
        }
        
        # Add the appropriate threshold parameter
        if score_type == 'structure':
            kwargs['rmsd_threshold'] = threshold
            kwargs['similarity_threshold'] = 0.6  # default
        else:
            kwargs['similarity_threshold'] = threshold
            kwargs['rmsd_threshold'] = 3.0  # default
        
        # Get updated predictions using integration functions
        network_kinases = predict_kinases_network(**kwargs)
        similar_sites = get_similar_sites(site_id, similarity_threshold=kwargs.get('similarity_threshold'), 
                                         rmsd_threshold=kwargs.get('rmsd_threshold'))
        heatmap_data = get_network_heatmap_data(**kwargs)
        family_distribution = get_kinase_family_distribution_network(**kwargs)
        
        # Return updated data
        return jsonify({
            'site_id': site_id,
            'score_type': score_type,
            'top_kinases': network_kinases,
            'heatmap': heatmap_data,
            'kinase_families': family_distribution,
            'site_count': len(similar_sites),
            'rmsd_threshold': kwargs.get('rmsd_threshold'),
            'similarity_threshold': kwargs.get('similarity_threshold')
        })
    except Exception as e:
        logger.error(f"Error updating network kinase predictions: {e}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/health/database', methods=['GET'])
def database_health():
    """Check database health."""
    health_status = check_database_health()
    return jsonify(health_status)

@app.route('/faq')
def faq():
    """Render the FAQ page."""
    return render_template('faq.html')

if __name__ == '__main__':
    app.run(debug=True)