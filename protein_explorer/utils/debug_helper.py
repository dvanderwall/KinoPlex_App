"""
Debug Helper module - Utilities for debugging phosphosite data.
"""

def debug_phosphosite_data(phosphosites):
    """
    Print debug information about phosphosite data.
    
    Args:
        phosphosites: List of phosphosite dictionaries
    """
    print(f"DEBUG: Found {len(phosphosites)} phosphosites")
    
    # Print information for the first few sites
    for i, site in enumerate(phosphosites[:5]):  # Show only first 5 for brevity
        print(f"DEBUG: Site {i+1}: {site.get('site', 'Unknown')}")
        print(f"  is_known: {site.get('is_known', 'Not set')} (type: {type(site.get('is_known')).__name__})")
        print(f"  is_known_phosphosite: {site.get('is_known_phosphosite', 'Not set')} (type: {type(site.get('is_known_phosphosite')).__name__})")
        print(f"  StructuralSimAvailable: {site.get('StructuralSimAvailable', 'Not set')} (type: {type(site.get('StructuralSimAvailable')).__name__})")
        
        # Check for presence of important fields
        important_fields = ['resno', 'siteType', 'motif', 'mean_plddt', 'site_plddt', 'nearby_count', 'surface_accessibility']
        for field in important_fields:
            print(f"  {field}: {'Present' if field in site else 'Missing'}")
    
    # Count sites that are marked as known
    known_count = sum(1 for site in phosphosites if site.get('is_known', False))
    print(f"DEBUG: {known_count} of {len(phosphosites)} sites are marked as known")
    
    # Count sites that have structural similarity data
    struct_sim_count = sum(1 for site in phosphosites if site.get('StructuralSimAvailable', False))
    print(f"DEBUG: {struct_sim_count} of {len(phosphosites)} sites have structural similarity data")