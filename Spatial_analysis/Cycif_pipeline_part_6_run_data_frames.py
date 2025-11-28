"""
CYCIF PIPELINE PART 6: FREQUENCY ANALYSIS
Version: 18.0 (Spatial metrics + CSV subsets + buffer + 2 compartments)

NOTE: This is a template script. Users must provide their own:
- Sample identifiers
- Group assignments (group_assignments.csv)
- Data paths
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path
from shapely.geometry import shape, Point
from shapely.vectorized import contains
import warnings
warnings.filterwarnings('ignore')
from datetime import datetime

print("=" * 80)
print("CYCIF PIPELINE PART 6: FREQUENCY ANALYSIS")
print("=" * 80)
print(f"Version: 18.0 (Spatial metrics + compartment analysis)")
print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print("=" * 80)

# ============================================================================
# CONFIGURATION - MODIFY THESE PATHS FOR YOUR DATA
# ============================================================================

BASE_DIR = Path(r'path/to/your/data')  # <- MODIFY THIS
OUTPUT_BASE = Path('output/analysis')
OUTPUT_CRYPT_IEL = OUTPUT_BASE / 'Crypt_IEL'
OUTPUT_LP = OUTPUT_BASE / 'LP'

EXCLUDED_SAMPLES = []
GROUPS_FILE = BASE_DIR / 'group_assignments.csv'  # Must contain: sample, group columns
PIXEL_SIZE_UM = 0.325
BUFFER_DISTANCE_PX = 31  # Buffer for compartment assignment

# Cell population definitions
CD45_POPULATION = 'All_Leukocytes'
EPITHELIAL_POPULATION = 'Epithelmarker_Positive'
CD4_SUBSET = 'CD4_T'
CD8_SUBSET = 'CD8_T'
TREG_SUBSET = 'Treg'
TISSUE_RESIDENT_SUBSET = 'Tissue_Resident_T'

# Available cell subsets for analysis
LEUKOCYTE_SUBSETS = [
    'All_Leukocytes', 'T_Cells', 'CD4_T', 'CD8_T', 'DN_T', 'DP_T',
    'Tissue_Resident_T', 'CD8_CD103_CD69_T', 'Treg',
    'CD4_T_naiv', 'CD8_T_naiv', 'TCRgd',
    'B_Cells', 'NK_Cells',
    'Macrophages', 'monocytes', 'Neutrophils', 'mDC',
    'prolif_Leuko'
]

NON_LEUKOCYTE_SUBSETS = [
    'Epithelmarker_Positive', 'prolif_Epi'
]

ALL_SUBSETS = LEUKOCYTE_SUBSETS + NON_LEUKOCYTE_SUBSETS

for dir_path in [OUTPUT_CRYPT_IEL, OUTPUT_LP]:
    dir_path.mkdir(parents=True, exist_ok=True)

print(f"\n✅ Configuration loaded")
print(f"   Total subsets: {len(ALL_SUBSETS)} ({len(LEUKOCYTE_SUBSETS)} leukocyte + {len(NON_LEUKOCYTE_SUBSETS)} non-leukocyte)")

# ============================================================================
# LOAD SAMPLE GROUPS
# ============================================================================

# NOTE: Users must create their own group_assignments.csv with columns:
# - sample: Sample identifier (e.g., 'SAMPLE_001', 'SAMPLE_002')
# - group: Group assignment (e.g., 'GROUP_A', 'GROUP_B')

# Example structure:
# sample,group
# SAMPLE_001,GROUP_A
# SAMPLE_002,GROUP_A
# SAMPLE_003,GROUP_B

groups_df = pd.read_csv(GROUPS_FILE)
sample_to_group = dict(zip(groups_df['sample'].astype(str), groups_df['group']))

# Define your sample list here
SAMPLES = list(groups_df['sample'].astype(str).unique())

print(f"✅ Sample groups loaded: {len(SAMPLES)} samples")

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def load_geojson(geojson_path):
    """Load GeoJSON file and extract crypt features."""
    with open(geojson_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    crypt_features = []
    for feat in data.get('features', []):
        props = feat.get('properties', {})
        classification = props.get('classification', {}).get('name', '').lower()
        name = props.get('name', '').lower()
        if classification == 'crypt' or ('crypt' in name and 'non' not in name):
            crypt_features.append(feat)
    return crypt_features

def calculate_polygon_area(feature):
    """Calculate polygon area in pixels and mm²."""
    poly = shape(feature['geometry'])
    area_pixels = poly.area
    area_mm2 = area_pixels * (PIXEL_SIZE_UM ** 2) / 1e6
    return area_pixels, area_mm2

def assign_cells_to_crypts_no_buffer(cells, crypt_features, x_col='X_centroid', y_col='Y_centroid'):
    """Assign cells to crypts without buffer zone."""
    cells = cells.copy()
    cells['crypt_id'] = ''
    cells['crypt_name'] = ''
    cells['crypt_index'] = -1
    
    for crypt_idx, feat in enumerate(crypt_features):
        poly = shape(feat['geometry'])
        mask = contains(poly, cells[x_col].values, cells[y_col].values)
        uuid = feat.get('id', f'unknown_{crypt_idx}')
        name = feat.get('properties', {}).get('name', f'Crypt_{crypt_idx}')
        unassigned_mask = (cells['crypt_id'] == '') & mask
        cells.loc[unassigned_mask, 'crypt_id'] = uuid
        cells.loc[unassigned_mask, 'crypt_name'] = name
        cells.loc[unassigned_mask, 'crypt_index'] = crypt_idx
    
    n_assigned = (cells['crypt_id'] != '').sum()
    print(f"   No buffer: {n_assigned:,}/{len(cells):,} cells ({n_assigned/len(cells)*100:.1f}%)")
    return cells

def assign_cells_to_crypts_with_buffer(cells, crypt_features, buffer_px, x_col='X_centroid', y_col='Y_centroid'):
    """Assign cells to crypts with buffer zone for IEL capture."""
    cells = cells.copy()
    cells['crypt_id'] = ''
    cells['crypt_name'] = ''
    cells['crypt_index'] = -1
    cells['distance_to_crypt'] = np.inf
    
    for crypt_idx, feat in enumerate(crypt_features):
        poly = shape(feat['geometry'])
        poly_buffered = poly.buffer(buffer_px)
        mask = contains(poly_buffered, cells[x_col].values, cells[y_col].values)
        uuid = feat.get('id', f'unknown_{crypt_idx}')
        name = feat.get('properties', {}).get('name', f'Crypt_{crypt_idx}')
        
        candidate_cells = cells[mask]
        if len(candidate_cells) > 0:
            points = [Point(x, y) for x, y in zip(candidate_cells[x_col], candidate_cells[y_col])]
            distances = np.array([poly.distance(pt) for pt in points])
            for idx, dist in zip(candidate_cells.index, distances):
                if cells.loc[idx, 'crypt_id'] == '' or dist < cells.loc[idx, 'distance_to_crypt']:
                    cells.loc[idx, 'crypt_id'] = uuid
                    cells.loc[idx, 'crypt_name'] = name
                    cells.loc[idx, 'crypt_index'] = crypt_idx
                    cells.loc[idx, 'distance_to_crypt'] = dist
    
    cells = cells.drop(columns=['distance_to_crypt'])
    n_assigned = (cells['crypt_id'] != '').sum()
    print(f"   With buffer ({buffer_px}px): {n_assigned:,}/{len(cells):,} cells ({n_assigned/len(cells)*100:.1f}%)")
    return cells

print("✅ Helper functions defined")

# ============================================================================
# MAIN ANALYSIS - IMPLEMENT YOUR ANALYSIS LOGIC HERE
# ============================================================================

print("\n" + "="*80)
print("NOTE: This is a template. Implement your analysis logic below.")
print("="*80)

# Example analysis workflow:
# 1. Load cell data for each sample
# 2. Load compartment segmentation (GeoJSON)
# 3. Assign cells to compartments
# 4. Calculate frequency metrics
# 5. Export results

print("\n✅ Template ready for customization")
