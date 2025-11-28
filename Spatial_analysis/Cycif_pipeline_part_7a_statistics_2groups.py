"""
CYCIF PIPELINE PART 7A: TWO-GROUP STATISTICS
Compares two experimental groups using non-parametric tests.

NOTE: This is a template script. Users must:
- Define their own group names
- Provide their own data files
- Customize metrics as needed
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path
import numpy as np

# Configuration
OUTPUT_BASE = Path('output/analysis/statistics_2groups')
OUTPUT_BASE.mkdir(parents=True, exist_ok=True)

# Create subfolders
(OUTPUT_BASE / 'CSV').mkdir(exist_ok=True)
(OUTPUT_BASE / 'PNG').mkdir(exist_ok=True)
(OUTPUT_BASE / 'SVG').mkdir(exist_ok=True)

print("="*80)
print("CYCIF PIPELINE PART 7A: TWO-GROUP STATISTICS")
print("="*80)

# ============================================================================
# CONFIGURATION - MODIFY FOR YOUR STUDY
# ============================================================================

# Define your group names (pseudonymized)
GROUP_A = 'GROUP_A'  # e.g., Control group
GROUP_B = 'GROUP_B'  # e.g., Experimental group

# Load data - MODIFY PATHS
# df_lp = pd.read_csv('path/to/lamina_propria_metrics.csv')
# df_crypt = pd.read_csv('path/to/crypt_level_metrics.csv')

# Example data loading (replace with your actual data)
print("\nNOTE: This is a template. Load your own data files.")
print("Required columns: 'sample', 'group', and metric columns")

# Define metrics to analyze
NORMALIZATIONS = ['per_CD45', 'per_Total', 'per_100_Epi', 'per_mm2']
SUBSETS = ['Treg', 'CD4_T', 'CD8_T', 'T_Cells', 'B_Cells', 'NK_Cells', 
           'Macrophages', 'DN_T', 'DP_T', 'Tissue_Resident_T']

def create_bar_plot_with_points(data, metric, tissue, group_a, group_b, output_prefix):
    """
    Create bar plot comparing two groups with individual data points.
    
    Parameters:
    -----------
    data : DataFrame
        Data containing 'group' column and metric values
    metric : str
        Column name of metric to plot
    tissue : str
        Tissue compartment name for title
    group_a, group_b : str
        Names of the two groups to compare
    output_prefix : str
        Prefix for output filename
    """
    
    # Get data for both groups
    data_a = data[data['group'] == group_a][metric].dropna()
    data_b = data[data['group'] == group_b][metric].dropna()
    
    if len(data_a) == 0 or len(data_b) == 0:
        return None
    
    # Statistical test (non-parametric)
    stat, pval = stats.mannwhitneyu(data_a, data_b, alternative='two-sided')
    
    # Create figure
    fig, ax = plt.subplots(figsize=(6, 5))
    
    # Calculate means and SEM
    groups = [group_a, group_b]
    means = [data_a.mean(), data_b.mean()]
    sems = [data_a.sem(), data_b.sem()]
    
    # Bar plot
    bars = ax.bar(groups, means, yerr=sems, capsize=5, 
                   color=['#3498db', '#e74c3c'], alpha=0.7, edgecolor='black')
    
    # Add individual points
    for i, (group_name, group_data) in enumerate([(group_a, data_a), (group_b, data_b)]):
        x = np.random.normal(i, 0.04, size=len(group_data))
        ax.scatter(x, group_data, alpha=0.6, s=50, color='black', zorder=3)
    
    # Add p-value annotation
    sig_text = f"p={pval:.4f}"
    if pval < 0.001:
        sig_text += " ***"
    elif pval < 0.01:
        sig_text += " **"
    elif pval < 0.05:
        sig_text += " *"
    else:
        sig_text += " ns"
    
    ax.text(0.95, 0.95, sig_text, transform=ax.transAxes, 
            ha='right', va='top', fontsize=12, 
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Labels
    metric_name = metric.replace('_', ' ').title()
    ax.set_ylabel(metric_name, fontsize=12)
    ax.set_title(f'{tissue}: {metric_name}', fontsize=14, fontweight='bold')
    ax.set_xlabel('')
    
    # Grid
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)
    
    plt.tight_layout()
    
    # Save
    png_path = OUTPUT_BASE / 'PNG' / f'{output_prefix}_{metric}_{tissue}.png'
    svg_path = OUTPUT_BASE / 'SVG' / f'{output_prefix}_{metric}_{tissue}.svg'
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(svg_path, bbox_inches='tight')
    plt.close()
    
    return {
        'metric': metric, 
        'tissue': tissue, 
        'p_value': pval, 
        f'{group_a}_mean': data_a.mean(), 
        f'{group_b}_mean': data_b.mean(),
        f'{group_a}_n': len(data_a), 
        f'{group_b}_n': len(data_b)
    }

# ============================================================================
# MAIN ANALYSIS - TEMPLATE
# ============================================================================

print("\n" + "="*80)
print("TEMPLATE: Implement your analysis below")
print("="*80)

# Example workflow:
# results = []
# for tissue, df in [('Compartment_1', df_compartment1), ('Compartment_2', df_compartment2)]:
#     for subset in SUBSETS:
#         for norm in NORMALIZATIONS:
#             metric = f'{subset}_{norm}'
#             if metric in df.columns:
#                 result = create_bar_plot_with_points(df, metric, tissue, GROUP_A, GROUP_B, 'comparison')
#                 if result:
#                     results.append(result)

# Save results
# df_results = pd.DataFrame(results)
# df_results.to_csv(OUTPUT_BASE / 'CSV' / 'statistical_results.csv', index=False)

print("\n✅ Template ready for customization")
print("="*80)
