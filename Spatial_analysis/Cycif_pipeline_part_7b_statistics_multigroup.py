"""
CYCIF PIPELINE PART 7B: MULTI-GROUP STATISTICS
Compares multiple experimental groups using Kruskal-Wallis test.

NOTE: This is a template script. Users must:
- Define their own group names
- Provide their own data files
- Customize metrics as needed
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import kruskal
from pathlib import Path

BASE_DIR = Path('output/analysis')
OUTPUT_DIR = BASE_DIR / 'statistics_multigroup'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("="*80)
print("CYCIF PIPELINE PART 7B: MULTI-GROUP STATISTICS")
print("="*80)

# ============================================================================
# CONFIGURATION - MODIFY FOR YOUR STUDY
# ============================================================================

# Define your groups (pseudonymized)
GROUPS = ['GROUP_A', 'GROUP_B', 'GROUP_C']
GROUP_COLORS = {
    'GROUP_A': '#2ecc71',  # Green
    'GROUP_B': '#e74c3c',  # Red
    'GROUP_C': '#3498db'   # Blue
}

# Key metrics to analyze
KEY_METRICS = [
    'Treg_per_CD45',
    'CD4_T_per_CD45',
    'CD8_T_per_CD45',
    'T_Cells_per_CD45',
    'CD45_per_mm2',
    'Treg_per_CD4'
]

# Load data - MODIFY PATHS
# df_lp = pd.read_csv(BASE_DIR / 'LP' / 'lamina_propria_metrics.csv')
# df_crypt = pd.read_csv(BASE_DIR / 'Crypt_IEL' / 'crypt_level_metrics.csv')

print("\nNOTE: This is a template. Load your own data files.")
print("Required columns: 'sample', 'group', and metric columns")

def create_multigroup_boxplot(df, metric, tissue, groups, group_colors, output_dir):
    """
    Create boxplot comparing multiple groups with Kruskal-Wallis test.
    
    Parameters:
    -----------
    df : DataFrame
        Data containing 'group' column and metric values
    metric : str
        Column name of metric to plot
    tissue : str
        Tissue compartment name for title
    groups : list
        List of group names to compare
    group_colors : dict
        Dictionary mapping group names to colors
    output_dir : Path
        Output directory for plots
    
    Returns:
    --------
    dict : Result dictionary with statistics
    """
    
    if metric not in df.columns:
        return None
    
    # Get data for each group
    group_data = {}
    for group in groups:
        data = df[df['group'] == group][metric].dropna()
        if len(data) >= 2:
            group_data[group] = data
    
    if len(group_data) < 2:
        return None
    
    # Kruskal-Wallis test (non-parametric ANOVA)
    H, p = kruskal(*group_data.values())
    
    # Create plot
    fig, ax = plt.subplots(figsize=(8, 5))
    
    positions = list(range(1, len(group_data) + 1))
    data_list = [group_data[g] for g in groups if g in group_data]
    colors_list = [group_colors[g] for g in groups if g in group_data]
    labels_list = [g for g in groups if g in group_data]
    
    # Box plot
    bp = ax.boxplot(data_list, positions=positions, widths=0.5, patch_artist=True,
                    showfliers=False)
    for patch, color in zip(bp['boxes'], colors_list):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    # Add individual points
    for i, (data, pos, color) in enumerate(zip(data_list, positions, colors_list)):
        ax.scatter([pos]*len(data), data, color=color, alpha=0.6, s=50, zorder=3)
    
    # P-value annotation
    sig = '***' if p < 0.001 else '**' if p < 0.01 else '*' if p < 0.05 else 'ns'
    ax.text(np.mean(positions), ax.get_ylim()[1]*0.95, 
            f'Kruskal-Wallis p={p:.4f} {sig}', 
            ha='center', fontsize=10, 
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.set_ylabel(metric.replace('_', ' '))
    ax.set_xticks(positions)
    ax.set_xticklabels(labels_list)
    ax.set_title(f'{tissue}: {metric}')
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    
    # Save
    output_path = output_dir / f'{tissue}_{metric}.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    # Compile results
    result = {
        'tissue': tissue,
        'metric': metric,
        'p_value': p,
        'significant': p < 0.05
    }
    
    for group in groups:
        if group in group_data:
            result[f'{group}_n'] = len(group_data[group])
            result[f'{group}_mean'] = group_data[group].mean()
            result[f'{group}_std'] = group_data[group].std()
    
    return result

# ============================================================================
# MAIN ANALYSIS - TEMPLATE
# ============================================================================

print("\n" + "="*80)
print("TEMPLATE: Implement your analysis below")
print("="*80)

# Example workflow:
# results = []
# for tissue, df in [('Compartment_1', df_comp1), ('Compartment_2', df_comp2)]:
#     for metric in KEY_METRICS:
#         result = create_multigroup_boxplot(df, metric, tissue, GROUPS, GROUP_COLORS, OUTPUT_DIR)
#         if result:
#             results.append(result)
#             sig = "***" if result['p_value'] < 0.001 else "**" if result['p_value'] < 0.01 else "*" if result['p_value'] < 0.05 else "ns"
#             print(f"  ✅ {metric}: p={result['p_value']:.4f} {sig}")

# Save results
# df_results = pd.DataFrame(results)
# df_results.to_csv(OUTPUT_DIR / 'statistics_summary.csv', index=False)

print("\n✅ Template ready for customization")
print("="*80)
