# Epoxy_CyNif Pipeline v1.0

🔬 **Epoxy-based Cyclic Nanobody Immunofluorescence (Epoxy_CyNif)** – End-to-end workflow for multicycle staining/imaging: illumination correction, stitching/registration/deconvolution+EDF, spillover/autofluorescence removal, segmentation, GUI-based QC (CyLinter), and spatial statistics.

📖 **[View Documentation](https://kaloshi.github.io/Epo_CyNif/)**

## TL;DR
- **Inputs:** Multicycle TIFF/CZI stacks plus marker CSV per sample.
- **Paths:** `BASE_EXPORT = C:\Users\researcher\data\Epoxy_CyNif\data\export\<sample>`.
- **Workflow:** Part 1 → 2 → 3 → 4a/4b → 5 (CyLinter QC) → 6/7 (Statistics).
- **Day 2+:** In Part 5, run cells 2–3 + extended GUI (Cell 12), delete checkpoints selectively, restart from desired module.
- **Outputs:** Segmentation/Processed/Clustering in `BASE_EXPORT`, CyLinter checkpoints/reports in `cylinter_output_prune_test/`, statistics in `python/analysis-V18/`.

## Pipeline Modules

### 1) Part 1 – Illumination/Stitching/Registration/Decon&EDF
`Illumination_correction,Stiching,Registration,Decon&EDF/Cycif_pipeline_part_1_stiching_REG_decon_EDF.ipynb`  
Illumination correction, stitching, registration, deconvolution/EDF; generates fused multichannel stacks.

### 2) Part 2 – Spillover Removal
`Spillover&AF_removal/Cycif_pipeline_part_2_Spillover.ipynb`  
Mutual information spillover estimation, 0–1 normalization, back-transformation to original data type.

### 3) Part 3 – Autofluorescence Removal (ACE)
`Spillover&AF_removal/Cycif_pipeline_part_3_AF_REV_ACE.ipynb`  
Autofluorescence correction on the image stack.

### 4a) Segmentation
`Segmentation/Cycif_pipeline_part_4a_segmentation.ipynb`  
Set `SAMPLE_ID`, load marker CSV from `BASE_EXPORT`, DAPI optimization (best preprocessing variant), InstanSeg nuclei/cells, feature exports.

### 4b) Batch Segmentation Helper
`Segmentation/Cycif_pipeline_part_4b_batch_MICROSAM_ROBUST_V8.ipynb`  
Batch/robust segmentation for multiple samples.

### 5) CyLinter QC (GUI)
`Cylinter/Cycif_pipeline_part_5_cylinter.ipynb`  
- **Day 1:** Complete run with marker selection + pipeline start.  
- **Day 2+:** Extended GUI, delete checkpoints, restart from module X.  
Uses `cylinter_config.yml` + `markers.csv` in notebook folder; checkpoints/reports in `cylinter_output_prune_test/`.

### 6/7) Spatial Analysis & Statistics (v18)
- `Spatial_analysis/Cycif_pipeline_part_6_run_data_frames.py`  
- `Spatial_analysis/Cycif_pipeline_part_7_statistics_v18_2groups_CORRECT.py`  
- `Spatial_analysis/Cycif_pipeline_part_7_statistics_healthy_vs_mutations_v18_SIMPLE.py`  

Frequency metrics and group statistics for LP/Crypt and mutation groups; expects inputs in `python/analysis-V18/`.

## Data Layout
- **Sample root:** `BASE_EXPORT = C:\Users\researcher\data\Epoxy_CyNif\data\export\<sample>`
- **Markers:** `markers_<SAMPLE_ID>.csv` in sample folder.
- **Outputs:**  
  - Segmentation/Processed/Clustering → `BASE_EXPORT`  
  - CyLinter Checkpoints/Reports → `cylinter_output_prune_test/`  
  - Statistics → `python/analysis-V18/`

## Workflow
- **Day 1 (Initial run):** Part 1 → Part 2 → Part 3 → Part 4a/4b → Part 5 (CyLinter) → Part 6/7.
- **Day 2+ / Error correction:**  
  - Restart Part 5 kernel, run cells 2–3.  
  - Open extended GUI (Cell 12), delete checkpoint(s), select start module, run.  
  - Interactive GUIs: `selectROIs`, `setContrast`, `gating`. Markers without threshold will be re-gated.

## Requirements
```
python 3.9.x–3.11.x
pandas>=2.2.0
numpy>=1.26.0
scipy>=1.12.0
scikit-image>=0.22.0
tifffile>=2024.1.0
shapely>=2.0.0
matplotlib>=3.8.0
seaborn>=0.13.0
ipywidgets>=8.1.0
pyyaml>=6.0.0
instanseg  # Model: fluorescence_nuclei_and_cells
cylinter   # CLI/GUI
napari     # Optional: visual inspection
```

See [requirements.txt](requirements.txt) for full dependency list.

## Installation
```bash
pip install -r requirements.txt
```

## QC / Checkpoints (CyLinter)
- **Checkpoints:** `cylinter_output_prune_test/checkpoints/*.parquet`
- **Gating thresholds:** `cylinter_report.yml` (persists across runs)
- `--module X` sets start point; end point is always pipeline end; checkpoints skip completed modules.

## Key Features
- Reconstructs high-dimensional cyclic IF nanobody stacks (Epoxy_CyNif) across all cycles.
- Removes technical artifacts (illumination, spillover, autofluorescence) for reproducible quantification.
- Performs robust cell/nucleus segmentation with DAPI optimization and exports feature tables.
- Provides interactive QC (CyLinter) with checkpointing for fast day-2+ iterations.
- Computes spatial frequencies and group statistics for downstream analyses.

## Citation
If you use this pipeline, please cite the following tools:

- **InstanSeg:** Goldsborough, T. et al. (2024) *InstanSeg: an embedding-based instance segmentation algorithm*. [arXiv:2408.15954](https://doi.org/10.48550/arXiv.2408.15954)
- **CyLinter:** Laboratory of Systems Pharmacology, Harvard Medical School. [GitHub](https://github.com/labsyspharm/cylinter)

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments
Third-party packages (InstanSeg, CyLinter, scikit-image, tifffile, shapely, pandas, numpy, seaborn, ipywidgets, pyyaml, napari) are used under their respective licenses (MIT/BSD/Apache).
