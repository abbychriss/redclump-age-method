This code base contains all the relevant python code used in the paper "On the Age Calibration of Open Clusters using Red Clump Stars" AJ 169:81 (2025) (https://doi.org/10.3847/1538-3881/ada1c4) written by Abby Chriss and Guy Worthey.

# WSU Open Cluster Age Determination Pipeline

## Overview
This codebase implements a complete pipeline for determining open cluster ages using multi-wavelength photometry (Gaia DR3, PanSTARRS) and isochrone fitting. The primary method uses **Red Clump Age Dating**—comparing color differences between red giant branches and red clumps across multiple stellar evolution models.

---

## File Organization

### **1. Data Acquisition & Cone Searches**
Scripts for retrieving cluster data from astronomical archives.

| File | Purpose |
|------|---------|
| `gaia_cone_search.py` | Generic Gaia DR3 cone search template |
| `RC_cone_searches.py` | Batch orchestration for 25+ clusters |

### **2. Cluster Membership Determination**
Statistical analysis using proper motions, parallax, and position.

| File | Purpose |
|------|---------|
| `RC_cluster_member_study.py` | Multi-cluster red clump analysis (25+ clusters) |
| `RC_cluster_member_study2.py` | Enhanced version with PanSTARRS cross-matching |

### **3. Cluster Filtering & Data Management**
Filter and process cluster information from catalogues.

| File | Functions | Purpose |
|------|-----------|---------|
| `Cantat-Gaudin cluster ages.py` | N/A | Filter CG clusters by age (8.3-9.3 log years) and distance (<4 kpc) |
| `cluster_filters_CG.py` | N/A | Apply parallax and age constraints to Cantat-Gaudin data |
| `cluster_filters_CG_2.py` | N/A | Advanced CG filtering with uncertainty handling |
| `star_filters.py` | N/A | Filter CG member stars by cluster parameters |
| `cross_ref_csv.py` | N/A | Cross-reference Gaia and PanSTARRS catalogues |

### **4. PanSTARRS Data Interface**
Scripts for querying and cross-matching with PanSTARRS photometry.

| File | Purpose |
|------|---------|
| `panstarrs_query.py` | HTTP-based PanSTARRS DR2 catalogue API interface |
| `panstarrs_cmd.py` | Generate curl commands for cross-matching |
| `panstarrs_crossmatch.py` | Create batch PanSTARRS query scripts |
| `panstarrs_photometry.py` | Incomplete PanSTARRS REST API implementation |

### **5. Color-Magnitude Diagram & Visualization**
Create publication-quality plots of cluster populations.

| File | Functions | Purpose |
|------|-----------|---------|
| `plot_cmd.py` | N/A | Generic CMD plotter with extinction/distance corrections |
| `plot_cmd_panstarrs.py` | N/A | PanSTARRS-only CMD (g-i vs g colors) |
| `plot_cmd_pan+gaia.py` | N/A | Combined Gaia+PanSTARRS CMD comparison |
| `plotgaia.py` | `fetchiso(iiso)`, `fetchbasti(iiso)` | Gaia HR diagram with isochrone overlays |
| `isochrone_tracker.py` | N/A | Interactive point-clicker for manual isochrone extraction |

### **6. Red Clump Analysis**
Scripts focused on red clump identification and properties.

| File | Functions | Purpose |
|------|-----------|---------|
| `plot_RCs.py` | `fetchiso(iiso)`, `check_cluster(name)` | Red clump identification and median color extraction |
| `plot_RCs_2.py` | `fetchiso(iiso)`, `check_cluster(name)` | Advanced red clump analysis with weighting |
| `RC_luminosity_pileup.py` | `fetchiso(iiso)`, `fetchbasti(iiso)` | Analyze red clump luminosity distribution |
| `plot_smc_cmd.py` | N/A | SMC cluster red clump analysis (B-I photometry) |

### **7. Red Giant Branch (RGB) Color Analysis**
Extract RGB colors for age determination.

| File | Functions | Purpose |
|------|-----------|---------|
| `plot_rgb_median_color.py` | `f(x1,y1,y)` | Extract RGB color at red clump magnitude (linear interpolation) |
| `plot rgb median color.py` | `f(x1,y1,y)` | Alternate RGB median calculation |
| `rgb_error.py` | `fetchiso(iiso)`, `fetchbasti(iiso)` | Quantify photometric measurement uncertainties in RGB |
| `rgb_error_fit.py` | N/A | Fit error functional forms for error propagation |

### **8. Color Transformation & Theoretical Comparisons**
Calculate color indices across isochrone models.

| File | Functions | Purpose |
|------|-----------|---------|
| `d_BPRP.py` | `fetchiso(iiso)`, `fetchbasti(iiso)`, `f(x1,y1,y)` | BP-RP color differences (red giant vs red clump) |
| `d_BPRP_2.py` | `fetchiso(iiso)`, `fetchbasti(iiso)`, `fetchparsec(iiso)` | Multi-model BP-RP comparison (Padova, BASTI, PARSEC) |
| `d_BR.py` | `func(x, a)`, `func(x, a, b)` | Linear regression: BP-RP to B-R color transformation |
| `iso_d_br.py` | `fetchiso(iiso)`, `fetchisomp(iiso)`, `fetchbasti(iiso)`, `fetchbasti6(iiso)`, `fetchparsec(iiso)`, `fetchparsecmp(iiso)` | B-R color differences across 6 isochrone models |
| `plot_gaia_isochrones.py` | `fetchiso(iiso)`, `fetchbasti(iiso)`, `fetchbasti6(iiso)`, `fetchparsec(iiso)` | Overlay Gaia-band isochrones on CMDs |

### **9. Isochrone Data Management**
Load and parse stellar evolution model data.

| File | Functions | Purpose |
|------|-----------|---------|
| `isochrones.py` | `fetchiso(iiso)` | Load isochrone files; parse block-structured Hess diagrams |

### **10. Statistical Analysis & Fitting**
Linear regression and probability analysis.

| File | Functions | Purpose |
|------|-----------|---------|
| `linear_reg_ngc7789.py` | `func(x, a, b, c)`, `func(x, a, b, c, d, e, f, g)` | Polynomial/linear fitting for age relationships |
| `prob_linear_reg.py` | N/A | Probabilistic membership + regression analysis across clusters |

---

## Key Functions Reference

### Isochrone Loading Functions
Used across multiple analysis scripts to retrieve stellar parameters and magnitudes:

- **`fetchiso(iiso)`** — Load Padova isochrone models (Marigo 2007)
- **`fetchisomp(iiso)`** — Load Padova metal-poor isochrone models
- **`fetchbasti(iiso)`** — Load BASTI stellar evolution models
- **`fetchbasti6(iiso)`** — Load BASTI v6 models
- **`fetchparsec(iiso)`** — Load PARSEC isochrone models
- **`fetchparsecmp(iiso)`** — Load PARSEC metal-poor models

### RGB Color Extraction
- **`f(x1,y1,y)`** — Linear interpolation to extract RGB color at specified magnitude

### Utility Functions
- **`check_cluster(name)`** — Verify cluster name in reference tables
- **`func(x, a, ...)`** — Various fitting functions (linear, polynomial)

---

## Data Pipeline Flow

```
1. Cone Searches (Gaia/PanSTARRS) 
   ↓
2. Cluster Filtering (Age, distance, parallax constraints)
   ↓
3. Membership Determination (Proper motion analysis)
   ↓
4. Cross-Catalogue Matching (Gaia + PanSTARRS)
   ↓
5. Color-Magnitude Diagram Creation
   ↓
6. Red Clump & RGB Identification
   ↓
7. Isochrone Fitting & Color Comparison
   ↓
8. Age Determination (Red Clump Method)
```

---

## Scientific Method

**Red Clump Age Dating:** The code measures the color difference between red clump stars and the red giant branch at the clump luminosity to test whether there is a correlation between the delta(color) and clump age.

Multiple stellar evolution models (Padova, BASTI, PARSEC) are compared to validate results and assess systematic uncertainties.

---

## Configuration Files & Data

- **`cluster_parameters_Kharchenko_2013.txt`** — Reference cluster parameters
- **`metallicity_data.txt`** — Cluster metallicity values
- **`lightcone.txt`** — Physics diagram data
- **`rup37_manual.txt`** — Ruprecht 37 raw analysis data
- **`Isochrones/`** folder — Pre-computed stellar evolution model data
- **`Red Clump Age Method/`** folders — LaTeX articles and analysis notebooks

