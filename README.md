# fNIRS Analysis Pipelines (MATLAB) 

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language](https://img.shields.io/badge/Matlab-R202X-orange.svg)](https://www.mathworks.com/products/matlab.html)

Welcome to this repository!

This repository contains the source code and analysis pipelines used in our research article. It provides the **10 processing pipelines** compared in our study, allowing for reproduction of our results on Time Series (Average) and GLM (Beta values) analyses.

### Reference Article

If you use this code or our findings, please cite our paper:

> **Optimizing short-channel regression in fNIRS: an empirical evaluation with ecological audiovisual stimuli**
>
> *Y. Lemaire, P. Barone, K. Strelnikov, O. Deguine, A. Caclin, J. Ginzburg.*
> Published in Neurophotonics

---

## Prerequisites & Requirements

To run these scripts successfully, you **must** have the following toolboxes installed in your MATLAB environment:

1.  **NIRS Brain AnalyzIR Toolbox** (Santosa et al.)
    * *Core requirement for GLM analyses.*
    * 📥 [Download here (GitHub)](https://github.com/huppertt/nirs-toolbox)
2.  **Homer2**
    * *Mandatory specifically for the **"Average"** (Time Series) analysis part.*
    * 📥 [Download here](https://homer-fnirs.org/)

---

## Repository Structure

The project is organized into two main directories corresponding to the two types of analyses performed:

### 1. `Average/` (Time Series Analysis)
Contains scripts for the 10 pipelines designed to analyze the temporal shape of fNIRS curves.

* **Goal:** Visualize the reconstructed hemodynamic response based on the pipeline used.
* **Main Script:** `Final_average_analysis_XXX.m` (Run this to execute analysis).
* **⚠️ Important Note:** The file `add_SC_reggressors_custom_function.m` contains custom helper functions. It **must** be placed in the **same directory** as the main script for the code to run.

### 2. `GLM/` (Beta Values)
Contains scripts for the 10 pipelines using the General Linear Model (GLM).

* **Goal:** Obtain "Beta" values (regression weights) based on the Hemodynamic Response Function (HRF).
* **Main Script:** `Main_GLM_XXX.m`.
* **Custom Modules:** The subfolder `GLM_custom_function/` contains modified versions of the native algorithms from the NIRS Brain AnalyzIR Toolbox.

---

## Installation: Specific Steps for GLM

🚨 **ATTENTION: Manual Configuration Required** 🚨

To ensure the GLM pipelines function correctly, you must manually integrate our custom modules into your local NIRS Toolbox installation.

The folder `GLM/GLM_custom_function` contains modified `.m` files (starting with `GLM_...` and `AddShortSeperationRegressors_...`).

**Procedure:**

1.  Locate the `GLM/GLM_custom_function` folder in this repository.
2.  **Copy all files** inside it (e.g., `GLM_custom.m`, `AddShortSeperationRegressors_custom.m`, etc.).
3.  **Paste them** into the source folder of your **NIRS Brain AnalyzIR Toolbox** installation at the following path:

```matlab
% Example Path
...\your_path_to\nirs-toolbox-master\+nirs\+modules\
```
    Without this step, Main_GLM.m will fail to find the custom analysis modules.

---

**Toolbox Citation** 

This work relies heavily on the NIRS Brain AnalyzIR Toolbox. Please also cite the original work of the toolbox developers:

    The NIRS brain AnalyzIR toolbox. Santosa, H., Zhai, X., Fishburn, F., & Huppert, T. (2018). Algorithms, 11(5), 73. https://doi.org/10.3390/a11050073

<p align="center"> <sub>This code is provided "as is" to support research reproducibility. Feel free to open an Issue if you encounter difficulties.</sub> </p>
