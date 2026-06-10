# Catchment Variability and Dynamic Signature Analysis

## Overview

This repository contains the code used to generate:

* Dynamic signatures using Wavelet Transform Coherence (WTC) analysis in MATLAB.
* Variability signatures using Python-based analysis.

The workflow reproduces the analyses performed in the study and can be adapted for datasets from other countries with minor modifications.

---

## Repository Structure

```text
.
├── Germany_WTC_Main_code.m
├── Germany_catchment_variability_signatures_code.py
└── README.md
```

---

## Prerequisites

### MATLAB

* MATLAB (tested with recent versions)
* Grinsted Wavelet Coherence Toolbox

### Python

* Python 3.8 or later

Required packages:

* Keras
* Pandas
* NumPy

Install dependencies using:

```bash
pip install keras pandas numpy
```

---

## Step 1: Install the Wavelet Coherence Toolbox

Download and install the Grinsted Wavelet Coherence Toolbox:

https://github.com/grinsted/wavelet-coherence

Add the toolbox directory to your MATLAB path before running the MATLAB scripts.

---

## Step 2: Prepare the Dataset

Download the provided dataset "https://doi.org/10.5281/zenodo.20630407" and extract the contents of `data.zip` into your working directory.

Ensure that all input files are available before executing the workflow.

---

## Step 3: Generate Dynamic Signatures

Run the MATLAB script:

```matlab
Germany_WTC_Main_code.m
```

This script performs Wavelet Transform Coherence (WTC) analysis and extracts dynamic signatures for countinuous scales.

The following temporal scale ranges are analyzed:

* 2-7 days
* 7-14 days
* 14-30 days
* 30-183 days
* 183-365 days

For each scale range, the lag and RMI signatures were calculated:

* RMI ratio
* Lag ratio

These metrics constitute the dynamic signatures used in the analysis.

---

## Step 4: Generate Variability Signatures

Run the Python script:

```bash
python Germany_catchment_variability_signatures_code.py
```

This workflow computes variability signatures from the hydrological input data.

---

## Step 5: Applying the Workflow to Other Countries

The current implementation is configured for the Germany dataset.

To apply the methodology to another country:

1. Replace the input datasets.
2. Update CSV file names.
3. Modify column names if the data structure differs.
4. Adjust input and output file paths.
5. Update output file names if required.

The underlying methodology and signature calculations remain unchanged.

---

## Notes

* Execute the MATLAB workflow before running the Python workflow.
* Ensure all required data files are present before execution.
* MATLAB and Python components are independent but follow a sequential workflow.
* Dynamic signatures are generated first using WTC analysis, followed by variability signature computation.

---

## Citation

If you use this repository in your research, please cite the associated publication.
