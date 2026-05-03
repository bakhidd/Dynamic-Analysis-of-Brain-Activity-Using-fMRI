# Create README.md file with proper markdown content

content = """# Dynamic Analysis of Brain Activity Using fMRI

This repository contains an implementation of a full pipeline for analysis and prediction of brain activity from fMRI data using nonlinear manifold learning and dynamical systems methods.

The project was developed as part of a Bachelor Thesis in Biomechanics and Medical Engineering, with a focus on brain signal analysis and data-driven modeling.

---

## Project Overview

The goal of this project is to:

- Extract low-dimensional representations of brain activity  
- Model temporal dynamics of fMRI signals  
- Predict future brain states under external stimuli  
- Reconstruct predicted signals back into the original fMRI space  

The approach combines machine learning, spectral methods, and dynamical systems theory.

---

## Pipeline Structure

### 1. fMRI Preprocessing

- Loading fMRI scans in NIfTI format  
- Brain parcellation using the AAL atlas  
- Extraction of ROI time series  
- Detrending, filtering, and normalization  

**Implemented in:**
testing.py

---

### 2. Manifold Learning (Diffusion Maps)

- Construction of similarity kernel  
- Markov normalization  
- Eigenvalue decomposition  
- Selection of parsimonious eigenvectors  

**Implemented in:**
reproduce_diffusion_maps_fmri.py

---

### 3. Reduced Order Modeling

#### Neural Networks (FNN)

- Separate feedforward network for each latent coordinate  
- Hyperparameter tuning using cross-validation  
- Iterative prediction of system dynamics  

#### Koopman Operator (EDMD)

- Linear approximation of nonlinear dynamics  
- Incorporation of external stimuli as control input  

#### Baseline

- Naive Random Walk (NRW)

---

### 4. Reconstruction to fMRI Space

- Geometric Harmonics for neural network predictions  
- Koopman modes for EDMD predictions  

---

### 5. Evaluation

The following metrics are used:

- Root Mean Squared Error (RMSE)  
- L2 norm  

Evaluation is performed on selected brain regions, primarily in the visual cortex.

---

## Output

The pipeline generates:

- Eigenspectrum plots  
- Latent space dynamics predictions  
- Reconstructed brain signals  
- Quantitative error comparisons  

Results are saved in:

/results

---

## Dataset

The project uses an fMRI Attention dataset with corresponding experimental conditions:

- attention  
- non-attention  
- static  
- fixation  

Stimuli are loaded from a MATLAB file:

conditions.mat
"""

