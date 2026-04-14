# Eddy-Centric Along-Track Altimetry Analysis

MATLAB tools for extracting, modeling, and analyzing mesoscale ocean eddies from satellite along-track altimetry data. This framework fits parametric eddy models (Gaussian and elliptical) to along-track sea surface height (SSH) observations, enabling high-resolution characterization of eddy structure and evolution.

## Associated Publication

> C.-Y. Ohh, P. Gaube, J. Early, B. Curtis, and J. Lilly, "Model-guided approach to access high-resolution mesoscale eddy dynamics from along-track altimetry," *Geophysical Research Letters* (in preparation).

## Getting Started

The repository includes four demo live scripts that illustrate the main workflows:

| Script | Description |
|--------|-------------|
| `test_analyticalEddy.mlx` | Synthetic eddy observations and Observing System Simulation Experiments (OSSEs) |
| `test_realEddy.mlx` | Fitting models to real satellite altimetry data (e.g., Agulhas Ring) |
| `test_QGModelEddy.mlx` | Eddies from quasigeostrophic model output |
| `test_multiMissions.mlx` | Multi-mission satellite composite analysis |

## Repository Structure

```
alongtrackFitting/
├── Core functions
│   ├── alongtrackFromLatLonDomain.m      # Extract along-track data within a domain
│   ├── extractAlongtrackLatLonEddyCenter.m
│   ├── findEddyCentroid.m                # Eddy center detection from SSH
│   ├── analyticalEddyModel.m            # Parametric eddy model (Gaussian/elliptical)
│   ├── FitAlongTrackXYToEddyModel.m     # Fit model to along-track data
│   ├── FitAlongTrackXYToEddyModelWindowed.m  # Time-windowed fitting
│   ├── composite2D.m                    # 2D radial composites
│   └── compute_model_error.m            # Model error metrics
├── Visualization
│   ├── plotAlongtrack.m
│   ├── visualizeEddyFit.m
│   ├── plotFitWindowed.m
│   └── makePropagatingVideo.m
├── Utilities
│   ├── latlon2xy_centered.m
│   ├── radialProfile.m
│   └── generateTrackingBounds.m
└── Demo scripts (*.mlx)
```

## Contact

C.-Y. Ohh — cohh@uw.edu
Applied Physics Laboratory, University of Washington
