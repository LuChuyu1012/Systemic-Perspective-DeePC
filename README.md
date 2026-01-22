[![MATLAB R2023a](https://img.shields.io/badge/MATLAB-R2023a-green)](https://nl.mathworks.com/products/new_products/release2023a.html)


# A Systemic Perspective on Input Design for Data-Driven Predictive Control
Reference MATLAB implementation accompanying "A Systemic Perspective on Input Design for Data-Driven Predictive Control", providing scripts for experiment generation, DeePC trajectory tracking, and result visualization.
## Prerequisites

You will need:

- `MATLAB` 
- `MATLAB Statistics Toolbox` (see [MathWorks product page](https://nl.mathworks.com/products/statistics.html))
- `CVX` (see [cvxr.com](https://cvxr.com/))
- `Gurobi` (see [gurobi.com](https://www.gurobi.com/))

## Installation

To clone this repository, see https://nl.mathworks.com/help/simulink/ug/clone-git-repository.html

## Overview

This repository implements a set of DeePC trajectory-tracking experiments on six discrete-time LTI systems.
For each system, DeePC is run with different input designs, and the resulting closed-loop tracking trajectories
are compared to highlight how the choice of input data affects tracking performance.

Performance is mainly evaluated using:
- Convergence time
- RMSE (root-mean-square tracking error)

The experiments are organized into two categories:
1. Noise-free setting (no measurement noise)
2. Noisy setting (with measurement noise)

All generated results, including figures and saved data files, can be found in the corresponding results folders.
