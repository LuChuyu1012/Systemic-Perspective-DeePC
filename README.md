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

## Benchmark case study: LTI reference tracking with measurement noise

In the experiments, a total of six LTI systems are considered. This section presents one representative system as a benchmark case study.
```math
\begin{equation*}
x(t+1)=Ax(t)+Bu(t),\qquad y(t)=Cx(t)+v(t).
\end{equation*}
```


Here, $x(t)\in\mathbb{R}^2$ is the state, $u(t)\in\mathbb{R}$ is the input, and $y(t)\in\mathbb{R}$ is the measured output.  
The term $v(t)$ denotes the measurement noise and is modeled as Gaussian white noise:
```math
\begin{equation*}
v(t)\sim\mathcal{N}(0,10^{-4}).
\end{equation*}
```


The system matrices used in this benchmark are:
```math
\begin{equation*}
A=
\begin{bmatrix}
0.7326 & -0.0861\\
0.1722 & 0.9909
\end{bmatrix},\qquad
B=
\begin{bmatrix}
0.0609\\
0.0064
\end{bmatrix},\qquad
C=
\begin{bmatrix}
0 & 1.4142
\end{bmatrix}.
\end{equation*}
```

To ensure a fair comparison across input-design cases, the measurement noise realization is kept identical in all experiments, and the input amplitude is adjusted such that the signal-to-noise ratio (SNR) is the same across cases. The measured output used to construct the Hankel matrix is

```math
y(t)=y_0(t)+v(t),
```
where $y_0(t)$ denotes the noise-free output. The SNR is defined as

```math
\begin{equation*}
\mathrm{SNR}_{\mathrm{dB}}
=10\log_{10}\!\left(\frac{\sum_{t=1}^{T} \left|y_0(t)\right|^2}{\sum_{t=1}^{T} \left|v(t)\right|^2}\right).
\end{equation*}
```

### Input-design cases

To study how the excitation content affects DeePC performance, four input-design cases are considered:

- **Case 0 (white noise):** u(k) is generated as Gaussian white noise.
- **Case 1–3 (multisine):** u(k) is generated as a periodic multisine signal with identical phase construction across cases, but with different excitation bands.

The spectra of the four inputs are reported below to highlight the excited frequency components.


### Reference tracking visualization

The following GIF showcases the reference tracking performance achieved by the DeePC formulation under the above four input-design cases.



