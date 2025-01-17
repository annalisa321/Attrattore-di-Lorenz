# Lorenz System Simulation & Ensemble Analysis

This repository contains two Python scripts for simulating the Lorenz system and conducting an ensemble analysis with perturbations, respectively. Both scripts numerically solve the Lorenz system using Euler's method and analyze the sensitivity of the system to small changes in initial conditions. Visualizations and error analysis (RMSE) are also included.

## Files

- **Lor.py**: Simulates the Lorenz system and analyzes the Root Mean Square Error (RMSE) for different perturbations.
- **Ensemble.py**: Performs an ensemble analysis on multiple simulations, comparing the results of different runs and analyzing the statistical properties of the system.

---

## 1. **Lor.py: Lorenz System Simulation and Perturbation Analysis**

### Overview

The `Lor.py` script simulates the Lorenz system with different initial conditions, applies small perturbations, and computes the Root Mean Square Error (RMSE) for each case. Several plots are generated to visualize the differences and errors.

### Key Steps

1. **Imports and Setup**: The script imports necessary libraries such as `numpy`, `matplotlib`, and `pandas`, and sets up the initial conditions and parameters for the Lorenz system simulation.

2. **Euler’s Method**: The script uses Euler’s method to numerically solve the differential equations of the Lorenz system for both the main system and the perturbed versions. The system is updated at each timestep.

3. **Perturbed Systems**: Perturbed versions of the system are computed with slight changes in initial conditions. The system is then simulated for each perturbation.

4. **RMSE Calculation**: The Root Mean Square Error (RMSE) is calculated between the reference system and the perturbed systems at each timestep.

5. **Visualization**: The script generates several plots:
   - **Lorenz Attractor Plot**: 2D plot of \(X\) vs \(Z\)
   - **3D Lorenz Plot**: 3D plot of \(X\), \(Y\), and \(Z\)
   - **Difference Between \(X(t)\)**: Plots the difference between the reference and perturbed \(X(t)\)
   - **RMSE Plots**: RMSE comparison for different perturbations
   - **Final Table**: Displays the times when RMSE exceeds 0.5 for each perturbation
