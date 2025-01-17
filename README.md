# Lorenz System Simulation

This Python project simulates the behavior of the Lorenz system using Euler's method and compares the results under various initial conditions. The Lorenz system is a set of three nonlinear ordinary differential equations often used to model chaotic systems.

## 1. Overview of the Code

The code solves the Lorenz system with the following equations:

\[
\frac{dX}{dt} = \sigma(Y - X)
\]

\[
\frac{dY}{dt} = rX - XY - Y
\]

\[
\frac{dZ}{dt} = XY - bZ
\]

Where:
- \(X\), \(Y\), and \(Z\) are the state variables.
- \(\sigma\), \(r\), and \(b\) are parameters of the system.

The code is structured into several key sections:
- **Initialization**: Initializes parameters, arrays, and initial values.
- **Euler's Method**: Numerically integrates the Lorenz system using Euler's method.
- **Noise Implementation**: Introduces small perturbations to the initial conditions for multiple simulations (100 instances).
- **Error Calculation**: Computes the Root Mean Square Error (RMSE) for field averages and mean RMSE across all simulations.
- **Plotting**: Visualizes the RMSE and compares the accuracy of the numerical solution.

## 2. Parameters

The main parameters for the Lorenz system are:

- **sigma** = 10
- **r** = 28
- **b** = 8/3
- **Dt** (time step) = 0.005
- **t_start** = 0
- **t_end** = 4

Additionally, a random perturbation (\(\epsilon\)) is added to the initial condition for each of the 100 simulations.

## 3. Key Variables

- `X_arr`, `Y_arr`, `Z_arr`: Arrays storing the values of \(X\), \(Y\), and \(Z\) at each time step.
- `R`: Stores the Euclidean distance between the current state and the reference state (initial conditions).
- `S`, `XP`, `YP`, `ZP`: Mean values across all simulations.
- `RMSE`: Root Mean Square Error at each time step.

## 4. Euler’s Method

Euler's method is used to compute the system's evolution at each time step:
- The system evolves with a time step of \(Dt = 0.005\) from \(t_start = 0\) to \(t_end = 4\) seconds.

## 5. Plotting

The final output is a plot that shows:

- **Green Line**: The mean RMSE across all simulations.
- **Blue Line**: The RMSE of the mean field (average of \(X\), \(Y\), and \(Z\)).
- **Red Line**: RMSE computed with a third case (S_4), which evaluates the variability across the simulations.

The plot provides insight into the convergence and error behavior as the system evolves over time.
