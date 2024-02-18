# Particle Dynamics Simulation

This Python script simulates the dynamics of particles in a 3D space considering electrostatic and gravitational forces between the particles.

## Table of Contents
- [Introduction](#introduction)
- [Requirements](#requirements)
- [Usage](#usage)
- [Input](#input)
- [Simulation Details](#simulation-details)
- [Graphs](#graphs)
- [External Functions](#external-functions)
- [Time Measurement](#time-measurement)
- [License](#license)

## Authors

- [@self-Puneet](https://github.com/self-Puneet)
- [@self-nasu](https://github.com/self-nasu)


## Introduction

The script performs a particle simulation over a specified time duration, calculating forces, velocities, displacements, and various physical quantities at each iteration. The simulation results are visualized through graphs using Matplotlib.

## Requirements

- Python 3.x
- NumPy
- Matplotlib

## Usage

1. Clone the repository:

   ```bash
   git clone https://github.com/Self-nasu/particle-simulator
   ```

2. Navigate to the project directory:

   ```bash
   cd your-repository
   ```

3. Run the simulation script:

   ```bash
   python particle_simulation.py
   ```

## Input

- The script prompts the user to input the number of particles, simulation time, and the number of frames.
- Initial dimensions and velocities of each particle are also inputted by the user.

## Simulation Details

- The simulation loop iterates through time steps, calculating forces, velocities, displacements, and updating particle positions.
- Various physical quantities, such as speed, linear momentum, kinetic energy, and potential energies, are calculated at each iteration.

## Graphs

- The script generates multiple graphs using Matplotlib to visualize the simulation results:
  - Particle coordinates
  - Force plots
  - Speed plots
  - Linear momentum plots
  - Kinetic energy plots
  - Total potential energy plots
  - Overall energy plots

## External Functions

- The script calls external functions (`hello`, `coordinates_plot`, `force_plot`, etc.) to handle the plotting of graphs. Ensure these functions are correctly implemented.

## Time Measurement

- The script measures and prints the total time taken for the simulation and the time taken for one iteration.

## License

This project is licensed under the [GNU GENERAL PUBLIC LICENSE](LICENSE).