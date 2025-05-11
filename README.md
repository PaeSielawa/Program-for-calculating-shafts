# Shaft Analysis Program

This MATLAB program performs comprehensive mechanical analysis of a shaft, including:
- Force reactions in XZ and XY planes
- Shear forces diagrams
- Bending moments calculations
- Torsional moment analysis
- Shaft diameter calculations
- Deflection analysis using Clebsch method
- Critical speed calculations
- Bearing selection
- Fatigue calculations

## Features

- Calculates and visualizes:
  - Shear forces in XZ and XY planes
  - Bending moments
  - Torsional moments
  - Reduced moments (using HMH hypothesis)
  - Shaft diameters
  - Deflection angles
- Performs bearing selection calculations
- Includes keyway calculations
- Calculates critical speeds
- Performs fatigue analysis

## Input Parameters

The main input parameters include:
- Shaft dimensions (a, b, c, d)
- Forces (F1y, F1z, F2z)
- Material properties (E, Re, Rm)
- Operating conditions (n - rpm)

## Outputs

- Multiple plots showing forces, moments and diameters
- Calculated minimum shaft diameters for different sections
- Bearing loads and selection guidelines
- Critical speed values
- Fatigue analysis results

## Requirements

- MATLAB (tested on version 2021 or newer)
- No additional toolboxes required

## Usage

1. Open the MATLAB environment
2. Run the script `JZPS_program.m`
3. Review the generated plots and console output

## License

This project was created by PS. All rights reserved.
