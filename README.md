# Bennu Asteroid Shape Model Visualization and STL Export

![Asteroid Bennu](https://www.nasa.gov/wp-content/uploads/2019/03/bennu-full-mosaic-center-crop-9-25-2018.jpg)

## Overview

This MATLAB project provides tools for analyzing, visualizing, and exporting 3D shape models of asteroid (101955) Bennu. It enables researchers and 3D printing enthusiasts to work with NASA's radar-derived shape model data in various formats.

## Features

- **Complete shape model visualization** with multiple viewing options
- **Geometric analysis** of the asteroid (volume, surface area, dimensions)
- **Export to STL format** for 3D printing and CAD applications
- **Multiple scaling options** for different uses
- **Comprehensive scientific metadata** and contextual information
- **Robust error handling** for reliable operation

## File Structure

- `bennu.m` - Main visualization and analysis script
- `export_stl.m` - Standalone STL export utility
- `data/` - Contains source data files in TAB and XML formats
  - `101955bennu.tab` - Primary shape model file (vertices and faces)
  - `pole.tab` - Pole orientation data
  - `rotate.tab` - Rotation period data
  - Additional XML documentation files
- `3d_model/` - Contains exported STL files
  - `bennu_original_scale.stl` - Model in original scale (kilometers)
  - `bennu_3d_print_10000.stl` - Model scaled for 3D printing (1:10,000)
  - `bennu_3d_print_1000.stl` - Larger 3D printing model (1:1,000)

## Usage

### Complete Analysis

Running the program will:
1. Load the asteroid shape model from TAB files
2. Calculate geometric properties
3. Create comprehensive visualizations
4. Export 3D models in STL format

### STL Export Only

This will create STL files at three different scales:
- Original scale (kilometers)
- 3D printing scale 1:10,000 (1 km = 100 mm)
- 3D printing scale 1:1,000 (1 km = 1000 mm)


## Data Source

The shape model data is derived from radar observations conducted by the Arecibo and Goldstone observatories in 1999 and 2005, as described in:

> Nolan, M.C., et al. (2013). "Shape model and surface properties of the OSIRIS-REx target Asteroid (101955) Bennu from radar and lightcurve observations." Icarus, 226(1), 629-640.

## Requirements

- MATLAB R2014b or newer
- No additional toolboxes required

## Acknowledgments

Shape model data courtesy of NASA Planetary Data System (PDS).
