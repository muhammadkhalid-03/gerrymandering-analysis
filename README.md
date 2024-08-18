# Gerrymandering Analysis using Metropolis-Hastings

## Project Overview

This project implements a version of the Metropolis-Hastings algorithm to perform simplified Gerrymandering analysis. It generates thousands of example districting plans on a square grid, which act as an ensemble for analysis.

## Files

- `hw7.py`: Main file containing the implementation of Metropolis-Hastings and Single Node Flip algorithm
- `gerry.py`: Helper functions for checks and plotting (provided by professor)
- `energy.py`: Implementation of energy functions for the sampling process (provided by professor)
- `.png`: Files containing plots before and after 10000 iterations

## Key Components

### Districting Plan Representation

- N x N matrix of 0s and 1s
- List of coordinate pairs (r, c) representing the border between districts

### Helper Functions

1. `isContiguous`: Checks if a given district is contiguous
2. `onBorder`: Checks if a specific precinct is on the border between districts
3. `countBorder`: Counts the number of 1s and 0s in the current border
4. `plotDistrict`: Plots the districting plan

### Main Functions to Implement

1. `updateBorder`: Updates the border after a node flip
2. `calcQ`: Calculates the transition probability Q(X → Y)
3. `metropolis`: Implements the Metropolis-Hastings algorithm

### Energy Functions

- Horizontal Symmetry
- Vertical Symmetry
- Horizontal Asymmetry
- Vertical Asymmetry
- Compactness (Isoperimetric ratio)

## Test Scenarios

- Uniform sampling: `filename = "uniform.png"`
- Left asymmetry: `filename = "left-asym.png"`
- Compactness: `filename = "compact.png"`
- Top-left priority: `filename = "top left.png"`

## Parameters

- N: Size of the N x N grid
- maxIters: Maximum number of iterations
- maxTime: Maximum computation time (in seconds)
- plotIters: Number of iterations between plots
- offset_thresh: Allowable offset between district sizes

## Output

- Generates plots of districting plans at specified intervals
- Saves images of the plans in the specified directory
