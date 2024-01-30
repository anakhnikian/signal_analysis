# README

This repository contains a collection of spectral analysis and connectivity measures for analyzing physical and biological time series. 

See documentation of individual functions for a full list of functionality and options

## Spectral analysis 
Functions for converting time series into frequency domain representations including spectral density matrices, which form the core component of connectivity estimation. This directory also contains functions for extracting single-channel power spectral densities from full spectral matrices, and for converting individual time series to psds on the fly.

## Connectivity analysis
Functions for extracting connectivity patterns among univariate and multivariate time series in the frequency domain. Currently implement functions are

-standard coherence

-lagged coherence (removes zero-lag effect)

-partial coherence (regresses out mediated connections)

## Demonstration

Run the program demo_model_analysis to return power spectral densities, coherence, and lagged coherence for a simulated ten node network.  This program will generate the figures provided in the images/coherence_ar2_60cyc directory, and provides examples of function calling and variable passing to demonstrate pipeline development using the provided functions. 

### Demo Data: The attached .mat file contains the simulated data in a times by trials by channels array. 10 channels each consisting of 100 2 sec epochs were generated. The simulated time series consist of a set of coupled AR(2) processes update parameters 0.55 and -0.81 and a simulated sampling rate of 200 Hz, yielding a peak frequency of 40 Hz and temporal decay multiplier equal to 0.9, yielding a narrowband model.  The following channels are coupled with a time lag of 1 data index (0.005 sec) with coupling constant 0.5

{1 2} {1 6} {5 9} {7 10}

Zero lag interactions among all channels were introduced by adding identical 60 Hz sine noise to all channels, simulating voltage measurements corrupted with 60 cycle power line noise. Note that this component of the variance is present in power spectra and standard coherence data, but is removed from the lagged coherence spectra.

## Dependencies

Data tapers require the Matlab signal processing toolbox. Internal implementations of each taper are on my to-do list.
