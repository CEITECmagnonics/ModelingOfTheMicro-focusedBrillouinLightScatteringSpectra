# ModelingOfTheMicro-focusedBrillouinLightScatteringSpectra
This repository have codes used for paper titled Modeling of the micro-focused Brillouin light scattering spectra. The scripts can be used to calculate the signal of thermal micro and BLS and of BLS on coherently excited spin waves.

# Master codes
These codes generates the presented results.

## BLSPtSweep.m
Code for calculating the thermal signal with Pt cover layer

## BLSsignalThicknessSweep.m
Code for calculating thermal spectra of the magnetic layers with different thicknesses

## FieldSweep.m
Code for calculating thermal spectra of the magnetic layers at different magnetic fields

## FillingFactorSWeep.m
Code for calculating thermal spectra of the magnetic layers with different filling factor of objective lens

## NASweepCoherent.m
Code for calculating detection sensitivity of uBLS to cohrently excited spin waves detected by objective lenses with different NA.

# Internal functions
These codes are used as an internal functions.

## focalFields.m
Code for calculating the electric fields from given objective lens. This function is called internally by the script. The incident polarization on the objective lens is assumed linear.

## focalFieldsAzm.m
Code for calculating the electric fields from given objective lens. This function is called internally by the script. The incident polarization on the objective lens is assumed azimuthal.

## focalFieldsRadial.m
Code for calculating the electric fields from given objective lens. This function is called internally by the script. The incident polarization on the objective lens is assumed radial.

## Fresnelq.m
Code for calculating Fresnel coefficients.

## GetBLSsignal.m
Code for calculating the thermal BLS signal.

## GetBLSsignalCoherent.m
Code for calculating the coherent BLS signal.

## GetBLSsignalCoverLayer.m
Code for calculating the coherent BLS signal with presence of platinum cover layer.

## GetBLSsignalCoverLayerSiO2.m
Code for calculating the coherent BLS signal with presence of SiO2 cover layer.

## MatrixMultiplier.m
Code for multipluication of the matrices

## multiWaitbar.m
Code for showing the progress of the calculation. Available from: https://www.mathworks.com/matlabcentral/fileexchange/26589-multiwaitbar. Author: Ben Tordoff

## rowwiseLast.m
Finds locations k(i) of first non-zero value in each row A(i,:), with k(i)=NaN if the entire row is zero.

## sphGreenFunction.m
Gives the spherical Green function

## SpinWaveGreen.m
Gives the density of states (imaginary part of the Green function) of spin waves in thin films using aproximation by Lorentzian oscillators.

## SpinWaveToolkit.py
Python function to give spin wave eigenstates. Here it is depricated in comparison to current development: https://github.com/CEITECmagnonics/SpinWaveToolkit

## SubstrateFun.m
Function to give characteristics of various substrates.


