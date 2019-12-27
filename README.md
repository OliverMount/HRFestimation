# HRFestimation
Estimation of resting state HRF and its parameters using Pseudo point process approach
The original MATLAB source codes from other authors can be found here 
https://github.com/brainlife/app-rsHRF

This page contains R codes written by us for the Impact of HRF deconvolution on rs-networks 

1. to compute the resting-state HRF and its parameters
   (The main R file to obtain the rs-HRF and its parameters is Give_HRF_para.R 
   The R functions required for the main R file are 
      a. EstiHRF.R
      b. Give_Events.R
      c. NeuralEvent.R
      d. find_peaks.R
      e. IBasis.R
      f. convolmtx.R
2. to compute the effective filter variances use Eff_Filt_Var.R
