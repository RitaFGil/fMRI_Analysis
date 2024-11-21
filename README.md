# fMRI_Analysis
FMRI Analysis with different stimulation conditions

This code calculates the averaged cycle (mean +/- sem) for a group of animals that underwent fMRI acquisitions with different conditions. 

Initial parameters include condition folder names and numbers, a reference animal to which the group was normalized to (if no normalization was performed just make sure you input a number higher than your total N), stimulation vector and repetition time.

Nifties are loaded and converted into matrices, regions of interest can be delineated per acuiqisiton, detrending is done with a polinomial fit to resting periods, and individual cycles are averaged per condition in the end.
____

For each run the "steady state" of the response is identified and calculated per condition. This involdes resampling of the datasets.

_____

Based on steady state quantifications, a neurometric curve is ploted for the different stimulation conditions.
