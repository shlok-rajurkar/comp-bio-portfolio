# Curve deconvolution script for fityk

Multi-step process for taking 1200 pt ion mobility data and exporting into deconvoluted curve data. Tools are generally flexible and should be used with consideration of their effects.

## File breakdown

Fityk folder contains a simple .fit script for placing and fitting Gaussians on spectra as well as a more complex .lua script that places a peak in each given range and has an option to export file into a given folder. The ranges for peaks and range of interest can be configured quite easily based on the required situations.  

JMP folder contains a python script to process a correlation matrix for spectral analysis.  

R folder contains some assorted data cleaning functions pertaining to empty values and file storage.
