# GARPTools
GARPTools: software for data prep and model evaluation for GARP species distribution modeling

The GARPTools package provides tools to prepare data for input into the desktop version of GARP 
(the genetic algorithm for rule-set prediction) and for the evaluation of the accuracy of output models. 
DesktopGARP is a software package for biodiversity and ecological research that allows the user to predict 
and analyze wild species distributions. GARP is a presence-only genetic algorithm that models species' potential 
geographic distributions through an iterative process of training and testing that occurs through resampling 
and replacement of input data.

In particular, the package GARPTools was designed to support:
•	Resampling of environmental layers to the same spatial scale and extent
•	Cropping of environmental layers to study area
•	Removal of duplicate presence points within a single raster cell
•	Splitting presence data into training (for input into DesktopGARP) and testing data points
•	Summation of GARP model outputs into a single best subset raster
•	Calculation of the confusion matrix, including model sensitivity and specificity
•	Assessment of model accuracy with AUC, commission, and omission metrics
•	Plotting of the Receiver Operating Characteristic curve
•	Plotting the range rules of variables output by DesktopGARP

