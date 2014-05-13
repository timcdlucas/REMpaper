REM paper
========

The analytical work (mostly) for REM paper.

## Background

#### Motivation

+ Camera traps and acoustic sensors are becoming cheap and widely used
+ Without individual recognition mark-recapture methods cannot be used
+ For widely moving animals, transect methods or occupancy models are unsuited
+ Currently Random Encounter Models do not account for complexities in detection

#### REM models

Ideal gas models assume that animals move randomly through space. If we examine our sensors motion relative to the animal population, we can work out the area covered by the sensor in a given period of time. This is the 'width' of the sensors detection region multiplied by average animal speed. If however the sensor has different widths depending on the angle of approach (i.e. a segment shape detection area of a camera trap) we must integrate the width for all approach angles, and then find the average width.

The case for detectors with detection angles &#8804;&#960;/2 has been examined by [Rowcliffe 2008](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2008.01473.x/abstract). I largely follow the notation within.


## The code

The main script is REM Analysis.py. This script solves the complicated REM models and uses a variety of tests to check that they are correct. It uses the symbolic algebra package sympy in python.

regionsPlot.R is a script for plotting where the different models are located in parameter space.

#### Models

This [pdf](imgs/regions.pdf) shows the qualitatively different model regions. Despite independant derivations, many of the models end up with the same solution. This is shown [here](imgs/equalRegions.pdf)

A number of functions are defined near the top of the script. These include a function for calculating the model solution named `calcModel`, a number of tests of the model (which are then all run from the function `allChecks` and some latex outputting utilities in `parseLaTeX`. 

Each model should then have

##### Model number and region
The model number relates to the numbers in the pdf and my notebook. Each model should have the all the inequalities that define the region in the top code block as well. 

##### Model description
A list that contains 4 colunmns. Each row relates to one integral and has 4 elements: the expression in the integral, the variable that is integrated and the lower and upper bound. All analysis and tests should derive from this list so its impossible to test one model and save output from a differently typed model.

##### Model solution
Saved to object `pModelNumber` the solved and simplified model

##### Conditions
A list of the conditions where the model is valid (i.e. the parameter space) is saved to `condModelNumber`.

##### Example values
A dictionary named `repModelNumber` with an example parameter set which is then used to test aspects of the models. The values in the dictionary should be tested against the conditions (or inequalities) in `condModelNumber`.)

##### Call to `allchecks`

This runs the following tests: 

###### Individual integrals
Total of each integral in the model solution should be positive. If this is not true either the expression in the integral is wrong or the limits are wrong or *wrong way around*.

###### Individual integrals.
Average of each integral (sum divided by difference between bounds) should be between 0 and 2r.

###### Potentiall a plot of the function. 
The results should always be increasing for either angle parameter

##### LaTeX
A LaTeX output by running `parseLaTeX`.



### Model Testing

I test each model against each other. Models that are adjacent in the modelRegions [pdf](ModelRegions.pdf) should be equal when we substitute the equation along which they are adjacent. These tests are outputted to `checksFile.tex`

All models that touch the line a=0 are tested with this value as the profile width should be 0 when a=0. Note that s=0 does not necessarily imply p=0.


### Plotting

A function is defined that calculates p given any combination of parameters. This is used to plot a full grid of values to check for odd behaviour.


### Output R function

Finally, an R function is saved from within the python script. This isn't fully automatic as the structure of the if statements and the conditions are hard to extract from the model solutions. However, automatically parsing the function solutions goes some way to avoiding typing errors.






