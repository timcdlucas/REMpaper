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

This script solves many of the more complicated REM models. It uses the symbolic algebra package sympy.

#### Models

This [pdf](ModelRegions.pdf) shows the qualitatively different model regions.

For most models in the pdf, there should be a section in the script that gives

##### Model number and region
The model number relates to the numbers in the pdf and my notebook. Each model should have the all the inequalities that define the region in the top code block as well. 

##### Model description
A list that contains 4 colunmns. Each row relates to one integral and has 4 elements: the expression in the integral, the variable that is integrated and the lower and upper bound. All analysis and tests should derive from this list so its impossible to test one model and save output from a differently typed model.

##### Model solution
Saved to object `pModelNumber` the solved and simplified model

##### Example values
A dictionary named `repModelNumber` with an example parameter set which is then used to test aspects of the models. The values in the dictionary should be tested against the inequalities.

##### Individual integrals
Total of each integral in the model solution should be positive. If this is not true either the expression in the integral is wrong or the limits are wrong or *wrong way around*.

##### Individual integrals.
Average of each integral (sum divided by difference between bounds) should be between 0 and 2r.

##### Potentiall a plot of the function. 
The results should always be increasing for either angle parameter

##### LaTeX
A LaTeX output that is currently appended to `latexOutput`



### Model Testing

I test each model against each other. Models that are adjacent in the modelRegions [pdf](ModelRegions.pdf) should be equal when we substitute the equation along which they are adjacent.







