REMpaper
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

The case for detectors with detection angles $<pi/2$ has been 
[Rowcliffe 2008](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2664.2008.01473.x/abstract)


## The code

This script solves many of the more complicated REM models. It uses the symbolic algebra 









