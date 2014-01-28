"""
Systematic analysis of REM models
Tim Lucas 
01/10/13
"""


from sympy import *
import numpy as np
import matplotlib.pyplot as pl
from datetime import datetime
import Image as Im


# Use LaTeX printing
from sympy import init_printing ;
init_printing()
# Make LaTeX output white. Because I use a dark theme
init_printing(forecolor="White") 


# Load symbols used for symbolic maths
s, a, r, g1, g2, g3, g4 = symbols('theta_s theta_a r gamma_1 gamma_2 gamma_3 gamma_4', positive=True)
r1 = {r:1} # useful for lots of checks


# A list that will recieve the latex ouput of each model

latexOutput = []
longLatexOutput = []

# Define functions to neaten up later code.

# Calculate the final profile averaged over pi.
def calcModel(model):
        x = pi**-1 * sum( [integrate(x[0], x[1:]) for x in model] ).simplify().trigsimp()
        return x

# Do the replacements fit within the area defined by the conditions? 
def confirmReplacements(conds, reps):
        if not all([c.subs(reps) for c in eval(conds)]):
                print('reps' + conds[4:] + ' incorrect')

# is average profile in range 0r-2r?
def profileRange(prof, reps):
        if not 0 <= eval(prof).subs(dict(reps, **r1)) <= 2:
                print('Total ' + prof + ' not in 0, 2r')

# Are the individuals integrals >0r
def intsPositive(model, reps):
        m = eval(model)
        for i in range(len(m)):
                if not integrate(m[i][0], m[i][1:]).subs(dict(reps, **r1)) > 0:
                    print('Integral ' + str(i+1) + ' in ' + model + ' is negative')

# Are the individual averaged integrals between 0 and  2r
def intsRange(model, reps):
        m = eval(model)
        for i in range(len(m)):
                if not 0 <= (integrate(m[i][0], m[i][1:])/(m[i][3]-m[i][2])).subs(dict(reps, **r1)) <= 2:
                        print('Integral ' + str(i+1) + ' in ' + model + ' has averaged integral outside 0<p<2r')

# Are the bounds the correct way around
def checkBounds(model, reps):
        m = eval(model)
        for i in range(len(m)):
                if not (m[i][3]-m[i][2]).subs(reps) > 0:
                        print('Bounds ' + str(i+1) + ' in ' + model + ' has lower bounds bigger than upper bounds')        

# create latex strings with the 1) the final calculated model and 2) the integral equation that defines it.
def parseLaTeX(prof):
        m = eval( 'm' + prof[1:] )    
        latexOutput.append(prof + ' = ' + latex(eval(prof)))
        longLatexOutput.append(prof + ' =\\frac{1}{\pi} \left(\int_{'+latex(m[0][2])+'}^{'+latex(m[0][3])+'}'+latex(m[0][0])+'\;\mathrm{d}'+latex(m[0][1])+'+\int_{'+latex(m[1][2])+'}^{'+latex(m[1][3])+'}'+latex(m[1][0])+'\;\mathrm{d}'+latex(m[1][1])+'+\int_{'+latex(m[2][2])+'}^{'+latex(m[2][3])+'}'+latex(m[2][0])+'\;\mathrm{d}'+latex(m[2][1])+'\\right)')

# Apply all checks.
def allChecks(prof):
        model = 'm' + prof[1:]
        reps = eval('rep' + prof[1:])
        conds = 'cond' + prof[1:]
        confirmReplacements(conds, reps)
        profileRange(prof, reps)
        intsPositive(model, reps)
        intsRange(model, reps)
        checkBounds(model, reps)

#########################################################
# 221 animal: a = 2*pi.  sensor: s > pi, a > 3pi - s  #
#########################################################



m221 = [ [2*r,                 g4, pi/2, s/2        ],
         [r + r*cos(g4 - s/2), g4, s/2,  pi         ],
         [r + r*cos(g4 + s/2), g4, pi,   2*pi-s/2   ],
         [2*r,                 g4, 2*pi-s/2, 3*pi/2 ] ]

# Replacement values in range
rep221 = {s:3*pi/2, a:2*pi} 

# Define conditions for model
cond221 = [pi <= s, a >= 3*pi - s]

# Calculate model, run checks, write output.
p221 = calcModel(m221)
allChecks('p221')
parseLaTeX('p221')

###############################################################################
# 222 animal: a > pi.  sensor: s > pi Condition: a < 3pi - s, a > 4pi - 2s  #
###############################################################################



m222 = [ [2*r,                 g4, pi/2, s/2        ],
         [r + r*cos(g4 - s/2), g4, s/2,  5*pi/2 - s/2 - a/2 ],
         [r + r*cos(g4 + s/2), g4, 5*pi/2 - s/2 - a/2,   2*pi-s/2 ],
         [2*r,                 g4, 2*pi-s/2, 3*pi/2 ] ]


# Replacement values in range
rep222 = {s:5*pi/3, a:4*pi/3} 

# Define conditions for model
cond222 = [pi <= s, a >= pi, a <= 3*pi - s, a >= 4*pi - 2*s]

# Calculate model, run checks, write output.
p222 = calcModel(m222)
allChecks('p222')
parseLaTeX('p222')


##################################################################
# 223 animal: a > pi.  sensor: s > pi Condition: a < 4pi - 2s  #
##################################################################



m223 = [ [2*r,                g4, pi/2, s/2        ],
        [r + r*cos(g4 - s/2), g4, s/2,  pi         ],
        [r + r*cos(g4 - s/2), g4, pi,   s/2 + pi/2 ],
        [r                  , g4, s/2 + pi/2,   5*pi/2 - s/2 - a/2 ],
        [r + r*cos(g4 + s/2), g4, 5*pi/2 - s/2 - a/2,   2*pi-s/2 ],
        [2*r,                 g4, 2*pi-s/2, 3*pi/2 ] ]


# Replacement values in range
rep223 = {s:5*pi/4-0.1, a:3*pi/2}

# Define conditions for model
cond223 = [pi <= s, a >= pi, a <= 4*pi - 2*s]

# Calculate model, run checks, write output.
p223 = calcModel(m223)
allChecks('p223')
parseLaTeX('p223')




########################################################
# 131 animal: a = 2*pi.   sensor:  pi/2 <= s <= pi      #
########################################################

m131 = [ [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2],
        [r - r*cos(g3 - s),    g3, 0, s - pi/2],
        [r,                    g3, s - pi/2, pi/2],
        [r - r*cos(g3),    g3, pi/2, s],
        [2*r*sin(s/2)*sin(g1),  g1, s/2, pi/2] ]

# Replacement values in range
rep131 = {s:3*pi/4} 

# Define conditions for model
cond131 = [pi/2 <= s, s <= pi]

# Calculate model, run checks, write output.
p131 = calcModel(m131)
allChecks('p131')
parseLaTeX('p131')



#################################################################################
# 231 animal: a > pi.  Sensor: pi/2 <= s <= pi. Condition: a/2 > pi - s/2  #
#################################################################################


m231 = [ [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2],
         [r - r*cos(g3-s),      g3, 0, s - pi/2],
         [r,                    g3, s - pi/2, pi/2],
         [r,                    g3, pi/2, 3*pi/2 - a/2],
         [r-r*cos(g3),          g3, 3*pi/2 - a/2, s],
         [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2 ] ]


rep231 = {s:3*pi/4, a:15*pi/8} # Replacement values in range

# Define conditions for model
cond231 = [a > pi, pi/2 <= s, s <= pi, a >= 3*pi - 2*s]

# Calculate model, run checks, write output.
p231 = calcModel(m231)
allChecks('p231')
parseLaTeX('p231')


#################################################################################
# 232 animal: a > pi.  Sensor: pi/2 <= s <= pi. Cond: 2pi - s < a < 3pi - 2s  #
#################################################################################


m232 = [ [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2],
         [r - r*cos(g3-s),      g3, 0, s - pi/2],
         [r,                    g3, s - pi/2, pi/2],
         [r,                    g3, pi/2, s],
         [r*cos(g1 - s/2),      g1, s/2, 3*pi/2 - a/2 - s/2],
         [2*r*sin(s/2)*sin(g1), g1, 3*pi/2 - a/2 - s/2, pi/2 ] ]


rep232 = {s:5*pi/8, a:6*pi/4} # Replacement values in range

# Define conditions for model
cond232 = [a > pi, pi/2 <= s, s <= pi, 2*pi - s <= a, a <= 3*pi - 2*s]

# Calculate model, run checks, write output.
p232 = calcModel(m232)
allChecks('p232')
parseLaTeX('p232')


#################################################################################
# 233 animal:  a > pi.  Sensor: pi/2 <= s <= pi. Condition: a <= 2pi - s      #
#################################################################################

m233 = [ [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2],
         [r - r*cos(g3-s),      g3, 0, s - pi/2],
         [r,                    g3, s - pi/2, pi/2],
         [r,                    g3, pi/2, s],
         [r*cos(g1 - s/2),      g1, s/2, a/2 + s/2 - pi/2] ]

rep233 = {s:3*pi/4, a:9*pi/8} # Replacement values in range

# Define conditions for model
cond233 = [a > pi,  pi/2 <= s, s <= pi, a <= 2*pi - s]

# Calculate model, run checks, write output.
p233 = calcModel(m233)
allChecks('p233')
parseLaTeX('p233')




###############################################################################
# 241 animal: a>pi.  Sensor: s <= pi/2. Condition: 2*pi - s < a             #
###############################################################################

m241 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s,          pi/2],
         [r,                    g3, 0,          s],
         [r*sin(g2),            g2, s,          pi/2],
         [r*cos(g1 - s/2),      g1, pi/2 - s/2, 3*pi/2 - s/2 - a/2],
         [2*r*sin(s/2)*sin(g1), g1, 3*pi/2 - s/2 - a/2, pi/2] ]


rep241 = {s:3*pi/8, a:29*pi/16} # Replacement values in range

# Define conditions for model
cond241 = [a >= pi, s <= pi/2, 2*pi - s <= a  ]

# Calculate model, run checks, write output.
p241 = calcModel(m241)
allChecks('p241')
parseLaTeX('p241')

####################################################################################
# 242 animal: a>pi.  Sensor: s <= pi/2. Condition:  2*pi - 2*s <= a <= 2*pi - s  #
####################################################################################

m242 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s,          pi/2],
         [r,                    g3, 0,          s],
         [r*sin(g2),            g2, s,          pi/2],
         [r*cos(g1 - s/2),      g1, pi/2 - s/2, a/2 + s/2 - pi/2] ]

rep242 = {s:3*pi/8, a:3*pi/2} # Replacement values in range

# Define conditions for model
cond242 = [a >= pi, s <= pi/2, 2*pi - 2*s <= a, a <= 2*pi - s]

# Calculate model, run checks, write output.
p242 = calcModel(m242)
allChecks('p242')
parseLaTeX('p242')


#####################################################################
# 243 animal: a>pi.  Sensor: s <= pi/2. Condition: a <= 2pi - 2s  #
#####################################################################

m243 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
          [r*sin(g2),            g2, s,          pi/2],
          [r,                    g3, 0,          s],
          [r*sin(g2),            g2, pi - a/2,   pi/2] ]


rep243 = {s:pi/9, a:10*pi/9} # Replacement values in range

# Define conditions for model
cond243 = [s <= pi/2, a >= pi, a <= 2*pi - 2*s]

# Calculate model, run checks, write output.
p243 = calcModel(m243)
allChecks('p243')
parseLaTeX('p243')



###############################################################################
# 321 animal: a <= pi.  Sensor: s > pi. Condition: a < s - pi, 4pi - 2s < a #
###############################################################################


m321 = [ [ 2*r*sin(a/2),                        g4, pi/2,         s/2 + pi/2 - a/2       ],
         [ r*sin(a/2) + r*sin(s/2 + pi/2 - g4), g4, s/2 + pi/2 - a/2, 5*pi/2 - a/2 - s/2  ], 
         [ 2*r*sin(a/2),                        g4, 5*pi/2 - a/2 - s/2,         3*pi/2] ]


rep321 = {s:19*pi/10, a:pi/2} # Replacement values in range

# Define conditions for model
cond321 = [a <= pi, s >= pi, a <= s - pi, a >= 4*pi - 2*s]

# Calculate model, run checks, write output.
p321 = calcModel(m321)
allChecks('p321')
parseLaTeX('p321')



#########################################################################################
# 322 animal: a <= pi.  Sensor: s > pi. Condition: a < s - pi, 2pi - s < a < 4pi - 2s #
#########################################################################################

m322 = [ [ 2*r*sin(a/2),                        g4, pi/2,               s/2 + pi/2 - a/2  ],
         [ r*sin(a/2) + r*cos(g4 - s/2),        g4, s/2 + pi/2 - a/2,   s/2 + pi/2        ],
         [ r*sin(a/2),                          g4, s/2 + pi/2,         5*pi/2 - a/2 - s/2],
         [ 2*r*sin(a/2),                        g4, 5*pi/2 - a/2 - s/2, 3*pi/2            ] ]

rep322 = {s:3*pi/2 + 0.1, a:pi/2} # Replacement values in range

# Define conditions for model
cond322 = [a <= pi, s >= pi, a <= s - pi, a >= 2*pi - s, a <= 4*pi - 2*s]

# Calculate model, run checks, write output.
p322 = calcModel(m322)
allChecks('p322')
parseLaTeX('p322')



###################################################################################
# 323 animal: a <= pi.  Sensor: s > pi. Condition: a <= s - pi and a < 2*pi - s #
###################################################################################

m323 = [ [ 2*r*sin(a/2),                       g4, pi/2,         s/2 + pi/2 - a/2       ],
         [ r*sin(a/2) + r*sin(s/2 + pi/2 - g4), g4, s/2 + pi/2 - a/2, s/2 + pi/2      ], 
         [ r*sin(a/2),                         g4, s/2 + pi/2,         s/2 + pi/2 + a/2] ]


rep323 = {s:3*pi/2, a:pi/3} # Replacement values in range


# Define conditions for model
cond323 = [a <= pi, s >= pi/2, a <= s - pi, a <= 2*pi - s]

# Calculate model, run checks, write output.
p323 = calcModel(m323)
allChecks('p323')
parseLaTeX('p323')


###############################################################################


"""
Complex profiles for a <= pi/2 
"""

# p-l-r for g1 profil. Calculated by AE in fig 22.4 minus AE in fig 22.3
p1 = (2*r*sin(s/4 - g1/2 + pi/4 + a/4)*sin(a/4 + pi/4 + g1/2 - s/4) - \
     2*r*sin((pi - a - 2*g1 + s)/4)*sin((pi - a + 2*g1 - s)/4)).simplify()

# p-l for g1 profiles
p2 = (2*r*sin(s/2)*sin(g1) - 2*r*sin((pi - a - 2*g1 + s)/4)*sin((pi - a + 2*g1 - s)/4)).simplify()

# p-l for g2 profile. 
p3 = (r*sin(g2) - (2*r*sin(g2/2 - a/4)*sin(pi/2 - g2/2 - a/4)).simplify()).trigsimp()



###########################################################################################
# 331 animal: a <= pi.  Sensor: pi/2 <= s <= pi. Condition: a >= s and a/2 >= s - pi/2 #
###########################################################################################


m331 =  [ [2*r*sin(s/2)*sin(g1),              g1, pi/2 - a/2 + s/2, pi/2            ],
          [p2,                                g1, s/2,              pi/2 - a/2 + s/2],
          [r*sin(a/2) + r*sin(s - pi/2 - g3), g3, 0,                s - pi/2        ],
          [r*sin(a/2),                        g3, s-pi/2,           s - pi/2 + a/2  ] ]


rep331 = {s:5*pi/8, a:6*pi/8} # Replacement values in range

# Define conditions for model
cond331 = [a <= pi, pi/2 <= s, s <= pi, a >= s, a/2 >= s - pi/2]

# Calculate model, run checks, write output.
p331 = calcModel(m331)
allChecks('p331')
parseLaTeX('p331')




##########################################################################################
# 332 animal: a <= pi.  Sensor: pi/2 <= s <= pi. Condition: a <= s and a/2 >= s- pi/2 #
##########################################################################################

m341 = [ [p1,         g1, pi/2 - s/2 + a/2, pi/2            ],
         [p2,         g1, pi/2 - s/2,       pi/2 - s/2 + a/2], # dissappears
         [p3,         g2, s,                pi/2            ],
         [r*sin(a/2), g3, 0,                a/2 + s - pi/2  ] ]


m332 =  [ [p2,                                g1, s/2,      pi/2           ],
          [r*sin(a/2) - r*cos(g3 - s),        g3, 0,        s - pi/2       ],
          [r*sin(a/2),                        g3, s - pi/2, s - pi/2 + a/2 ] ]


rep332 = {s:7*pi/8, a:7*pi/8} # Replacement values in range

# Define conditions for model
cond332 = [a <= pi, pi/2 <= s, s <= pi, a/2 <= s/2, a/2 >= s - pi/2]

# Calculate model, run checks, write output.
p332 = calcModel(m332)
allChecks('p332')
parseLaTeX('p332')







##########################################################################################
# 333 animal: a <= pi.  Sensor: pi/2 <= s <= pi. Condition: a <= s and a/2 <= s- pi/2 #
##########################################################################################



m333 = [ [p2,                                g1, s/2,            pi/2           ],
          [2*r*sin(a/2),                      g3, 0,              s - pi/2 - a/2 ],
          [r*sin(a/2) + r*sin(s - pi/2 - g3), g3, s - pi/2 - a/2, s - pi/2       ],
          [r*sin(a/2),                        g3, s - pi/2,       s - pi/2 + a/2 ] ]

rep333 = {s:7*pi/8, a:2*pi/8} # Replacement values in range

# Define conditions for model
cond333 = [a <= pi, pi/2 <= s, s <= pi, a/2 <= s/2, a/2 <= s - pi/2]

# Calculate model, run checks, write output.
p333 = calcModel(m333)
allChecks('p333')
parseLaTeX('p333')





##################################################################################
# 341 animal: a <= pi.  Sensor: s <= pi/2. Condition: a > pi - 2s &  a <= s      #
##################################################################################


m341 = [ [p1,         g1, pi/2 - s/2 + a/2, pi/2            ],
         [p2,         g1, pi/2 - s/2,       pi/2 - s/2 + a/2],
         [p3,         g2, s,                pi/2            ],
         [r*sin(a/2), g3, 0,                a/2 + s - pi/2  ] ]

rep341 = {s:pi/2-0.1, a:pi/4} # Replacement values in range

# Define conditions for model
cond341 = [a <= pi,  s <= pi/2,  a >= pi - 2*s,  a <= s]

# Calculate model, run checks, write output.
p341 = calcModel(m341)
allChecks('p341')
parseLaTeX('p341')


######################################################################################
# 342 animal: a <= pi.  Sensor: s <= pi/2. Condition: a > pi - 2s &  s <= a <= 2s  #
######################################################################################


m342 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 + s/2 - a/2, pi/2            ],
         [p2,                   g1, pi/2 - s/2,       pi/2 + s/2 - a/2],
         [p3,                   g2, s,                pi/2        ],
         [r*sin(a/2),           g3, 0,                a/2 + s -pi/2   ] ]


rep342 = {s:pi/2-0.1, a:pi/2} # Replacement values in range

# define conditions for model
cond342 = [a <= pi,  s <= pi/2,  a >= pi - 2*s,  s <= a, a <= 2*s]


# Calculate model, run checks, write output.
p342 = calcModel(m342)
allChecks('p342')
parseLaTeX('p342')





##################################################################################
# 343 animal: a <= pi.  Sensor: s <= pi/2. Condition: a > pi - 2s &  a > 2s      #
##################################################################################



m343 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2            ],
         [r*sin(g2),            g2, s,          a/2             ],
         [p3,                   g2, a/2,        pi/2            ],
         [r*sin(a/2),           g3, 0,          a/2 + s -pi/2   ] ]


rep343 = {s:pi/4, a:3*pi/4} # Replacement values in range


# Define conditions for model
cond343 = [a <= pi,  s <= pi/2,  a >= pi - 2*s,  a > 2*s]

# Calculate model, run checks, write output.
p343 = calcModel(m343)
allChecks('p343')
parseLaTeX('p343')




###################################################################################
# 344 animal: a <= pi.  Sensor: s <= pi/2. Condition: a <= pi - 2s & a <= s     #
###################################################################################

m344 = [ [p1, g1, pi/2 - s/2 + a/2, pi/2            ],
         [p2, g1, pi/2 - s/2,       pi/2 - s/2 + a/2],
         [p3, g2, s,                s + a/2         ] ]


rep344 = {s:2*pi/8, a:pi/8} # Replacement values in range

# Define conditions for model
cond344 = [a <= pi, s <= pi/2, a <= pi - 2*s, a <= s]

# Calculate model, run checks, write output.
p344 = calcModel(m344)
allChecks('p344')
parseLaTeX('p344')




#######################################################################################
# 345 animal: a <= pi.  Sensor: s <= pi/2. Condition: a <= pi - 2s & s <= a <= 2s   #
#######################################################################################


m345 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 + s/2 - a/2, pi/2            ],
         [p2,                   g1, pi/2 - s/2,       pi/2 + s/2 - a/2],
         [p3,                   g2, s,                s + a/2         ] ]

rep345 = {s:2*pi/8, a:pi/2-0.1} # Replacement values in range

# Define conditions for model
cond345 = [a <= pi, s <= pi/2, a <= pi - 2*s, s <= a, a <= 2*s]

# Calculate model, run checks, write output.
p345 = calcModel(m345)
allChecks('p345')
parseLaTeX('p345')



##################################################################################
# 346 animal: a <= pi.  Sensor: s <= pi/2. Condition: a <= pi - 2s &  2s <= a      #
##################################################################################

m346 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2    ],
         [r*sin(g2),            g2, s,          a/2     ],
         [p3,                   g2, a/2,        s + a/2 ] ]


rep346 = {s:1*pi/8, a:pi/2} # Replacement values in range

# Define conditions for model
cond346 = [a <= pi,  s <= pi/2,  a <= pi - 2*s,  2*s <= a]

# Calculate model, run checks, write output.
p346 = calcModel(m346)
allChecks('p346')
parseLaTeX('p346')




##################################################################################################################

####################
## Run tests     ###
####################

# create gas model object
gas = 2*r
p311 = 2*r*sin(a/2)
p141 = r*(2+s)/pi


# for each model run through every adjacent model. 
# Contains duplicatea but better for avoiding missed comparisons.
# Also contains replacement s->a and a->s just in case. 


allComps = [
['gas', 'p221', {s:2*pi}],
['gas', 'p311', {a:pi}],

['p221', 'gas', {s:2*pi}],
['p221', 'p131', {s:pi}],
['p221', 'p222',{a:3*pi-s}],
['p221', 'p222',{s:3*pi-a}],

['p222', 'p221',{a:3*pi-s}],
['p222', 'p221',{s:3*pi-a}],
['p222', 'p223',{a:4*pi-2*s}],
['p222', 'p223',{s:2*pi-a/2}],
['p222', 'p321',{a:pi}],

['p223', 'p222',{a:4*pi-2*s}],
['p223', 'p222',{s:2*pi-a/2}],
['p223', 'p322',{a:pi}],
['p223', 'p231',{s:pi}],

['p131','p221', {s:pi}],
['p131','p231',{a:2*pi}],

['p231','p223',{s:pi}],
['p231','p232',{a:3*pi-2*s}],
['p231','p232',{s:3*pi/2-a/2}],
['p231','p131',{a:2*pi}],

['p232','p241',{s:pi/2}],
['p232','p233',{a:2*pi-s}],
['p232','p233',{s:2*pi-a}],
['p232','p231',{a:3*pi-2*s}],
['p232','p231',{s:3*pi/2-a/2}],

['p233','p242',{s:pi/2}],
['p233','p232',{s:2*pi-a}],
['p233','p232',{a:2*pi-s}],
['p233','p331',{a:pi}],

['p141','p131', {s:pi/2}],
['p141','p241',{a:2*pi}],

['p241','p141',{a:2*pi}],
['p241','p242',{a:2*pi-s}],
['p241','p242',{s:2*pi-a}],
['p241','p232',{s:pi/2}],

['p242','p241',{a:2*pi-s}], 
['p242','p241',{s:2*pi-a}],
['p242','p243',{s:pi-a/2}],
['p242','p243',{a:2*pi-2*s}],
['p241','p233',{s:pi/2}],

['p243','p242',{s:2*pi-2*a}],
['p243','p242',{a:2*pi-2*s}],
['p243','p343',{a:pi}],

['p311','p321',{s:2*pi}],
['p311','gas',{a:pi}],

['p321','p322',{s:2*pi-a/2}],
['p321','p322',{a:4*pi-2*s}],
['p321','p311',{s:2*pi}],
['p321','p222',{a:pi}],

['p322','p321',{a:4*pi-2*s}],
['p322','p321',{s:2*pi-a/2}],
['p322','p323',{a:2*pi-s}],
['p322','p323',{s:2*pi-a}],
['p322','p223',{a:pi}],

['p323','p322',{s:2*pi-a}],
['p323','p322',{a:2*pi-s}],
['p323','p333',{s:pi}],

['p331','p342',{s:pi/2}],
['p331','p332',{a:s}],
['p331','p332',{s:a}],
['p331','p233',{a:pi}],

['p332','p331',{a:s}],
['p332','p331',{s:a}],
['p332','p341',{s:pi/2}],
['p332','p333',{a:2*s-pi}],
['p332','p333',{s:a/2+pi/2}],

['p333','p332',{s:a/2+pi/2}],
['p333','p332',{a:2*s-pi}],
['p333','p323',{s:pi}],


['p341','p344',{a:pi-2*s}],
['p341','p344',{s:pi/2-a/2}],
['p341','p342',{s:a}],
['p341','p342',{a:s}],
['p341','p332',{s:pi/2}],

['p342','p341',{s:a}],
['p342','p341',{a:s}],
['p342','p345',{s:pi/2-a/2}],
['p342','p345',{a:pi-2*s}],
['p342','p343',{a:2*s}],
['p342','p343',{s:a/2}],
['p342','p331',{s:pi/2}],

['p343','p346',{s:pi/2-a/2}],
['p343','p346',{a:pi-2*s}],
['p343','p342',{a:2*s}],
['p343','p342',{s:a/2}],
['p343','p243',{a:pi}],


['p344','p345',{s:a}],
['p344','p345',{a:s}],
['p344','p341',{s:pi/2-a/2}],
['p344','p341',{a:pi-2*s}],

['p345','p344',{a:s}],
['p345','p344',{s:a}],
['p345','p346',{a:2*s}],
['p345','p346',{s:a/2}],
['p345','p342',{a:pi-2*s}],
['p345','p342',{s:pi/2-a/2}],

['p346','p345',{a:2*s}],
['p346','p345',{s:a/2}],
['p346','p343',{a:pi-2*s}],
['p346','p343',{s:pi/2-a/2}]
]


# List of regions that cover a=0. Should equal 0 when a=0.
zeroRegions = ['p346', 'p345', 'p344', 'p341', 'p332', 'p333', 'p323', 'p322', 'p321', 'p311']

# Run through all the comparisons. Need simplify(). Even together() gives some false negatives.

checkFile = open('/home/tim/Dropbox/PhD/Analysis/REM-chapter/checksFile.tex','w')

checkFile.write('All checks evaluated.\nTim Lucas - ' + str(datetime.now()) + '\n')
for i in range(len(allComps)):
        if (eval(allComps[i][0]).subs(allComps[i][2]) - eval(allComps[i][1]).subs(allComps[i][2])).simplify() == 0:
                checkFile.write(str(i) + ': ' + allComps[i][0]+ ' and ' +allComps[i][1]+': OK\n')
        else:
                checkFile.write(str(i) + ': ' + allComps[i][0]+ ' and ' +allComps[i][1]+': Incorrect\n')

for i in range(len(zeroRegions)):
        if eval(zeroRegions[i]).subs({a:0}).simplify() == 0:
                checkFile.write(zeroRegions[i] + ' at a=0: OK\n')
        else:
                checkFile.write(zeroRegions[i] + ' at a=0: Incorrect\n')

checkFile.close()


# And print to terminal
#for i in range(len(allComps)):
#        if not (eval(allComps[i][0]).subs(allComps[i][2]) - eval(allComps[i][1]).subs(allComps[i][2])).simplify() == 0:
#               print allComps[i][0] + ' and ' + allComps[i][1]+': Incorrect\n'

#####################################
## Check some that don't work well ##
#####################################

xRange = np.arange(0,pi/2, 0.01)
y332Range = [p332.subs({r:1, s:pi/2, a:i}).n() for i in xRange]
plot332 = pl.plot(xRange, y332Range)
pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p332Profile.pdf')
pl.close()

y341Range = [p341.subs({r:1, s:pi/2, a:i}).n() for i in xRange]
plot341 = pl.plot(xRange, y341Range)
pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p341Profile.pdf')
pl.close()



#pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p221Profile.pdf')
#pl.close()



######################
### Write output   ###
######################

# write out full latex model solutions and model statements

latexFile = open('/home/tim/Dropbox/PhD/Analysis/REM-chapter/ModelSolutions.tex', 'w')

latexFile.write('% LaTeX output. Solutions of all REM models.\n' + '%Tim Lucas - ' + str(datetime.now()) + '\n')
for i in range(len(latexOutput)):
        latexFile.write( '\\[' + latexOutput[i] + '\\]\n')

latexFile.close()

latexFile = open('/home/tim/Dropbox/PhD/Analysis/REM-chapter/ModelDefinitions.tex', 'w')

latexFile.write('% LaTeX output. Definitions of all REM models.\n' + '%Tim Lucas - ' + str(datetime.now()) + '\n')
for i in range(len(longLatexOutput)):
        latexFile.write( '\\[' + longLatexOutput[i] + '\\]\n')

latexFile.close()





