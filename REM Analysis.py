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

#########################################################
# 221 animal: a = 2*pi.  sensor: s > pi, a > 3pi - s  #
#########################################################



m221 = [ [2*r,                 g4, pi/2, s/2        ],
         [r + r*cos(g4 - s/2), g4, s/2,  pi         ],
         [r + r*cos(g4 + s/2), g4, pi,   2*pi-s/2   ],
         [2*r,                 g4, 2*pi-s/2, 3*pi/2 ] ]

p221 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m221] ).simplify().trigsimp()




# Replacement values in range
rep221 = {s:3*pi/2, a:2*pi} 

# Define conditions for model
cond221 = [pi <= s, a >= 3*pi - s]
# Confirm replacements
if not all([c.subs(rep221) for c in cond221]):
        print('rep221 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p221.subs(dict(rep221, **r1)) <= 2:
        print('Total p221 not in 0, 2r')

# Are the individuals integrals >0r
for i in range(len(m221)):
        if not integrate(m221[i][0], m221[i][1:]).subs(dict(rep221, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p221 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m221)):
        if not 0 <= (integrate(m221[i][0], m221[i][1:])/(m221[i][3]-m221[i][2])).subs(dict(rep221, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p221 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m221)):
        if not (m221[i][3]-m221[i][2]).subs(rep221) > 0:
                print('Bounds ' + str(i+1) + ' in p221 has lower bounds bigger than upper bounds')
        
# Plot function

xRange = np.arange(pi,2*pi, 0.01)
yRange = [p221.subs({r:1, s:i}).n() for i in xRange]
#plot221 = pl.plot(xRange, yRange)
#pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p221Profile.pdf')
#pl.close()

# LaTeX output

latexOutput.append('p221 &= ' + latex(p221))

longLatexOutput.insert(0,'p221=&\\frac{1}{\pi} \left(2\int_{'+latex(m221[0][2])+'}^{'+latex(m221[0][3])+'}'+latex(m221[0][0])+'\;\mathrm{d}'+latex(m221[0][1])+'+2\int_{'+latex(m221[1][2])+'}^{'+latex(m221[1][3])+'}'+latex(m221[1][0])+'\;\mathrm{d}'+latex(m221[1][1])+'\\right)')



###############################################################################
# 222 animal: a > pi.  sensor: s > pi Condition: a < 3pi - s, a > 4pi - 2s  #
###############################################################################



m222 = [ [2*r,                 g4, pi/2, s/2        ],
         [r + r*cos(g4 - s/2), g4, s/2,  5*pi/2 - s/2 - a/2 ],
         [r + r*cos(g4 + s/2), g4, 5*pi/2 - s/2 - a/2,   2*pi-s/2 ],
         [2*r,                 g4, 2*pi-s/2, 3*pi/2 ] ]

p222 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m222] ).simplify().trigsimp()

# Replacement values in range
rep222 = {s:5*pi/3, a:4*pi/3} 

# Define conditions for model
cond222 = [pi <= s, a >= pi, a <= 3*pi - s, a >= 4*pi - 2*s]
## Confirm replacementsa
if not all([c.subs(rep222) for c in cond222]):
        print('rep222 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p222.subs(dict(rep222, **r1)) <= 2:
        print('Total p222 not in 0, 2r')

# Are the individuals integrals >0r
for i in range(len(m222)):
        if not integrate(m222[i][0], m222[i][1:]).subs(dict(rep222, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p222 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m222)):
        if not 0 <= (integrate(m222[i][0], m222[i][1:])/(m222[i][3]-m222[i][2])).subs(dict(rep222, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p222 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m222)):
        if not (m222[i][3]-m222[i][2]).subs(rep222) > 0:
                print('Bounds ' + str(i+1) + ' in p222 has lower bounds bigger than upper bounds')
        
# Plot function

xRange = np.arange(pi,2*pi, 0.01)
yRange = [p222.subs({r:1, s:i}).n() for i in xRange]
#plot222 = pl.plot(xRange, yRange)
#pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p222Profile.pdf')
#pl.close()

# LaTeX output

latexOutput.append('p222 &= ' + latex(p222))
longLatexOutput.append('p222=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m222]).lstrip('+') + '\\right)')


##################################################################
# 223 animal: a > pi.  sensor: s > pi Condition: a < 4pi - 2s  #
##################################################################



m223 = [ [2*r,                g4, pi/2, s/2        ],
        [r + r*cos(g4 - s/2), g4, s/2,  pi         ],
        [r + r*cos(g4 - s/2), g4, pi,   s/2 + pi/2 ],
        [r                  , g4, s/2 + pi/2,   5*pi/2 - s/2 - a/2 ],
        [r + r*cos(g4 + s/2), g4, 5*pi/2 - s/2 - a/2,   2*pi-s/2 ],
        [2*r,                 g4, 2*pi-s/2, 3*pi/2 ] ]

p223 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m223] ).simplify().trigsimp()

# Replacement values in range
rep223 = {s:5*pi/4-0.1, a:3*pi/2}

# Define conditions for model
cond223 = [pi <= s, a >= pi, a <= 4*pi - 2*s]
## Confirm replacementsa
if not all([c.subs(rep223) for c in cond223]):
    print('rep223 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p223.subs(dict(rep223, **r1)) <= 2:
    print('Total p223 not in 0, 2r')

# Are the individuals integrals >0r
for i in range(len(m223)):
    if not integrate(m223[i][0], m223[i][1:]).subs(dict(rep223, **r1)) > 0:
        print('Integral ' + str(i+1) + ' in p223 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m223)):
    if not 0 <= (integrate(m223[i][0], m223[i][1:])/(m223[i][3]-m223[i][2])).subs(dict(rep223, **r1)) <= 2:
        print('Integral ' + str(i+1) + ' in p223 has averaged integral outside 0<p<2r')

# Are the bounds the correct way around
for i in range(len(m223)):
    if not (m223[i][3]-m223[i][2]).subs(rep223) > 0:
        print('Bounds ' + str(i+1) + ' in p223 has lower bounds bigger than upper bounds')

# Plot function

xRange = np.arange(pi,2*pi, 0.01)
yRange = [p223.subs({r:1, s:i}).n() for i in xRange]
#plot223 = pl.plot(xRange, yRange)
#pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p223Profile.pdf')
#pl.close()

# LaTeX output

latexOutput.append('p223 &= ' + latex(p223))
longLatexOutput.append('p223=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m223]).lstrip('+') + '\\right)')




########################################################
# 131 animal: a = 2*pi.   sensor:  pi/2 <= s <= pi      #
########################################################

m131 = [ [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2],
        [r - r*cos(g3 - s),    g3, 0, s - pi/2],
        [r,                    g3, s - pi/2, pi/2],
        [r - r*cos(g3),    g3, pi/2, s],
        [2*r*sin(s/2)*sin(g1),  g1, s/2, pi/2] ]

p131 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m131] ).simplify().trigsimp()


# Replacement values in range
rep131 = {s:3*pi/4} 

# Define conditions for model
cond131 = [pi/2 <= s, s <= pi]
# Confirm replacements
if not all([c.subs(rep131) for c in cond131]):
        print('rep131 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p131.subs(dict(rep131, **r1)) <= 2:
        print('Total p131 not in 0, 2r')

# Are the individuals integrals >0r
for i in range(len(m131)):
        if not integrate(m131[i][0], m131[i][1:]).subs(dict(rep131, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p131 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m131)):
        if not 0 <= (integrate(m131[i][0], m131[i][1:])/(m131[i][3]-m131[i][2])).subs(dict(rep131, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p131 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m131)):
        if not (m131[i][3]-m131[i][2]).subs(rep131) > 0:
                print('Bounds ' + str(i+1) + ' in p131 has lower bounds bigger than upper bounds')        


# Plot function

xRange = np.arange(pi/2,pi, 0.01)
yRange = [p131.subs({r:1, s:i}).n() for i in xRange]
#plot131 = pl.plot(xRange, yRange)
#pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p131Profile.pdf')
#pl.close()

# LaTeX output

latexOutput.append('p131 &= ' + latex(p131))
longLatexOutput.append('p131=&\\frac{1}{\pi} \left(\int_{'+latex(m131[0][2])+'}^{'+latex(m131[0][3])+'}'+latex(m131[0][0])+'\;\mathrm{d}'+latex(m131[0][1])+'+\int_{'+latex(m131[1][2])+'}^{'+latex(m131[1][3])+'}'+latex(m131[1][0])+'\;\mathrm{d}'+latex(m131[1][1])+'+\int_{'+latex(m131[2][2])+'}^{'+latex(m131[2][3])+'}'+latex(m131[2][0])+'\;\mathrm{d}'+latex(m131[2][1])+'\\right)')



#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Two models in 7.8.
#animal: pi <= a .  Sensor: pi/2 <= s <= pi
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#################################################################################
# 233 animal:  a > pi.  Sensor: pi/2 <= s <= pi. Condition: a <= 2pi - s      #
#################################################################################

m233 = [ [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2],
         [r - r*cos(g3-s),      g3, 0, s - pi/2],
         [r,                    g3, s - pi/2, pi/2],
         [r,                    g3, pi/2, s],
         [r*cos(g1 - s/2),      g1, s/2, a/2 + s/2 - pi/2],
         [0,                    g1, a/2 + s/2 - pi/2, pi/2 ] ]

p233 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m233] ).simplify().trigsimp()

rep233 = {s:3*pi/4, a:9*pi/8} # Replacement values in range

# Define conditions for model
cond233 = [a > pi,  pi/2 <= s, s <= pi, a <= 2*pi - s]
# Confirm replacements
if not all([c.subs(rep233) for c in cond233]):
        print('rep233 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p233.subs(dict(rep233, **r1)) <= 2:
        print('Total p233 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m233)):
        if not integrate(m233[i][0], m233[i][1:]).subs(dict(rep233, **r1)) >= 0:
                print('Integral ' + str(i+1) + ' in p233 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m233)):
        if not 0 <= (integrate(m233[i][0], m233[i][1:])/(m233[i][3]-m233[i][2])).subs(dict(rep233, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p233 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m233)):
        if not (m233[i][3]-m233[i][2]).subs(rep233) >= 0:
                print('Bounds ' + str(i+1) + ' in p233 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p233 &= ' + latex(p233))
longLatexOutput.append('p233=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m233]).lstrip('+') + '\\right)')



#################################################################################
# 232 animal: a > pi.  Sensor: pi/2 <= s <= pi. Cond: 2pi - s < a < 3pi - 2s  #
#################################################################################


m232 = [ [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2],
         [r - r*cos(g3-s),      g3, 0, s - pi/2],
         [r,                    g3, s - pi/2, pi/2],
         [r,                    g3, pi/2, s],
         [r*cos(g1 - s/2),      g1, s/2, 3*pi/2 - a/2 - s/2],
         [2*r*sin(s/2)*sin(g1), g1, 3*pi/2 - a/2 - s/2, pi/2 ] ]


p232 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m232] ).simplify().trigsimp()



rep232 = {s:5*pi/8, a:6*pi/4} # Replacement values in range

# Define conditions for model
cond232 = [a > pi, pi/2 <= s, s <= pi, 2*pi - s <= a, a <= 3*pi - 2*s]
# Confirm replacements
if not all([c.subs(rep232) for c in cond232]):
        print('rep232 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p232.subs(dict(rep232, **r1)) <= 2:
        print('Total p232 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m232)):
        if not integrate(m232[i][0], m232[i][1:]).subs(dict(rep232, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p232 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m232)):
        if not 0 <= (integrate(m232[i][0], m232[i][1:])/(m232[i][3]-m232[i][2])).subs(dict(rep232, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p232 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m232)):
        if not (m232[i][3]-m232[i][2]).subs(rep232) > 0:
                print('Bounds ' + str(i+1) + ' in p232 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p232 &= ' + latex(p232))
longLatexOutput.append('p232=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m232]).lstrip('+') + '\\right)')




#################################################################################
# 231 animal: a > pi.  Sensor: pi/2 <= s <= pi. Condition: a/2 > pi - s/2  #
#################################################################################


m231 = [ [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2],
         [r - r*cos(g3-s),      g3, 0, s - pi/2],
         [r,                    g3, s - pi/2, pi/2],
         [r,                    g3, pi/2, 3*pi/2 - a/2],
         [r-r*cos(g3),          g3, 3*pi/2 - a/2, s],
         [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2 ] ]


p231 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m231] ).simplify().trigsimp()



rep231 = {s:3*pi/4, a:15*pi/8} # Replacement values in range

# Define conditions for model
cond231 = [a > pi, pi/2 <= s, s <= pi, a >= 3*pi - 2*s]

# Confirm replacements
if not all([c.subs(rep231) for c in cond231]):
        print('rep231 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p231.subs(dict(rep231, **r1)) <= 2:
        print('Total p231 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m231)):
        if not integrate(m231[i][0], m231[i][1:]).subs(dict(rep231, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p231 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m231)):
        if not 0 <= (integrate(m231[i][0], m231[i][1:])/(m231[i][3]-m231[i][2])).subs(dict(rep231, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p231 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m231)):
        if not (m231[i][3]-m231[i][2]).subs(rep231) > 0:
                print('Bounds ' + str(i+1) + ' in p231 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p231 &= ' + latex(p231))
longLatexOutput.append('p231=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m231]).lstrip('+') + '\\right)')




#####################################################################
# 243 animal: a>pi.  Sensor: s <= pi/2. Condition: a <= 2pi - 2s  #
#####################################################################

m243 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
          [r*sin(g2),            g2, s,          pi/2],
          [r,                    g3, 0,          s],
          [r*sin(g2),            g2, pi - a/2,   pi/2] ]


p243 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m243] ).simplify().trigsimp()


rep243 = {s:pi/9, a:10*pi/9} # Replacement values in range

# Define conditions for model
cond243 = [s <= pi/2, a >= pi, a <= 2*pi - 2*s]
# Confirm replacements
if not all([c.subs(rep243) for c in cond243]):
        print('rep243 incorrect')
# is average profile in range 0r-2r?
if not 0 <= p243.subs(dict(rep243, **r1)) <= 2:
        print('Total p243 not in 0, 2r')

# Are the individuals integrals >0r
for i in range(len(m243)):
        if not integrate(m243[i][0], m243[i][1:]).subs(dict(rep243, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p243 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m243)):
        if not 0 <= (integrate(m243[i][0], m243[i][1:])/(m243[i][3]-m243[i][2])).subs(dict(rep243, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p243 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m243)):
        if not (m243[i][3]-m243[i][2]).subs(rep243) > 0:
                print('Bounds ' + str(i+1) + ' in p243 has lower bounds bigger than upper bounds')        


# LaTeX output

latexOutput.append('p243 &= ' + latex(p243))
longLatexOutput.append('p243=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m243]).lstrip('+') + '\\right)')


####################################################################################
# 242 animal: a>pi.  Sensor: s <= pi/2. Condition:  2*pi - 2*s <= a <= 2*pi - s  #
####################################################################################

m242 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s,          pi/2],
         [r,                    g3, 0,          s],
         [r*sin(g2),            g2, s,          pi/2],
         [r*cos(g1 - s/2),      g1, pi/2 - s/2, a/2 + s/2 - pi/2] ]


p242 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m242] ).simplify().trigsimp()


rep242 = {s:3*pi/8, a:3*pi/2} # Replacement values in range

# Define conditions for model
cond242 = [a >= pi, s <= pi/2, 2*pi - 2*s <= a, a <= 2*pi - s]
# Confirm replacements
if not all([c.subs(rep242) for c in cond242]):
        print('rep242 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p242.subs(dict(rep242, **r1)) <= 2:
        print('Total p242 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m242)):
        if not integrate(m242[i][0], m242[i][1:]).subs(dict(rep242, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p242 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m242)):
        if not 0 <= (integrate(m242[i][0], m242[i][1:])/(m242[i][3]-m242[i][2])).subs(dict(rep242, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p242 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m242)):
        if not (m242[i][3]-m242[i][2]).subs(rep242) > 0:
                print('Bounds ' + str(i+1) + ' in p242 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p242 &= ' + latex(p242))
longLatexOutput.append('p242 &= \\frac{1}{\pi} \left(\int_{'+latex(m242[0][2])+'}^{'+latex(m242[0][3])+'}'+latex(m242[0][0])+'\;\mathrm{d}'+latex(m242[0][1])+'+\int_{'+latex(m242[1][2])+'}^{'+latex(m242[1][3])+'}'+latex(m242[1][0])+'\;\mathrm{d}'+latex(m242[1][1])+'+\int_{'+latex(m242[2][2])+'}^{'+latex(m242[2][3])+'}'+latex(m242[2][0])+'\;\mathrm{d}'+latex(m242[2][1])+'+\int_{'+latex(m242[3][2])+'}^{'+latex(m242[3][3])+'}'+latex(m242[3][0])+'\;\mathrm{d}'+latex(m242[3][1])+'\\right)')


###############################################################################
# 241 animal: a>pi.  Sensor: s <= pi/2. Condition: 2*pi - s < a             #
###############################################################################

m241 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s,          pi/2],
         [r,                    g3, 0,          s],
         [r*sin(g2),            g2, s,          pi/2],
         [r*cos(g1 - s/2),      g1, pi/2 - s/2, 3*pi/2 - s/2 - a/2],
         [2*r*sin(s/2)*sin(g1), g1, 3*pi/2 - s/2 - a/2, pi/2] ]

p241 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m241] ).simplify().trigsimp()


rep241 = {s:3*pi/8, a:29*pi/16} # Replacement values in range

# Define conditions for model
cond241 = [a >= pi, s <= pi/2, 2*pi - s <= a  ]
# Confirm replacements
if not all([c.subs(rep241) for c in cond241]):
        print('rep241 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p241.subs(dict(rep241, **r1)) <= 2:
        print('Total p241 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m241)):
        if not integrate(m241[i][0], m241[i][1:]).subs(dict(rep241, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p241 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m241)):
        if not 0 <= (integrate(m241[i][0], m241[i][1:])/(m241[i][3]-m241[i][2])).subs(dict(rep241, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p241 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m241)):
        if not (m241[i][3]-m241[i][2]).subs(rep241) > 0:
                print('Bounds ' + str(i+1) + ' in p241 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p241 &= ' + latex(p241))
longLatexOutput.append('p241 &= \\frac{1}{\pi} \left(\int_{'+latex(m241[0][2])+'}^{'+latex(m241[0][3])+'}'+latex(m241[0][0])+'\;\mathrm{d}'+latex(m241[0][1])+'+\int_{'+latex(m241[1][2])+'}^{'+latex(m241[1][3])+'}'+latex(m241[1][0])+'\;\mathrm{d}'+latex(m241[1][1])+'+\int_{'+latex(m241[2][2])+'}^{'+latex(m241[2][3])+'}'+latex(m241[2][0])+'\;\mathrm{d}'+latex(m241[2][1])+'+\int_{'+latex(m241[3][2])+'}^{'+latex(m241[3][3])+'}'+latex(m241[3][0])+'\;\mathrm{d}'+latex(m241[3][1])+'\\right)')



###############################################################################
# 321 animal: a <= pi.  Sensor: s > pi. Condition: a < s - pi, 4pi - 2s < a #
###############################################################################


m321 = [ [ 2*r*sin(a/2),                        g4, pi/2,         s/2 + pi/2 - a/2       ],
         [ r*sin(a/2) + r*sin(s/2 + pi/2 - g4), g4, s/2 + pi/2 - a/2, 5*pi/2 - a/2 - s/2  ], 
         [ 2*r*sin(a/2),                        g4, 5*pi/2 - a/2 - s/2,         3*pi/2] ]

p321 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m321] ).simplify().trigsimp()



rep321 = {s:19*pi/10, a:pi/2} # Replacement values in range

# Define conditions for model
cond321 = [a <= pi, s >= pi, a <= s - pi, a >= 4*pi - 2*s]
# Confirm replacements
if not all([c.subs(rep321) for c in cond321]):
        print('rep321 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p321.subs(dict(rep321, **r1)) <= 2:
        print('Total p321 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m321)):
        if not integrate(m321[i][0], m321[i][1:]).subs(dict(rep321, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p321 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m321)):
        if not 0 <= (integrate(m321[i][0], m321[i][1:])/(m321[i][3]-m321[i][2])).subs(dict(rep321, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p321 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m321)):
        if not (m321[i][3]-m321[i][2]).subs(rep321) > 0:
                print('Bounds ' + str(i+1) + ' in p321 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p321 &= ' + latex(p321))
longLatexOutput.append('p321=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m321]).lstrip('+') + '\\right)')


#########################################################################################
# 322 animal: a <= pi.  Sensor: s > pi. Condition: a < s - pi, 2pi - s < a < 4pi - 2s #
#########################################################################################




m322 = [ [ 2*r*sin(a/2),                        g4, pi/2,               s/2 + pi/2 - a/2  ],
         [ r*sin(a/2) + r*cos(g4 - s/2),        g4, s/2 + pi/2 - a/2,   s/2 + pi/2        ],
         [ r*sin(a/2),                          g4, s/2 + pi/2,         5*pi/2 - a/2 - s/2],
         [ 2*r*sin(a/2),                        g4, 5*pi/2 - a/2 - s/2, 3*pi/2            ] ]


p322 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m322] ).simplify().trigsimp()



rep322 = {s:3*pi/2 + 0.1, a:pi/2} # Replacement values in range

# Define conditions for model
cond322 = [a <= pi, s >= pi, a <= s - pi, a >= 2*pi - s, a <= 4*pi - 2*s]
# Confirm replacements
if not all([c.subs(rep322) for c in cond322]):
        print('rep322 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p322.subs(dict(rep322, **r1)) <= 2:
        print('Total p322 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m322)):
        if not integrate(m322[i][0], m322[i][1:]).subs(dict(rep322, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p322 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m322)):
        if not 0 <= (integrate(m322[i][0], m322[i][1:])/(m322[i][3]-m322[i][2])).subs(dict(rep322, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p322 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m322)):
        if not (m322[i][3]-m322[i][2]).subs(rep322) > 0:
                print('Bounds ' + str(i+1) + ' in p322 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p322 &= ' + latex(p322))
longLatexOutput.append('p322=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m322]).lstrip('+') + '\\right)')




###################################################################################
# 323 animal: a <= pi.  Sensor: s > pi. Condition: a <= s - pi and a < 2*pi - s #
###################################################################################

m323 = [ [ 2*r*sin(a/2),                       g4, pi/2,         s/2 + pi/2 - a/2       ],
         [ r*sin(a/2) + r*sin(s/2 + pi/2 - g4), g4, s/2 + pi/2 - a/2, s/2 + pi/2      ], 
         [ r*sin(a/2),                         g4, s/2 + pi/2,         s/2 + pi/2 + a/2] ]
p323 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m323] ).simplify().trigsimp()



rep323 = {s:3*pi/2, a:pi/3} # Replacement values in range


# Define conditions for model
cond323 = [a <= pi, s >= pi/2, a <= s - pi, a <= 2*pi - s]
# Confirm replacements
if not all([c.subs(rep323) for c in cond323]):
        print('rep323 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p323.subs(dict(rep323, **r1)) <= 2:
        print('Total p323 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m323)):
        if not integrate(m323[i][0], m323[i][1:]).subs(dict(rep323, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p323 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m323)):
        if not 0 <= (integrate(m323[i][0], m323[i][1:])/(m323[i][3]-m323[i][2])).subs(dict(rep323, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p323 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m323)):
        if not (m323[i][3]-m323[i][2]).subs(rep323) > 0:
                print('Bounds ' + str(i+1) + ' in p323 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p323 &= ' + latex(p323))
longLatexOutput.append('p323=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m323]).lstrip('+') + '\\right)')





###############################################################################


"""
Complex profiles for a <= pi/2 
"""

# p-l-r for g1 profil. Calculated by AE in fig 22.4 minus AE in fig 22.3
p1 = (2*r*sin(s - 3*pi/2 + g1)*sin((g1 - pi/2 - s + a)/2) - \
     2*r*sin((pi - a - 2*g1 + s)/4)*sin((pi - a + 2*g1 - s)/4)).simplify()

# p-l for g1 profiles
p2 = (2*r*sin(s/2)*sin(g1) - 2*r*sin((pi - a - 2*g1 + s)/4)*sin((pi - a + 2*g1 - s)/4)).simplify()

# p-l for g2 profile. 
p3 = (r*sin(g2) - 2*r*sin(g2/2 - a/4)*sin((pi - a + 2*g2 - s)/4)).simplify()



###########################################################################################
# 331 animal: a <= pi.  Sensor: pi/2 <= s <= pi. Condition: a >= s and a/2 >= s - pi/2 #
###########################################################################################


m331 = [ [2*r*sin(s/2)*sin(g1),              g1, pi/2 - a/2 + s/2, pi/2            ],
          [p2,                                g1, s/2,              pi/2 - a/2 + s/2],
          [r*sin(a/2) + r*sin(s - pi/2 - g3), g3, 0,                s - pi/2        ],
          [r*sin(a/2),                        g3, s-pi/2,           s - pi/2 + a/2  ] ]

p331 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m331] ).simplify().trigsimp()



rep331 = {s:5*pi/8, a:6*pi/8} # Replacement values in range

# Define conditions for model
cond331 = [a <= pi, pi/2 <= s, s <= pi, a >= s, a/2 >= s - pi/2]
# Confirm replacements
if not all([c.subs(rep331) for c in cond331]):
        print('rep331 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p331.subs(dict(rep331, **r1)) <= 2:
        print('Total p331 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m331)):
        if not integrate(m331[i][0], m331[i][1:]).subs(dict(rep331, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p331 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m331)):
        if not 0 <= (integrate(m331[i][0], m331[i][1:])/(m331[i][3]-m331[i][2])).subs(dict(rep331, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p331 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m331)):
        if not (m331[i][3]-m331[i][2]).subs(rep331) > 0:
                print('Bounds ' + str(i+1) + ' in p331 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p331 &= ' + latex(p331))
longLatexOutput.append('p331=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m331]).lstrip('+') + '\\right)')





##########################################################################################
# 332 animal: a <= pi.  Sensor: pi/2 <= s <= pi. Condition: a <= s and a/2 >= s- pi/2 #
##########################################################################################


m332 = [ [p2,                              g1, s/2,      pi/2           ],
          [r*sin(a/2) + r*sin(s - pi/2 - g3), g3, 0,        s - pi/2       ],
          [r*sin(a/2),                      g3, s - pi/2, s - pi/2 + a/2 ] ]

p332 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m332] ).simplify().trigsimp()


rep332 = {s:7*pi/8, a:7*pi/8} # Replacement values in range

# Define conditions for model
cond332 = [a <= pi, pi/2 <= s, s <= pi, a/2 <= s/2, a/2 >= s - pi/2]
# Confirm replacements
if not all([c.subs(rep332) for c in cond332]):
        print('rep332 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p332.subs(dict(rep332, **r1)) <= 2:
        print('Total p332 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m332)):
        if not integrate(m332[i][0], m332[i][1:]).subs(dict(rep332, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p332 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m332)):
        if not 0 <= (integrate(m332[i][0], m332[i][1:])/(m332[i][3]-m332[i][2])).subs(dict(rep332, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p332 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m332)):
        if not (m332[i][3]-m332[i][2]).subs(rep332) > 0:
                print('Bounds ' + str(i+1) + ' in p332 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p332 &= ' + latex(p332))
longLatexOutput.append('p332=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m332]).lstrip('+') + '\\right)')








##########################################################################################
# 333 animal: a <= pi.  Sensor: pi/2 <= s <= pi. Condition: a <= s and a/2 <= s- pi/2 #
##########################################################################################


m333 = [ [p2,                                g1, s/2,            pi/2           ],
          [2*r*sin(a/2),                      g3, 0,              s - pi/2 - a/2 ],
          [r*sin(a/2) + r*sin(s - pi/2 - g3), g3, s - pi/2 - a/2, s - pi/2       ],
          [r*sin(a/2),                        g3, s - pi/2,       s - pi/2 + a/2 ] ]

p333 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m333] ).simplify()


rep333 = {s:7*pi/8, a:2*pi/8} # Replacement values in range

# Define conditions for model
cond333 = [a <= pi, pi/2 <= s, s <= pi, a/2 <= s/2, a/2 <= s - pi/2]
# Confirm replacements
if not all([c.subs(rep333) for c in cond333]):
        print('rep333 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p333.subs(dict(rep333, **r1)) <= 2:
        print('Total p333 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m333)):
        if not integrate(m333[i][0], m333[i][1:]).subs(dict(rep333, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p333 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m333)):
        if not 0 <= (integrate(m333[i][0], m333[i][1:])/(m333[i][3]-m333[i][2])).subs(dict(rep333, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p333 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m333)):
        if not (m333[i][3]-m333[i][2]).subs(rep333) > 0:
                print('Bounds ' + str(i+1) + ' in p333 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p333 &= ' + latex(p333))
longLatexOutput.append('p333=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m333]).lstrip('+') + '\\right)')

###################################################################################
# 344 animal: a <= pi.  Sensor: s <= pi/2. Condition: a <= pi - 2s & a <= s     #
###################################################################################

m344 = [ [p1, g1, pi/2 - s/2 + a/2, pi/2            ],
         [p2, g1, pi/2 - s/2,       pi/2 - s/2 + a/2],
         [p3, g2, s,                s + a/2         ] ]

p344 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m344] ).simplify().trigsimp()


rep344 = {s:2*pi/8, a:pi/8} # Replacement values in range

# Define conditions for model
cond344 = [a <= pi, s <= pi/2, a <= pi - 2*s, a <= s]
# Confirm replacements
if not all([c.subs(rep344) for c in cond344]):
        print('rep344 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p344.subs(dict(rep344, **r1)) <= 2:
        print('Total p344 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m344)):
        if not integrate(m344[i][0], m344[i][1:]).subs(dict(rep344, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p344 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m344)):
        if not 0 <= (integrate(m344[i][0], m344[i][1:])/(m344[i][3]-m344[i][2])).subs(dict(rep344, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p344 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m344)):
        if not (m344[i][3]-m344[i][2]).subs(rep344) > 0:
                print('Bounds ' + str(i+1) + ' in p344 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p344 &= ' + latex(p344))
longLatexOutput.append('p344=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m344]).lstrip('+') + '\\right)')


#######################################################################################
# 345 animal: a <= pi.  Sensor: s <= pi/2. Condition: a <= pi - 2s & s <= a <= 2s   #
#######################################################################################


m345 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 + s/2 - a/2, pi/2            ],
         [p2,                   g1, pi/2 - s/2,       pi/2 + s/2 - a/2],
         [p3,                   g2, s,                s + a/2         ] ]

p345 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m345] ).simplify().trigsimp()



rep345 = {s:2*pi/8, a:pi/2-0.1} # Replacement values in range

# Define conditions for model
cond345 = [a <= pi, s <= pi/2, a <= pi - 2*s, s <= a, a <= 2*s]
# Confirm replacements
if not all([c.subs(rep345) for c in cond345]):
        print('rep345 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p345.subs(dict(rep345, **r1)) <= 2:
        print('Total p345 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m345)):
        if not integrate(m345[i][0], m345[i][1:]).subs(dict(rep345, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p345 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m345)):
        if not 0 <= (integrate(m345[i][0], m345[i][1:])/(m345[i][3]-m345[i][2])).subs(dict(rep345, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p345 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m345)):
        if not (m345[i][3]-m345[i][2]).subs(rep345) > 0:
                print('Bounds ' + str(i+1) + ' in p345 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p345 &= ' + latex(p345))
longLatexOutput.append('p345=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m345]).lstrip('+') + '\\right)')


##################################################################################
# 346 animal: a <= pi.  Sensor: s <= pi/2. Condition: a <= pi - 2s &  2s <= a      #
##################################################################################

# This DOES = 343, just can't show it.

m346 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2    ],
         [r*sin(g2),            g2, s,          a/2     ],
         [p3,                   g2, a/2,        s + a/2 ] ]

p346 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m346] ).simplify().trigsimp()


rep346 = {s:1*pi/8, a:pi/2} # Replacement values in range

# Define conditions for model
cond346 = [a <= pi,  s <= pi/2,  a <= pi - 2*s,  2*s <= a]
# Confirm replacements
if not all([c.subs(rep346) for c in cond346]):
        print('rep346 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p346.subs(dict(rep346, **r1)) <= 2:
        print('Total p346 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m346)):
        if not integrate(m346[i][0], m346[i][1:]).subs(dict(rep346, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p346 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m346)):
        if not 0 <= (integrate(m346[i][0], m346[i][1:])/(m346[i][3]-m346[i][2])).subs(dict(rep346, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p346 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m346)):
        if not (m346[i][3]-m346[i][2]).subs(rep346) > 0:
                print('Bounds ' + str(i+1) + ' in p346 has lower bounds bigger than upper bounds')        


pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p221Profile.pdf')
pl.close()

# LaTeX output

latexOutput.append('p346 &= ' + latex(p346))
longLatexOutput.append('p346=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m346]).lstrip('+') + '\\right)')

##################################################################################
# 341 animal: a <= pi.  Sensor: s <= pi/2. Condition: a > pi - 2s &  a <= s       #
##################################################################################


m341 = [ [p1,         g1, pi/2 - s/2 + a/2, pi/2            ],
         [p2,         g1, pi/2 - s/2,       pi/2 - s/2 + a/2],
         [p3,         g2, s,                pi/2            ],
         [r*sin(a/2), g3, 0,                a/2 + s - pi/2  ] ]

p341 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m341] ).simplify().trigsimp()

rep341 = {s:pi/2-0.1, a:pi/4} # Replacement values in range

# Define conditions for model
cond341 = [a <= pi,  s <= pi/2,  a >= pi - 2*s,  a <= s]
# Confirm replacements
if not all([c.subs(rep341) for c in cond341]):
        print('rep341 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p341.subs(dict(rep341, **r1)) <= 2:
        print('Total p341 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m341)):
        if not integrate(m341[i][0], m341[i][1:]).subs(dict(rep341, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p341 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m341)):
        if not 0 <= (integrate(m341[i][0], m341[i][1:])/(m341[i][3]-m341[i][2])).subs(dict(rep341, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p341 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m341)):
        if not (m341[i][3]-m341[i][2]).subs(rep341) > 0:
                print('Bounds ' + str(i+1) + ' in p341 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p341 &= ' + latex(p341))
longLatexOutput.append('p341=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m341]).lstrip('+') + '\\right)')

######################################################################################
# 342 animal: a <= pi.  Sensor: s <= pi/2. Condition: a > pi - 2s &  s <= a <= 2s  #
######################################################################################


m342 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 + s/2 - a/2, pi/2            ],
         [p2,                   g1, pi/2 - s/2,       pi/2 + s/2 - a/2],
         [p3,                   g2, s,                pi/2        ],
         [r*sin(a/2),           g3, 0,                a/2 + s -pi/2   ] ]

p342 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m342] ).simplify().trigsimp()



rep342 = {s:pi/2-0.1, a:pi/2} # Replacement values in range

# Confirm replacements
if not (a.subs(rep342) <= pi and s.subs(rep342) <= pi/2 and a.subs(rep342) >= pi - 2*s.subs(rep342) and s.subs(rep342) <= a.subs(rep342) <= 2*s.subs(rep342)):
        print('rep342 incorrect')

# Define conditions for model

cond342 = [a <= pi,  s <= pi/2,  a >= pi - 2*s,  s <= a, a <= 2*s]
# Confirm replacements
if not all([c.subs(rep342) for c in cond342]):
        print('rep342 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p342.subs(dict(rep342, **r1)) <= 2:
        print('Total p342 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m342)):
        if not integrate(m342[i][0], m342[i][1:]).subs(dict(rep342, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p342 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m342)):
        if not 0 <= (integrate(m342[i][0], m342[i][1:])/(m342[i][3]-m342[i][2])).subs(dict(rep342, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p342 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m342)):
        if not (m342[i][3]-m342[i][2]).subs(rep342) > 0:
                print('Bounds ' + str(i+1) + ' in p342 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p342 &= ' + latex(p342))
longLatexOutput.append('p342=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m342]).lstrip('+') + '\\right)')





##################################################################################
# 343 animal: a <= pi.  Sensor: s <= pi/2. Condition: a > pi - 2s &  a > 2s      #
##################################################################################



m343 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2            ],
         [r*sin(g2),            g2, s,          a/2             ],
         [p3,                   g2, a/2,        pi/2            ],
         [r*sin(a/2),           g3, 0,          a/2 + s -pi/2   ] ]

p343 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m343] ).simplify().trigsimp()



rep343 = {s:pi/4, a:3*pi/4} # Replacement values in range


# Define conditions for model
cond343 = [a <= pi,  s <= pi/2,  a >= pi - 2*s,  a > 2*s]
# Confirm replacements
if not all([c.subs(rep343) for c in cond343]):
        print('rep343 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p343.subs(dict(rep343, **r1)) <= 2:
        print('Total p343 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m343)):
        if not integrate(m343[i][0], m343[i][1:]).subs(dict(rep343, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p343 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m343)):
        if not 0 <= (integrate(m343[i][0], m343[i][1:])/(m343[i][3]-m343[i][2])).subs(dict(rep343, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p343 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m343)):
        if not (m343[i][3]-m343[i][2]).subs(rep343) > 0:
                print('Bounds ' + str(i+1) + ' in p343 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p343 &= ' + latex(p343))
longLatexOutput.append('p343=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m343]).lstrip('+') + '\\right)')



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


# I'm pretty sure 346 and 343 are currently equal, even though this can't find that. 
# p344 not right for values close to zero?

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
['p222', 'p321',{s:pi}],

['p223', 'p222',{a:4*pi-2*s}],
['p223', 'p222',{s:2*pi-a/2}],
['p223', 'p232',{a:pi}],
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
['p141','p241',{s:2*pi}],

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

['p321','p322',{s:a-pi}],
['p321','p322',{a:4*pi-2*s}],
['p321','p311',{s:2*pi-a/2}],
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
['p333','p326',{s:pi}],


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

# Run through all the comparisons. Need simplify(). Even together() gives some false negatives.

checkFile = open('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/checksFile.tex','w')

checkFile.write('All checks evaluated.\nTim Lucas - ' + str(datetime.now()) + '\n')
checkFile.write('I\'m pretty sure 346 and 343 are currently equal, even though this can\'t find that.')
for i in range(len(allComps)):
        if (eval(allComps[i][0]).subs(allComps[i][2]) - eval(allComps[i][1]).subs(allComps[i][2])).simplify() == 0:
                checkFile.write(str(i) + ': ' + allComps[i][0]+ ' and ' +allComps[i][1]+': OK\n')
        else:
                checkFile.write(str(i) + ': ' + allComps[i][0]+ ' and ' +allComps[i][1]+': Incorrect\n')

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

latexFile = open('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/ModelSolutions.tex', 'w')

latexFile.write('% LaTeX output. Solutions of all REM models.\n' + '%Tim Lucas - ' + str(datetime.now()) + '\n')
for i in range(len(latexOutput)):
        latexFile.write( '\\[' + latexOutput[i] + '\\]\n')

latexFile.close()

latexFile = open('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/ModelDefinitions.tex', 'w')

latexFile.write('% LaTeX output. Definitions of all REM models.\n' + '%Tim Lucas - ' + str(datetime.now()) + '\n')
for i in range(len(longLatexOutput)):
        latexFile.write( '\\[' + longLatexOutput[i] + '\\]\n')

latexFile.close()





