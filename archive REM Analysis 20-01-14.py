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

###################################################
# 7.3.1 animal: a = 2*pi.  sensor: s > pi         #
###################################################



m731 = [ [2*r,                 g4, pi/2, s/2        ],
         [r + r*cos(g4 - s/2), g4, s/2,  pi         ],
         [r + r*cos(g4 + s/2), g4, pi,   2*pi-s/2   ],
         [2*r,                 g4, 2*pi-s/2, 3*pi/2 ] ]

p731 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m731] ).simplify().trigsimp()




# Replacement values in range
rep731 = {s:3*pi/2, a:2*pi} 

# Define conditions for model
cond731 = [pi <= s]
# Confirm replacements
if not all([c.subs(rep731) for c in cond731]):
        print('rep731 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p731.subs(dict(rep731, **r1)) <= 2:
        print('Total p731 not in 0, 2r')

# Are the individuals integrals >0r
for i in range(len(m731)):
        if not integrate(m731[i][0], m731[i][1:]).subs(dict(rep731, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p731 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m731)):
        if not 0 <= (integrate(m731[i][0], m731[i][1:])/(m731[i][3]-m731[i][2])).subs(dict(rep731, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p731 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m731)):
        if not (m731[i][3]-m731[i][2]).subs(rep731) > 0:
                print('Bounds ' + str(i+1) + ' in p731 has lower bounds bigger than upper bounds')
        
# Plot function

xRange = np.arange(pi,2*pi, 0.01)
yRange = [p731.subs({r:1, s:i}).n() for i in xRange]
#plot731 = pl.plot(xRange, yRange)
#pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p731Profile.pdf')
#pl.close()

# LaTeX output

latexOutput.append('p731 &= ' + latex(p731))

longLatexOutput.insert(0,'p731=&\\frac{1}{\pi} \left(2\int_{'+latex(m731[0][2])+'}^{'+latex(m731[0][3])+'}'+latex(m731[0][0])+'\;\mathrm{d}'+latex(m731[0][1])+'+2\int_{'+latex(m731[1][2])+'}^{'+latex(m731[1][3])+'}'+latex(m731[1][0])+'\;\mathrm{d}'+latex(m731[1][1])+'\\right)')



##################################################################
# 7.3.2 animal: a > pi.  sensor: s > pi Condition: a < 5pi - 2s  #
##################################################################



m732 = [ [2*r,                 g4, pi/2, s/2        ],
         [r + r*cos(g4 - s/2), g4, s/2,  pi         ],
         [r + r*cos(g4 - s/2), g4, pi,   5*pi/2 - s/2 - a/2 ],
         [r + r*cos(g4 + s/2), g4, 5*pi/2 - s/2 - a/2,   2*pi-s/2 ],
         [2*r,                 g4, 2*pi-s/2, 3*pi/2 ] ]

p732 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m732] ).simplify().trigsimp()

# Replacement values in range
rep732 = {s:5*pi/4-0.1, a:3*pi/2} 

# Define conditions for model
cond732 = [pi <= s, a >= pi, a <= 5*pi - 2*s]
## Confirm replacementsa
if not all([c.subs(rep732) for c in cond732]):
        print('rep732 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p732.subs(dict(rep732, **r1)) <= 2:
        print('Total p732 not in 0, 2r')

# Are the individuals integrals >0r
for i in range(len(m732)):
        if not integrate(m732[i][0], m732[i][1:]).subs(dict(rep732, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p732 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m732)):
        if not 0 <= (integrate(m732[i][0], m732[i][1:])/(m732[i][3]-m732[i][2])).subs(dict(rep732, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p732 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m732)):
        if not (m732[i][3]-m732[i][2]).subs(rep732) > 0:
                print('Bounds ' + str(i+1) + ' in p732 has lower bounds bigger than upper bounds')
        
# Plot function

xRange = np.arange(pi,2*pi, 0.01)
yRange = [p732.subs({r:1, s:i}).n() for i in xRange]
#plot732 = pl.plot(xRange, yRange)
#pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p732Profile.pdf')
#pl.close()

# LaTeX output

latexOutput.append('p732 &= ' + latex(p732))
longLatexOutput.append('p732=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m732]).lstrip('+') + '\\right)')



########################################################
# 7.5. animal: a = 2*pi.   sensor:  pi/2 <= s <= pi      #
########################################################

m75 = [ [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2],
        [r - r*cos(g3 - s),    g3, 0, s - pi/2],
        [r,                    g3, s - pi/2, pi/2],
        [r - r*cos(g3),    g3, pi/2, s],
        [2*r*sin(s/2)*sin(g1),  g1, s/2, pi/2] ]

p75 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m75] ).simplify().trigsimp()


# Replacement values in range
rep75 = {s:3*pi/4} 

# Define conditions for model
cond75 = [pi/2 <= s, s <= pi]
# Confirm replacements
if not all([c.subs(rep75) for c in cond75]):
        print('rep75 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p75.subs(dict(rep75, **r1)) <= 2:
        print('Total p75 not in 0, 2r')

# Are the individuals integrals >0r
for i in range(len(m75)):
        if not integrate(m75[i][0], m75[i][1:]).subs(dict(rep75, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p75 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m75)):
        if not 0 <= (integrate(m75[i][0], m75[i][1:])/(m75[i][3]-m75[i][2])).subs(dict(rep75, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p75 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m75)):
        if not (m75[i][3]-m75[i][2]).subs(rep75) > 0:
                print('Bounds ' + str(i+1) + ' in p75 has lower bounds bigger than upper bounds')        


# Plot function

xRange = np.arange(pi/2,pi, 0.01)
yRange = [p75.subs({r:1, s:i}).n() for i in xRange]
#plot75 = pl.plot(xRange, yRange)
#pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p75Profile.pdf')
#pl.close()

# LaTeX output

latexOutput.append('p75 &= ' + latex(p75))
longLatexOutput.append('p75=&\\frac{1}{\pi} \left(\int_{'+latex(m75[0][2])+'}^{'+latex(m75[0][3])+'}'+latex(m75[0][0])+'\;\mathrm{d}'+latex(m75[0][1])+'+\int_{'+latex(m75[1][2])+'}^{'+latex(m75[1][3])+'}'+latex(m75[1][0])+'\;\mathrm{d}'+latex(m75[1][1])+'+\int_{'+latex(m75[2][2])+'}^{'+latex(m75[2][3])+'}'+latex(m75[2][0])+'\;\mathrm{d}'+latex(m75[2][1])+'\\right)')




#####################################################################
# 7.7.1 animal: a>pi.  Sensor: s <= pi/2. Condition: a <= 2pi - 2s  #
#####################################################################

m771 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
          [r*sin(g2),            g2, s,          pi/2],
          [r,                    g3, 0,          s],
          [r*sin(g2),            g2, pi - a/2,   pi/2] ]


p771 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m771] ).simplify().trigsimp()


rep771 = {s:pi/9, a:10*pi/9} # Replacement values in range

# Define conditions for model
cond771 = [s <= pi/2, a >= pi, a <= 2*pi - 2*s]
# Confirm replacements
if not all([c.subs(rep771) for c in cond771]):
        print('rep771 incorrect')
# is average profile in range 0r-2r?
if not 0 <= p771.subs(dict(rep771, **r1)) <= 2:
        print('Total p771 not in 0, 2r')

# Are the individuals integrals >0r
for i in range(len(m771)):
        if not integrate(m771[i][0], m771[i][1:]).subs(dict(rep771, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p771 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m771)):
        if not 0 <= (integrate(m771[i][0], m771[i][1:])/(m771[i][3]-m771[i][2])).subs(dict(rep771, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p771 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m771)):
        if not (m771[i][3]-m771[i][2]).subs(rep771) > 0:
                print('Bounds ' + str(i+1) + ' in p771 has lower bounds bigger than upper bounds')        


# LaTeX output

latexOutput.append('p771 &= ' + latex(p771))
longLatexOutput.append('p771=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m771]).lstrip('+') + '\\right)')


####################################################################################
# 7.7.2 animal: a>pi.  Sensor: s <= pi/2. Condition:  2*pi - 2*s <= a <= 2*pi - s  #
####################################################################################

m772 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s,          pi/2],
         [r,                    g3, 0,          s],
         [r*sin(g2),            g2, s,          pi/2],
         [r*cos(g1 - s/2),      g1, pi/2 - s/2, a/2 + s/2 - pi/2] ]


p772 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m772] ).simplify().trigsimp()


rep772 = {s:3*pi/8, a:3*pi/2} # Replacement values in range

# Define conditions for model
cond772 = [a >= pi, s <= pi/2, 2*pi - 2*s <= a, a <= 2*pi - s]
# Confirm replacements
if not all([c.subs(rep772) for c in cond772]):
        print('rep772 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p772.subs(dict(rep772, **r1)) <= 2:
        print('Total p772 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m772)):
        if not integrate(m772[i][0], m772[i][1:]).subs(dict(rep772, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p772 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m772)):
        if not 0 <= (integrate(m772[i][0], m772[i][1:])/(m772[i][3]-m772[i][2])).subs(dict(rep772, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p772 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m772)):
        if not (m772[i][3]-m772[i][2]).subs(rep772) > 0:
                print('Bounds ' + str(i+1) + ' in p772 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p772 &= ' + latex(p772))
longLatexOutput.append('p772 &= \\frac{1}{\pi} \left(\int_{'+latex(m772[0][2])+'}^{'+latex(m772[0][3])+'}'+latex(m772[0][0])+'\;\mathrm{d}'+latex(m772[0][1])+'+\int_{'+latex(m772[1][2])+'}^{'+latex(m772[1][3])+'}'+latex(m772[1][0])+'\;\mathrm{d}'+latex(m772[1][1])+'+\int_{'+latex(m772[2][2])+'}^{'+latex(m772[2][3])+'}'+latex(m772[2][0])+'\;\mathrm{d}'+latex(m772[2][1])+'+\int_{'+latex(m772[3][2])+'}^{'+latex(m772[3][3])+'}'+latex(m772[3][0])+'\;\mathrm{d}'+latex(m772[3][1])+'\\right)')


###############################################################################
# 7.7.3 animal: a>pi.  Sensor: s <= pi/2. Condition: 2*pi - s < a             #
###############################################################################

m773 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s,          pi/2],
         [r,                    g3, 0,          s],
         [r*sin(g2),            g2, s,          pi/2],
         [r*cos(g1 - s/2),      g1, pi/2 - s/2, 3*pi/2 - s/2 - a/2],
         [2*r*sin(s/2)*sin(g1), g1, 3*pi/2 - s/2 - a/2, pi/2] ]

p773 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m773] ).simplify().trigsimp()


rep773 = {s:3*pi/8, a:29*pi/16} # Replacement values in range

# Define conditions for model
cond773 = [a >= pi, s <= pi/2, 2*pi - s <= a  ]
# Confirm replacements
if not all([c.subs(rep773) for c in cond773]):
        print('rep773 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p773.subs(dict(rep773, **r1)) <= 2:
        print('Total p773 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m773)):
        if not integrate(m773[i][0], m773[i][1:]).subs(dict(rep773, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p773 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m773)):
        if not 0 <= (integrate(m773[i][0], m773[i][1:])/(m773[i][3]-m773[i][2])).subs(dict(rep773, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p773 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m773)):
        if not (m773[i][3]-m773[i][2]).subs(rep773) > 0:
                print('Bounds ' + str(i+1) + ' in p773 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p773 &= ' + latex(p773))
longLatexOutput.append('p773 &= \\frac{1}{\pi} \left(\int_{'+latex(m773[0][2])+'}^{'+latex(m773[0][3])+'}'+latex(m773[0][0])+'\;\mathrm{d}'+latex(m773[0][1])+'+\int_{'+latex(m773[1][2])+'}^{'+latex(m773[1][3])+'}'+latex(m773[1][0])+'\;\mathrm{d}'+latex(m773[1][1])+'+\int_{'+latex(m773[2][2])+'}^{'+latex(m773[2][3])+'}'+latex(m773[2][0])+'\;\mathrm{d}'+latex(m773[2][1])+'+\int_{'+latex(m773[3][2])+'}^{'+latex(m773[3][3])+'}'+latex(m773[3][0])+'\;\mathrm{d}'+latex(m773[3][1])+'\\right)')





#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
#Two models in 7.8.
#animal: pi <= a .  Sensor: pi/2 <= s <= pi
#"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#################################################################################
# 7.8.1 animal:  a > pi.  Sensor: pi/2 <= s <= pi. Condition: a <= 2pi - s      #
#################################################################################

m781 = [ [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2],
         [r - r*cos(g3-s),      g3, 0, s - pi/2],
         [r,                    g3, s - pi/2, pi/2],
         [r,                    g3, pi/2, s],
         [r*cos(g1 - s/2),      g1, s/2, a/2 + s/2 - pi/2],
         [0,                    g1, a/2 + s/2 - pi/2, pi/2 ] ]

p781 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m781] ).simplify().trigsimp()

rep781 = {s:3*pi/4, a:9*pi/8} # Replacement values in range

# Define conditions for model
cond781 = [a > pi,  pi/2 <= s, s <= pi, a <= 2*pi - s]
# Confirm replacements
if not all([c.subs(rep781) for c in cond781]):
        print('rep781 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p781.subs(dict(rep781, **r1)) <= 2:
        print('Total p781 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m781)):
        if not integrate(m781[i][0], m781[i][1:]).subs(dict(rep781, **r1)) >= 0:
                print('Integral ' + str(i+1) + ' in p781 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m781)):
        if not 0 <= (integrate(m781[i][0], m781[i][1:])/(m781[i][3]-m781[i][2])).subs(dict(rep781, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p781 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m781)):
        if not (m781[i][3]-m781[i][2]).subs(rep781) >= 0:
                print('Bounds ' + str(i+1) + ' in p781 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p781 &= ' + latex(p781))
longLatexOutput.append('p781=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m781]).lstrip('+') + '\\right)')



#################################################################################
# 7.8.2 animal: a > pi.  Sensor: pi/2 <= s <= pi. Cond: 2pi - s < a < 3pi - 2s  #
#################################################################################


m782 = [ [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2],
         [r - r*cos(g3-s),      g3, 0, s - pi/2],
         [r,                    g3, s - pi/2, pi/2],
         [r,                    g3, pi/2, s],
         [r*cos(g1 - s/2),      g1, s/2, 3*pi/2 - a/2 - s/2],
         [2*r*sin(s/2)*sin(g1), g1, 3*pi/2 - a/2 - s/2, pi/2 ] ]


p782 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m782] ).simplify().trigsimp()



rep782 = {s:5*pi/8, a:6*pi/4} # Replacement values in range

# Define conditions for model
cond782 = [a > pi, pi/2 <= s, s <= pi, 2*pi - s <= a, a <= 3*pi - 2*s]
# Confirm replacements
if not all([c.subs(rep782) for c in cond782]):
        print('rep782 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p782.subs(dict(rep782, **r1)) <= 2:
        print('Total p782 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m782)):
        if not integrate(m782[i][0], m782[i][1:]).subs(dict(rep782, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p782 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m782)):
        if not 0 <= (integrate(m782[i][0], m782[i][1:])/(m782[i][3]-m782[i][2])).subs(dict(rep782, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p782 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m782)):
        if not (m782[i][3]-m782[i][2]).subs(rep782) > 0:
                print('Bounds ' + str(i+1) + ' in p782 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p782 &= ' + latex(p782))
longLatexOutput.append('p782=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m782]).lstrip('+') + '\\right)')




#################################################################################
# 7.8.3 animal: a > pi.  Sensor: pi/2 <= s <= pi. Condition: a/2 > pi - s/2  #
#################################################################################


m783 = [ [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2],
         [r - r*cos(g3-s),      g3, 0, s - pi/2],
         [r,                    g3, s - pi/2, pi/2],
         [r,                    g3, pi/2, 3*pi/2 - a/2],
         [r-r*cos(g3),          g3, 3*pi/2 - a/2, s],
         [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2 ] ]


p783 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m783] ).simplify().trigsimp()



rep783 = {s:3*pi/4, a:15*pi/8} # Replacement values in range

# Define conditions for model
cond783 = [a > pi, pi/2 <= s, s <= pi, a >= 3*pi - 2*s]

# Confirm replacements
if not all([c.subs(rep783) for c in cond783]):
        print('rep783 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p783.subs(dict(rep783, **r1)) <= 2:
        print('Total p783 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m783)):
        if not integrate(m783[i][0], m783[i][1:]).subs(dict(rep783, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p783 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m783)):
        if not 0 <= (integrate(m783[i][0], m783[i][1:])/(m783[i][3]-m783[i][2])).subs(dict(rep783, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p783 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m783)):
        if not (m783[i][3]-m783[i][2]).subs(rep783) > 0:
                print('Bounds ' + str(i+1) + ' in p783 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p783 &= ' + latex(p783))
longLatexOutput.append('p783=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m783]).lstrip('+') + '\\right)')






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
# 7.10.1 animal: a <= pi.  Sensor: pi/2 <= s <= pi. Condition: a >= s and a/2 >= s - pi/2 #
###########################################################################################


m7101 = [ [2*r*sin(s/2)*sin(g1),              g1, pi/2 - a/2 + s/2, pi/2            ],
          [p2,                                g1, s/2,              pi/2 - a/2 + s/2],
          [r*sin(a/2) + r*sin(s - pi/2 - g3), g3, 0,                s - pi/2        ],
          [r*sin(a/2),                        g3, s-pi/2,           s - pi/2 + a/2  ] ]

p7101 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m7101] ).simplify().trigsimp()



rep7101 = {s:5*pi/8, a:6*pi/8} # Replacement values in range

# Define conditions for model
cond7101 = [a <= pi, pi/2 <= s, s <= pi, a >= s, a/2 >= s - pi/2]
# Confirm replacements
if not all([c.subs(rep7101) for c in cond7101]):
        print('rep7101 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p7101.subs(dict(rep7101, **r1)) <= 2:
        print('Total p7101 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m7101)):
        if not integrate(m7101[i][0], m7101[i][1:]).subs(dict(rep7101, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p7101 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m7101)):
        if not 0 <= (integrate(m7101[i][0], m7101[i][1:])/(m7101[i][3]-m7101[i][2])).subs(dict(rep7101, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p7101 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m7101)):
        if not (m7101[i][3]-m7101[i][2]).subs(rep7101) > 0:
                print('Bounds ' + str(i+1) + ' in p7101 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p7101 &= ' + latex(p7101))
longLatexOutput.append('p7101=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m7101]).lstrip('+') + '\\right)')





##########################################################################################
# 7.10.2 animal: a <= pi.  Sensor: pi/2 <= s <= pi. Condition: a <= s and a/2 >= s- pi/2 #
##########################################################################################


m7102 = [ [p2,                              g1, s/2,      pi/2           ],
          [r*sin(a/2) + r*sin(s - pi/2 - g3), g3, 0,        s - pi/2       ],
          [r*sin(a/2),                      g3, s - pi/2, s - pi/2 + a/2 ] ]

p7102 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m7102] ).simplify().trigsimp()


rep7102 = {s:7*pi/8, a:7*pi/8} # Replacement values in range

# Define conditions for model
cond7102 = [a <= pi, pi/2 <= s, s <= pi, a/2 <= s/2, a/2 >= s - pi/2]
# Confirm replacements
if not all([c.subs(rep7102) for c in cond7102]):
        print('rep7102 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p7102.subs(dict(rep7102, **r1)) <= 2:
        print('Total p7102 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m7102)):
        if not integrate(m7102[i][0], m7102[i][1:]).subs(dict(rep7102, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p7102 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m7102)):
        if not 0 <= (integrate(m7102[i][0], m7102[i][1:])/(m7102[i][3]-m7102[i][2])).subs(dict(rep7102, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p7102 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m7102)):
        if not (m7102[i][3]-m7102[i][2]).subs(rep7102) > 0:
                print('Bounds ' + str(i+1) + ' in p7102 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p7102 &= ' + latex(p7102))
longLatexOutput.append('p7102=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m7102]).lstrip('+') + '\\right)')








##########################################################################################
# 7.10.3 animal: a <= pi.  Sensor: pi/2 <= s <= pi. Condition: a <= s and a/2 <= s- pi/2 #
##########################################################################################


m7103 = [ [p2,                                g1, s/2,            pi/2           ],
          [2*r*sin(a/2),                      g3, 0,              s - pi/2 - a/2 ],
          [r*sin(a/2) + r*sin(s - pi/2 - g3), g3, s - pi/2 - a/2, s - pi/2       ],
          [r*sin(a/2),                        g3, s - pi/2,       s - pi/2 + a/2 ] ]

p7103 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m7103] ).simplify()


rep7103 = {s:7*pi/8, a:2*pi/8} # Replacement values in range

# Define conditions for model
cond7103 = [a <= pi, pi/2 <= s, s <= pi, a/2 <= s/2, a/2 <= s - pi/2]
# Confirm replacements
if not all([c.subs(rep7103) for c in cond7103]):
        print('rep7103 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p7103.subs(dict(rep7103, **r1)) <= 2:
        print('Total p7103 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m7103)):
        if not integrate(m7103[i][0], m7103[i][1:]).subs(dict(rep7103, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p7103 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m7103)):
        if not 0 <= (integrate(m7103[i][0], m7103[i][1:])/(m7103[i][3]-m7103[i][2])).subs(dict(rep7103, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p7103 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m7103)):
        if not (m7103[i][3]-m7103[i][2]).subs(rep7103) > 0:
                print('Bounds ' + str(i+1) + ' in p7103 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p7103 &= ' + latex(p7103))
longLatexOutput.append('p7103=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m7103]).lstrip('+') + '\\right)')

###################################################################################
# 7.9.1 animal: a <= pi.  Sensor: s <= pi/2. Condition: a <= pi - 2s & a <= s     #
###################################################################################

m791 = [ [p1, g1, pi/2 - s/2 + a/2, pi/2            ],
         [p2, g1, pi/2 - s/2,       pi/2 - s/2 + a/2],
         [p3, g2, s,                s + a/2         ] ]

p791 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m791] ).simplify().trigsimp()


rep791 = {s:2*pi/8, a:pi/8} # Replacement values in range

# Define conditions for model
cond791 = [a <= pi, s <= pi/2, a <= pi - 2*s, a <= s]
# Confirm replacements
if not all([c.subs(rep791) for c in cond791]):
        print('rep791 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p791.subs(dict(rep791, **r1)) <= 2:
        print('Total p791 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m791)):
        if not integrate(m791[i][0], m791[i][1:]).subs(dict(rep791, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p791 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m791)):
        if not 0 <= (integrate(m791[i][0], m791[i][1:])/(m791[i][3]-m791[i][2])).subs(dict(rep791, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p791 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m791)):
        if not (m791[i][3]-m791[i][2]).subs(rep791) > 0:
                print('Bounds ' + str(i+1) + ' in p791 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p791 &= ' + latex(p791))
longLatexOutput.append('p791=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m791]).lstrip('+') + '\\right)')


#######################################################################################
# 7.9.2 animal: a <= pi.  Sensor: s <= pi/2. Condition: a <= pi - 2s & s <= a <= 2s   #
#######################################################################################


m792 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 + s/2 - a/2, pi/2            ],
         [p2,                   g1, pi/2 - s/2,       pi/2 + s/2 - a/2],
         [p3,                   g2, s,                s + a/2         ] ]

p792 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m792] ).simplify().trigsimp()



rep792 = {s:2*pi/8, a:pi/2-0.1} # Replacement values in range

# Define conditions for model
cond792 = [a <= pi, s <= pi/2, a <= pi - 2*s, s <= a, a <= 2*s]
# Confirm replacements
if not all([c.subs(rep792) for c in cond792]):
        print('rep792 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p792.subs(dict(rep792, **r1)) <= 2:
        print('Total p792 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m792)):
        if not integrate(m792[i][0], m792[i][1:]).subs(dict(rep792, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p792 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m792)):
        if not 0 <= (integrate(m792[i][0], m792[i][1:])/(m792[i][3]-m792[i][2])).subs(dict(rep792, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p792 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m792)):
        if not (m792[i][3]-m792[i][2]).subs(rep792) > 0:
                print('Bounds ' + str(i+1) + ' in p792 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p792 &= ' + latex(p792))
longLatexOutput.append('p792=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m792]).lstrip('+') + '\\right)')


##################################################################################
# 7.9.3 animal: a <= pi.  Sensor: s <= pi/2. Condition: a <= pi - 2s &  2s <= a      #
##################################################################################

# This DOES = 796, just can't show it.

m793 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2    ],
         [r*sin(g2),            g2, s,          a/2     ],
         [p3,                   g2, a/2,        s + a/2 ] ]

p793 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m793] ).simplify().trigsimp()


rep793 = {s:1*pi/8, a:pi/2} # Replacement values in range

# Define conditions for model
cond793 = [a <= pi,  s <= pi/2,  a <= pi - 2*s,  2*s <= a]
# Confirm replacements
if not all([c.subs(rep793) for c in cond793]):
        print('rep793 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p793.subs(dict(rep793, **r1)) <= 2:
        print('Total p793 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m793)):
        if not integrate(m793[i][0], m793[i][1:]).subs(dict(rep793, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p793 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m793)):
        if not 0 <= (integrate(m793[i][0], m793[i][1:])/(m793[i][3]-m793[i][2])).subs(dict(rep793, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p793 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m793)):
        if not (m793[i][3]-m793[i][2]).subs(rep793) > 0:
                print('Bounds ' + str(i+1) + ' in p793 has lower bounds bigger than upper bounds')        


pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p731Profile.pdf')
pl.close()

# LaTeX output

latexOutput.append('p793 &= ' + latex(p793))
longLatexOutput.append('p793=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m793]).lstrip('+') + '\\right)')

##################################################################################
# 7.9.4 animal: a <= pi.  Sensor: s <= pi/2. Condition: a > pi - 2s &  a <= s       #
##################################################################################


m794 = [ [p1,         g1, pi/2 - s/2 + a/2, pi/2            ],
         [p2,         g1, pi/2 - s/2,       pi/2 - s/2 + a/2],
         [p3,         g2, s,                pi/2            ],
         [r*sin(a/2), g3, 0,                a/2 + s - pi/2  ] ]

p794 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m794] ).simplify().trigsimp()

rep794 = {s:pi/2-0.1, a:pi/4} # Replacement values in range

# Define conditions for model
cond794 = [a <= pi,  s <= pi/2,  a >= pi - 2*s,  a <= s]
# Confirm replacements
if not all([c.subs(rep794) for c in cond794]):
        print('rep794 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p794.subs(dict(rep794, **r1)) <= 2:
        print('Total p794 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m794)):
        if not integrate(m794[i][0], m794[i][1:]).subs(dict(rep794, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p794 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m794)):
        if not 0 <= (integrate(m794[i][0], m794[i][1:])/(m794[i][3]-m794[i][2])).subs(dict(rep794, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p794 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m794)):
        if not (m794[i][3]-m794[i][2]).subs(rep794) > 0:
                print('Bounds ' + str(i+1) + ' in p794 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p794 &= ' + latex(p794))
longLatexOutput.append('p794=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m794]).lstrip('+') + '\\right)')

######################################################################################
# 7.9.5 animal: a <= pi.  Sensor: s <= pi/2. Condition: a > pi - 2s &  s <= a <= 2s  #
######################################################################################


m795 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 + s/2 - a/2, pi/2            ],
         [p2,                   g1, pi/2 - s/2,       pi/2 + s/2 - a/2],
         [p3,                   g2, s,                pi/2        ],
         [r*sin(a/2),           g3, 0,                a/2 + s -pi/2   ] ]

p795 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m795] ).simplify().trigsimp()



rep795 = {s:pi/2-0.1, a:pi/2} # Replacement values in range

# Confirm replacements
if not (a.subs(rep795) <= pi and s.subs(rep795) <= pi/2 and a.subs(rep795) >= pi - 2*s.subs(rep795) and s.subs(rep795) <= a.subs(rep795) <= 2*s.subs(rep795)):
        print('rep795 incorrect')

# Define conditions for model

cond795 = [a <= pi,  s <= pi/2,  a >= pi - 2*s,  s <= a, a <= 2*s]
# Confirm replacements
if not all([c.subs(rep795) for c in cond795]):
        print('rep795 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p795.subs(dict(rep795, **r1)) <= 2:
        print('Total p795 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m795)):
        if not integrate(m795[i][0], m795[i][1:]).subs(dict(rep795, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p795 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m795)):
        if not 0 <= (integrate(m795[i][0], m795[i][1:])/(m795[i][3]-m795[i][2])).subs(dict(rep795, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p795 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m795)):
        if not (m795[i][3]-m795[i][2]).subs(rep795) > 0:
                print('Bounds ' + str(i+1) + ' in p795 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p795 &= ' + latex(p795))
longLatexOutput.append('p795=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m795]).lstrip('+') + '\\right)')





##################################################################################
# 7.9.6 animal: a <= pi.  Sensor: s <= pi/2. Condition: a > pi - 2s &  a > 2s      #
##################################################################################



m796 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2            ],
         [r*sin(g2),            g2, s,          a/2             ],
         [p3,                   g2, a/2,        pi/2            ],
         [r*sin(a/2),           g3, 0,          a/2 + s -pi/2   ] ]

p796 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m796] ).simplify().trigsimp()



rep796 = {s:pi/4, a:3*pi/4} # Replacement values in range


# Define conditions for model
cond796 = [a <= pi,  s <= pi/2,  a >= pi - 2*s,  a > 2*s]
# Confirm replacements
if not all([c.subs(rep796) for c in cond796]):
        print('rep796 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p796.subs(dict(rep796, **r1)) <= 2:
        print('Total p796 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m796)):
        if not integrate(m796[i][0], m796[i][1:]).subs(dict(rep796, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p796 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m796)):
        if not 0 <= (integrate(m796[i][0], m796[i][1:])/(m796[i][3]-m796[i][2])).subs(dict(rep796, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p796 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m796)):
        if not (m796[i][3]-m796[i][2]).subs(rep796) > 0:
                print('Bounds ' + str(i+1) + ' in p796 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p796 &= ' + latex(p796))
longLatexOutput.append('p796=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m796]).lstrip('+') + '\\right)')


###################################################################################
# 7.6.1 animal: a <= pi.  Sensor: s > pi. Condition: a <= s - pi and a < 2*pi - s #
###################################################################################

m761 = [ [2*r*sin(a/2),               g4, pi/2,             pi/2 + s/2  ],
         [r*sin(a/2), g4, pi/2 + s/2, a/2 + pi/2 + s/2       ] ]

p761 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m761] ).simplify().trigsimp()



rep761 = {s:3*pi/2, a:pi/3} # Replacement values in range


# Define conditions for model
cond761 = [a <= pi, s >= pi/2, a <= s - pi, a <= 2*pi - s]
# Confirm replacements
if not all([c.subs(rep761) for c in cond761]):
        print('rep761 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p761.subs(dict(rep761, **r1)) <= 2:
        print('Total p761 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m761)):
        if not integrate(m761[i][0], m761[i][1:]).subs(dict(rep761, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p761 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m761)):
        if not 0 <= (integrate(m761[i][0], m761[i][1:])/(m761[i][3]-m761[i][2])).subs(dict(rep761, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p761 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m761)):
        if not (m761[i][3]-m761[i][2]).subs(rep761) > 0:
                print('Bounds ' + str(i+1) + ' in p761 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p761 &= ' + latex(p761))
longLatexOutput.append('p761=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m761]).lstrip('+') + '\\right)')



##################################################################################
# 7.6.2 animal: a <= pi.  Sensor: s > pi. Condition: a > s - pi and a < 2*pi - s #
##################################################################################


m762 = [ [ 2*r*sin(a/2),                       g4, pi/2,         s/2 + pi/2 - a/2       ],
         [ r*sin(a/2) + r*sin(s/2 + pi/2 - g4), g4, s/2 + pi/2 - a/2, s/2 + pi/2      ], 
         [ r*sin(a/2),                         g4, s/2 + pi/2,         s/2 + pi/2 + a/2] ]

p762 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m762] ).simplify().trigsimp()



rep762 = {s:5*pi/4, a:pi/2} # Replacement values in range

# Define conditions for model
cond762 = [a <= pi, s >= pi, a >= s - pi, a <= 2*pi - s]
# Confirm replacements
if not all([c.subs(rep762) for c in cond762]):
        print('rep762 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p762.subs(dict(rep762, **r1)) <= 2:
        print('Total p762 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m762)):
        if not integrate(m762[i][0], m762[i][1:]).subs(dict(rep762, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p762 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m762)):
        if not 0 <= (integrate(m762[i][0], m762[i][1:])/(m762[i][3]-m762[i][2])).subs(dict(rep762, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p762 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m762)):
        if not (m762[i][3]-m762[i][2]).subs(rep762) > 0:
                print('Bounds ' + str(i+1) + ' in p762 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p762 &= ' + latex(p762))
longLatexOutput.append('p762=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m762]).lstrip('+') + '\\right)')


#########################################################################################
# 7.6.3 animal: a <= pi.  Sensor: s > pi. Condition: a < s - pi, 2pi - s < a < 4pi - 2s #
#########################################################################################


m763 = [ [ 2*r*sin(a/2), g4, pi/2,               s/2 + pi/2             ],
         [ r*sin(a/2),   g4, s/2 + pi/2,         5*pi/2 - a/2 - s/2     ], 
         [ 2*r*sin(a/2), g4, 5*pi/2 - a/2 - s/2, 3*pi/2]                ]

p763 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m763] ).simplify().trigsimp()



rep763 = {s:3*pi/2 + 0.1, a:pi/2} # Replacement values in range

# Define conditions for model
cond763 = [a <= pi, s >= pi, a <= s - pi, a >= 2*pi - s, a <= 4*pi - 2*s]
# Confirm replacements
if not all([c.subs(rep763) for c in cond763]):
        print('rep763 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p763.subs(dict(rep763, **r1)) <= 2:
        print('Total p763 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m763)):
        if not integrate(m763[i][0], m763[i][1:]).subs(dict(rep763, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p763 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m763)):
        if not 0 <= (integrate(m763[i][0], m763[i][1:])/(m763[i][3]-m763[i][2])).subs(dict(rep763, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p763 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m763)):
        if not (m763[i][3]-m763[i][2]).subs(rep763) > 0:
                print('Bounds ' + str(i+1) + ' in p763 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p763 &= ' + latex(p763))
longLatexOutput.append('p763=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m763]).lstrip('+') + '\\right)')

#########################################################################################
# 7.6.4 animal: a <= pi.  Sensor: s > pi. Condition: a > s - pi, 2pi - s < a < 4pi - 2s #
#########################################################################################




m764 = [ [ 2*r*sin(a/2),                        g4, pi/2,               s/2 + pi/2 - a/2  ],
         [ r*sin(a/2) + r*cos(g4 - s/2), g4, s/2 + pi/2 - a/2,   s/2 + pi/2        ], 
         [ r*sin(a/2),                          g4, s/2 + pi/2,         5*pi/2 - a/2 - s/2],
         [ 2*r*sin(a/2),                        g4, 5*pi/2 - a/2 - s/2, 3*pi/2            ] ]

p764 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m764] ).simplify().trigsimp()



rep764 = {s:5*pi/4, a:7*pi/8} # Replacement values in range

# Define conditions for model
cond764 = [a <= pi, s >= pi, a >= s - pi,  a >= 2*pi - s, a <= 4*pi - 2*s]
# Confirm replacements
if not all([c.subs(rep764) for c in cond764]):
        print('rep764 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p764.subs(dict(rep764, **r1)) <= 2:
        print('Total p764 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m764)):
        if not integrate(m764[i][0], m764[i][1:]).subs(dict(rep764, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p764 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m764)):
        if not 0 <= (integrate(m764[i][0], m764[i][1:])/(m764[i][3]-m764[i][2])).subs(dict(rep764, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p764 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m764)):
        if not (m764[i][3]-m764[i][2]).subs(rep764) > 0:
                print('Bounds ' + str(i+1) + ' in p764 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p764 &= ' + latex(p764))
longLatexOutput.append('p764=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m764]).lstrip('+') + '\\right)')

###############################################################################
# 7.6.5 animal: a <= pi.  Sensor: s > pi. Condition: a < s - pi, 4pi - 2s < a #
###############################################################################


m765 = [ [ 2*r*sin(a/2), g4, pi/2, 3*pi/2  ]]

p765 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m765] ).simplify().trigsimp()



rep765 = {s:19*pi/10, a:pi/2} # Replacement values in range

# Define conditions for model
cond765 = [a <= pi, s >= pi, a <= s - pi, a >= 4*pi - 2*s]
# Confirm replacements
if not all([c.subs(rep765) for c in cond765]):
        print('rep765 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p765.subs(dict(rep765, **r1)) <= 2:
        print('Total p765 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m765)):
        if not integrate(m765[i][0], m765[i][1:]).subs(dict(rep765, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p765 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m765)):
        if not 0 <= (integrate(m765[i][0], m765[i][1:])/(m765[i][3]-m765[i][2])).subs(dict(rep765, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p765 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m765)):
        if not (m765[i][3]-m765[i][2]).subs(rep765) > 0:
                print('Bounds ' + str(i+1) + ' in p765 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p765 &= ' + latex(p765))
longLatexOutput.append('p765=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m765]).lstrip('+') + '\\right)')


###################################################################################
# 7.6.6 animal: a <= pi.  Sensor: s > pi. Condition: a > s - pi and 4*pi - 2s < a #
###################################################################################


m766 = [ [ 2*r*sin(a/2),                        g4, pi/2,         s/2 + pi/2 - a/2       ],
         [ r*sin(a/2) + r*sin(s/2 + pi/2 - g4), g4, s/2 + pi/2 - a/2, 5*pi/2 - a/2 - s/2  ], 
         [ 2*r*sin(a/2),                        g4, 5*pi/2 - a/2 - s/2,         3*pi/2] ]

p766 = pi**-1 * sum( [integrate(x[0], x[1:]) for x in m766] ).simplify().trigsimp()



rep766 = {s:15*pi/8-0.1, a:7*pi/8} # Replacement values in range

# Define conditions for model
cond766 = [a <= pi, s >= pi, a >= s - pi, a >= 4*pi - 2*s]
# Confirm replacements
if not all([c.subs(rep766) for c in cond766]):
        print('rep766 incorrect')


# is average profile in range 0r-2r?
if not 0 <= p766.subs(dict(rep766, **r1)) <= 2:
        print('Total p766 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m766)):
        if not integrate(m766[i][0], m766[i][1:]).subs(dict(rep766, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p766 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m766)):
        if not 0 <= (integrate(m766[i][0], m766[i][1:])/(m766[i][3]-m766[i][2])).subs(dict(rep766, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p766 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m766)):
        if not (m766[i][3]-m766[i][2]).subs(rep766) > 0:
                print('Bounds ' + str(i+1) + ' in p766 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.append('p766 &= ' + latex(p766))
longLatexOutput.append('p766=&\\frac{1}{\pi} \left(' + ' '.join(['+\int_{'+latex(x[2])+'}^{'+latex(x[3])+'}'+latex(x[0])+'\;\mathrm{d}'+latex(x[1]) for x in m766]).lstrip('+') + '\\right)')




##################################################################################################################

####################
## Run tests     ###
####################

# create gas model object
gas = 2*r
p74 = 2*r*sin(a/2)
p711 = r*(2+s)/pi


# for each model run through every adjacent model. 
# Contains duplicatea but better for avoiding missed comparisons.
# Also contains replacement s->a and a->s just in case. 


# I'm pretty sure 793 and 796 are currently equal, even though this can't find that. 
# p791 not right for values close to zero?

allComps = [
['gas', 'p731', {s:2*pi}],
['gas', 'p74', {a:pi}],
['p731', 'gas', {s:2*pi}],
['p731', 'p75', {s:pi}],
['p731', 'p732',{a:3*pi-s}],
['p731', 'p732',{s:3*pi-a}],
['p732', 'p731',{a:3*pi-s}],
['p732', 'p731',{s:3*pi-a}],
['p732', 'p764',{a:pi}],
['p732', 'p766',{a:pi}],
['p732', 'p782',{s:pi}],
['p75', 'p731', {s:pi, a:2*pi}],
['p75','p711',{s:pi/2,a:2*pi}],
['p75','p782',{a:2*pi}],
['p773','p711',{a:2*pi}],
['p773','p772',{a:2*pi-s}],
['p773','p772',{s:2*pi-a}],
['p773','p782',{s:pi/2}],
['p772','p773',{a:2*pi-s}], 
['p772','p773',{s:2*pi-a}],
['p772','p771',{s:pi-a/2}],
['p772','p771',{a:2*pi-2*s}],
['p771','p772',{s:2*pi-2*a}],
['p771','p772',{a:2*pi-2*s}],
['p771','p796',{a:pi}],
['p783','p732',{s:pi}],
['p783','p782',{a:3*pi-2*s}],
['p783','p782',{s:3*pi/2-a/2}],
['p783','p75',{a:2*pi}],
['p782','p773',{s:pi/2}],
['p782','p781',{a:2*pi-s}],
['p782','p781',{s:2*pi-a}],
['p782','p783',{a:3*pi-2*s}],
['p782','p783',{s:3*pi/2-a/2}],
['p781','p772',{s:pi/2}],
['p781','p782',{s:2*pi-a}],
['p781','p782',{a:2*pi-s}],
['p781','p7101',{a:pi}],
['p791','p792',{s:a}],
['p791','p792',{a:s}],
['p791','p794',{s:pi/2-a/2}],
['p791','p794',{a:pi-2*s}],
['p792','p791',{a:s}],
['p792','p791',{s:a}],
['p792','p793',{a:2*s}],
['p792','p793',{s:a/2}],
['p792','p795',{a:pi-2*s}],
['p792','p795',{s:pi/2-a/2}],
['p793','p792',{a:2*s}],
['p793','p792',{s:a/2}],
['p793','p796',{a:pi-2*s}],
['p793','p796',{s:pi/2-a/2}],
['p794','p791',{a:pi-2*s}],
['p794','p791',{s:pi/2-a/2}],
['p794','p795',{s:a}],
['p794','p795',{a:s}],
['p794','p7102',{s:pi/2}],
['p795','p794',{s:a}],
['p795','p794',{a:s}],
['p795','p792',{s:pi/2-a/2}],
['p795','p792',{a:pi-2*s}],
['p795','p796',{a:2*s}],
['p795','p796',{s:a/2}],
['p795','p7101',{s:pi/2}],
['p796','p793',{s:pi/2-a/2}],
['p796','p793',{a:pi-2*s}],
['p796','p795',{a:2*s}],
['p796','p795',{s:a/2}],
['p796','p771',{a:pi}],
['p7101','p795',{s:pi/2}],
['p7101','p7102',{a:s}],
['p7101','p7102',{s:a}],
['p7101','p781',{a:pi}],
['p7102','p7101',{a:s}],
['p7102','p7101',{s:a}],
['p7102','p794',{s:pi/2}],
['p7102','p7103',{a:2*s-pi}],
['p7102','p7103',{s:a/2+pi/2}],
['p7103','p7102',{s:a/2+pi/2}],
['p7103','p7102',{a:2*s-pi}],
['p7103','p762',{s:pi}],
['p761','p762',{a:s-pi}],
['p761','p762',{s:a+pi}],
['p761','p763',{s:2*pi-a}],
['p761','p763',{a:2*pi-s}],
['p762','p761',{a:s-pi}],
['p762','p761',{s:a+pi}],
['p762','p764',{s:2*pi-a}],
['p762','p764',{a:2*pi-s}],
['p762','p7103',{s:pi}],
['p763','p761',{s:2*pi-a}],
['p763','p761',{a:2*pi-s}],
['p763','p764',{a:s-pi}],
['p763','p764',{s:a+pi}],
['p763','p765',{a:4*pi-2*s}],
['p763','p765',{s:2*pi-a/2}],
['p764','p762',{s:2*pi-a}],
['p764','p762',{a:2*pi-s}],
['p764','p763',{a:s-pi}],
['p764','p763',{s:a+pi}],
['p764','p766',{a:4*pi-2*s}],
['p764','p766',{s:2*pi-a/2}],
['p764','p732',{a:pi}],
['p765','p763',{a:4*pi-2*s}],
['p765','p763',{s:2*pi-a/2}],
['p765','p766',{s:2*pi-a}],
['p765','p766',{a:2*pi-s}],
['p765','p74',{s:2*pi}],
['p766','p764',{a:4*pi-2*s}],
['p766','p764',{s:2*pi-a/2}],
['p766','p765',{s:a+pi}],
['p766','p765',{a:s-pi}],
['p766','p732',{a:pi}],
['p74','p765',{s:2*pi}],
['p74','gas',{a:pi, s:2*pi}]
]

# Run through all the comparisons. Need simplify(). Even together() gives some false negatives.

checkFile = open('/home/tim/Dropbox/PhD/Analysis/REM-chapter/checksFile.tex','w')

checkFile.write('All checks evaluated.\nTim Lucas - ' + str(datetime.now()) + '\n')
checkFile.write('I\'m pretty sure 793 and 796 are currently equal, even though this can\'t find that.')
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
y7102Range = [p7102.subs({r:1, s:pi/2, a:i}).n() for i in xRange]
plot7102 = pl.plot(xRange, y7102Range)
pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p7102Profile.pdf')
pl.close()

y794Range = [p794.subs({r:1, s:pi/2, a:i}).n() for i in xRange]
plot794 = pl.plot(xRange, y794Range)
pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p794Profile.pdf')
pl.close()



#pl.savefig('/home/tim/Dropbox/PhD/Analysis/REM-chapter/imgs/p731Profile.pdf')
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
        latexFile.write( '\\[' + latexOutput[i] + '\\]\n')

latexFile.close()

'''
############################
### Viz regions         ####
############################

# Find the unique set of model names
# Crappy coding. But uses all comparison list
allModels = list(set([x[1] for x in allComps]))

# convert to name of condition objects
allConds = [ 'cond' + x[1:] for x in allModels]
# Remove some that don't have conditions (i.e. simple models)
allConds.remove('cond74')
allConds.remove('condas')
allConds.remove('cond711')




[eval(x) for x in allConds]
eval(allConds[0])



parInc = 0.1 # parameter increment
aGrid = sGrid = np.arange(parInc,2*pi,parInc)

modelRegions = np.zeros((len(aGrid), len(sGrid)), dtype=np.int8 )


for m in range(len(allConds)):
        for i in range(len(aGrid)):
                for j in range(len(sGrid)):
                        if all([eval(allConds[m])[c].subs({s:sGrid[j], a:aGrid[i]}) for c in range(len(eval(allConds[m])))]):
                                modelRegions[i,j] = m



rescaled = (255.0 / modelRegions.max() * (modelRegions - modelRegions.min())).astype(np.uint8)

im = Im.fromarray(rescaled)

im.save('/home/tim/Dropbox/test.png')



'''
















