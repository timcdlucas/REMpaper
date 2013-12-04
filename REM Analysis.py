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
import ImageDraw as ImD

# Use LaTeX printing
%load_ext sympyprinting 
# Make LaTeX output white. Because I use a dark theme
init_printing(forecolor="White") 


# Load symbols used for symbolic maths
s, a, r, g1, g2, g3, g4 = symbols('theta_s theta_a r gamma_1 gamma_2 gamma_3 gamma_4', positive=True)
r1 = {r:1} # useful for lots of checks


# A list that will recieve the latex ouput of each model

latexOutput = []
longLatexOutput = []

###################################################
# 7.3. animal: a > pi.  sensor: s > pi            #
###################################################

m73 = [ [2*r,                 g4, pi/2, s/2],
        [r + r*cos(g4 - s/2), g4, s/2,  pi ] ]

p73 = ((2*integrate(m73[0][0], m73[0][1:]) + 2*integrate(m73[1][0], m73[1][1:]))/pi).simplify()


# Replacement values in range
rep73 = {s:3*pi/2} 

# Define conditions for model
cond73 = [pi <= s]
# Confirm replacements
if not all([c.subs(rep73) for c in cond73]):
        print('rep73 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p73.subs(dict(rep73, **r1)) <= 2:
        print('Total p73 not in 0, 2r')

# Are the individuals integrals >0r
for i in range(len(m73)):
        if not integrate(m73[i][0], m73[i][1:]).subs(dict(rep73, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p73 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m73)):
        if not 0 <= (integrate(m73[i][0], m73[i][1:])/(m73[i][3]-m73[i][2])).subs(dict(rep73, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p73 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m73)):
        if not (m73[i][3]-m73[i][2]).subs(rep73) > 0:
                print('Bounds ' + str(i+1) + ' in p73 has lower bounds bigger than upper bounds')
        
# Plot function

xRange = np.arange(pi,2*pi, 0.01)
yRange = [p73.subs({r:1, s:i}).n() for i in xRange]
plot73 = pl.plot(xRange, yRange)
pl.savefig('/home/tim/Dropbox/phd/Analysis/REM-chapter/imgs/p73Profile.pdf')
pl.close()

# LaTeX output

latexOutput.insert(0,'p73 &= ' + latex(p73))
longLatexOutput.insert(0,'p73=\\frac{1}{\pi} \left(2\int_{'+latex(m73[0][2])+'}^{'+latex(m73[0][3])+'}'+latex(m73[0][0])+'\;\mathrm{d}'+latex(m73[0][1])+'+2\int_{'+latex(m73[1][2])+'}^{'+latex(m73[1][3])+'}'+latex(m73[1][0])+'\;\mathrm{d}'+latex(m73[1][1])+'\\right)')


########################################################
# 7.5. animal: a = 2*pi.   sensor:  pi/2 <= s <= pi      #
########################################################

m75 = [ [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2],
        [r - r*cos(g3 - s),    g3, 0, s - pi/2],
        [r,                    g3, s - pi/2, pi/2] ]

p75 = pi**-1 * (2*integrate(m75[0][0], m75[0][1:]) + 2*integrate(m75[1][0], m75[1][1:]) + integrate(m75[2][0], m75[2][1:])).trigsimp().simplify()


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
plot75 = pl.plot(xRange, yRange)
pl.savefig('/home/tim/Dropbox/phd/Analysis/REM-chapter/imgs/p75Profile.pdf')
pl.close()

# LaTeX output

latexOutput.insert(1,'p75 &= ' + latex(p75))
longLatexOutput.insert(1,'p75=\\frac{1}{\pi} \left(\int_{'+latex(m75[0][2])+'}^{'+latex(m75[0][3])+'}'+latex(m75[0][0])+'\;\mathrm{d}'+latex(m75[0][1])+'+\int_{'+latex(m75[1][2])+'}^{'+latex(m75[1][3])+'}'+latex(m75[1][0])+'\;\mathrm{d}'+latex(m75[1][1])+'+\int_{'+latex(m75[2][2])+'}^{'+latex(m75[2][3])+'}'+latex(m75[2][0])+'\;\mathrm{d}'+latex(m75[2][1])+'\\right)')


"""
Multiple models in 7.7.
animal: a>pi.  Sensor: s <= pi/2
"""


###################################################################
# 7.7.1 animal: a>pi.  Sensor: s <= pi/2. Condition: a/2 <= pi - s  #
###################################################################

m771 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s,          pi/2],
         [r*sin(g2),            g2, s,          pi - s/2],
         [r,                    g3, 0,          s] ]


p771 = pi**-1 * (integrate(m771[0][0], m771[0][1:]) + integrate(m771[1][0], m771[1][1:]) + integrate(m771[2][0], m771[2][1:]) + integrate(m771[3][0], m771[3][1:])).simplify().trigsimp()


rep771 = {s:pi/9, a:10*pi/9} # Replacement values in range

# Define conditions for model
cond771 = [s <= pi/2, a >= pi, a/2 <= pi - s]
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

latexOutput.insert(2,'p771 &= ' + latex(p771))
longLatexOutput.insert(2,'p771 &= \\frac{1}{\pi} \left(\int_{'+latex(m771[0][2])+'}^{'+latex(m771[0][3])+'}'+latex(m771[0][0])+'\;\mathrm{d}'+latex(m771[0][1])+'+\int_{'+latex(m771[1][2])+'}^{'+latex(m771[1][3])+'}'+latex(m771[1][0])+'\;\mathrm{d}'+latex(m771[1][1])+'+\int_{'+latex(m771[2][2])+'}^{'+latex(m771[2][3])+'}'+latex(m771[2][0])+'\;\mathrm{d}'+latex(m771[2][1])+'+\int_{'+latex(m771[3][2])+'}^{'+latex(m771[3][3])+'}'+latex(m771[3][0])+'\;\mathrm{d}'+latex(m771[3][1])+'\\right)')


###############################################################################
# 7.7.2 animal: a>pi.  Sensor: s <= pi/2. Condition:  pi - s <= a/2 <= pi - s/2  #
###############################################################################

m772 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s,          pi/2],
         [r*cos(g1 - s/2),      g1, pi/2 - s/2, a/2 + s/2 - pi/2],
         [r,                    g3, 0,          s] ]

p772 = pi**-1 * (integrate(m772[0][0], m772[0][1:]) + 2*integrate(m772[1][0], m772[1][1:]) + integrate(m772[2][0], m772[2][1:]) + integrate(m772[3][0], m772[3][1:])).simplify().trigsimp()


rep772 = {s:3*pi/8, a:3*pi/2} # Replacement values in range

# Define conditions for model
cond772 = [a >= pi, s <= pi/2, pi - s <= a/2, a/2 <= pi - s/2]
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

latexOutput.insert(3,'p772 &= ' + latex(p772))
longLatexOutput.insert(3,'p772 &= \\frac{1}{\pi} \left(\int_{'+latex(m772[0][2])+'}^{'+latex(m772[0][3])+'}'+latex(m772[0][0])+'\;\mathrm{d}'+latex(m772[0][1])+'+\int_{'+latex(m772[1][2])+'}^{'+latex(m772[1][3])+'}'+latex(m772[1][0])+'\;\mathrm{d}'+latex(m772[1][1])+'+\int_{'+latex(m772[2][2])+'}^{'+latex(m772[2][3])+'}'+latex(m772[2][0])+'\;\mathrm{d}'+latex(m772[2][1])+'+\int_{'+latex(m772[3][2])+'}^{'+latex(m772[3][3])+'}'+latex(m772[3][0])+'\;\mathrm{d}'+latex(m772[3][1])+'\\right)')



###############################################################################
# 7.7.3 animal: a>pi.  Sensor: s <= pi/2. Condition: pi - s/2 <= a/2            #
###############################################################################

m773 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s,          pi/2],
         [r*cos(g1 - s/2),      g1, pi/2 - s/2, 3*pi/2 - a/2 + s/2],
         [r,                    g3, 0,          s   ] ]

p773 = pi**-1 * (integrate(m773[0][0], m773[0][1:]) + 2*integrate(m773[1][0], m773[1][1:]) + integrate(m773[2][0], m773[2][1:]) + integrate(m773[3][0], m773[3][1:])).simplify().trigsimp()


rep773 = {s:3*pi/8, a:29*pi/16} # Replacement values in range

# Define conditions for model
cond773 = [a >= pi, s <= pi/2, pi - s/2 <= a/2]
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

latexOutput.insert(4,'p773 &= ' + latex(p773))
longLatexOutput.insert(4,'p773 &= \\frac{1}{\pi} \left(\int_{'+latex(m773[0][2])+'}^{'+latex(m773[0][3])+'}'+latex(m773[0][0])+'\;\mathrm{d}'+latex(m773[0][1])+'+\int_{'+latex(m773[1][2])+'}^{'+latex(m773[1][3])+'}'+latex(m773[1][0])+'\;\mathrm{d}'+latex(m773[1][1])+'+\int_{'+latex(m773[2][2])+'}^{'+latex(m773[2][3])+'}'+latex(m773[2][0])+'\;\mathrm{d}'+latex(m773[2][1])+'+\int_{'+latex(m773[3][2])+'}^{'+latex(m773[3][3])+'}'+latex(m773[3][0])+'\;\mathrm{d}'+latex(m773[3][1])+'\\right)')





"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Two models in 7.8.
animal: pi <= a .  Sensor: pi/2 <= s <= pi
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#################################################################################
# 7.8.1 animal:  a > pi.  Sensor: pi/2 <= s <= pi. Condition: a/2 <= pi - s/2   #
#################################################################################

m781 = [ [2*r*sin(a/2),         g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s/2,        pi/2],
         [r*cos(g1 - s/2),      g1, pi/2 - s/2, a/2 + s/2 -pi/2],
         [r,                    g3, 0,          s   ] ]

p781 = pi**-1 * (integrate(m781[0][0], m781[0][1:]) + integrate(m781[1][0], m781[1][1:]) + integrate(m781[2][0], m781[2][1:]) + integrate(m781[3][0], m781[3][1:])).simplify()


rep781 = {s:3*pi/4, a:15*pi/8} # Replacement values in range

# Define conditions for model
cond781 = [a > pi,  s <= pi/2, a/2 <= pi - s/2]
# Confirm replacements
if not all([c.subs(rep781) for c in cond781]):
        print('rep781 incorrect')

# is average profile in range 0r-2r?
if not 0 <= p781.subs(dict(rep781, **r1)) <= 2:
        print('Total p781 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m781)):
        if not integrate(m781[i][0], m781[i][1:]).subs(dict(rep781, **r1)) > 0:
                print('Integral ' + str(i+1) + ' in p781 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m781)):
        if not 0 <= (integrate(m781[i][0], m781[i][1:])/(m781[i][3]-m781[i][2])).subs(dict(rep781, **r1)) <= 2:
                print('Integral ' + str(i+1) + ' in p781 has averaged integral outside 0<p<2r')
                
# Are the bounds the correct way around
for i in range(len(m781)):
        if not (m781[i][3]-m781[i][2]).subs(rep781) > 0:
                print('Bounds ' + str(i+1) + ' in p781 has lower bounds bigger than upper bounds')        

# LaTeX output

latexOutput.insert(5,'p781 &= ' + latex(p781))
longLatexOutput.insert(5,'p781 &= \\frac{1}{\pi} \left(\int_{'+latex(m781[0][2])+'}^{'+latex(m781[0][3])+'}'+latex(m781[0][0])+'\;\mathrm{d}'+latex(m781[0][1])+'+\int_{'+latex(m781[1][2])+'}^{'+latex(m781[1][3])+'}'+latex(m781[1][0])+'\;\mathrm{d}'+latex(m781[1][1])+'+\int_{'+latex(m781[2][2])+'}^{'+latex(m781[2][3])+'}'+latex(m781[2][0])+'\;\mathrm{d}'+latex(m781[2][1])+'+\int_{'+latex(m781[3][2])+'}^{'+latex(m781[3][3])+'}'+latex(m781[3][0])+'\;\mathrm{d}'+latex(m781[3][1])+'\\right)')



#################################################################################
# 7.8.2 animal: a > pi.  Sensor: pi/2 <= s <= pi. Condition: a/2 > pi - s/2  #
#################################################################################


m782 = [ [2*r*sin(a/2),         g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s/2,        pi/2],
         [r*cos(g1 - s/2),      g1, pi/2 - s/2, 3*pi/2 + s/2 - a/2],
         [2*r*sin(s/2)*sin(g1), g1, pi/2,       3*pi/2 + s/2 - a/2],
         [r,                    g3, 0,          s   ] ]

p782 = pi**-1 * (integrate(m782[0][0], m782[0][1:]) + integrate(m782[1][0], m782[1][1:]) + integrate(m782[2][0], m782[2][1:]) + integrate(m782[3][0], m782[3][1:]) + integrate(m782[4][0], m782[4][1:])).simplify()


rep782 = {s:3*pi/4, a:15*pi/8} # Replacement values in range

# Define conditions for model
cond782 = [a > pi, s <= pi/2, a/2 >= pi - s/2]
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

latexOutput.insert(6,'p782 &= ' + latex(p782))
longLatexOutput.insert(6,'p782 &= \\frac{1}{\pi} \left(\int_{'+latex(m782[0][2])+'}^{'+latex(m782[0][3])+'}'+latex(m782[0][0])+'\;\mathrm{d}'+latex(m782[0][1])+'+\int_{'+latex(m782[1][2])+'}^{'+latex(m782[1][3])+'}'+latex(m782[1][0])+'\;\mathrm{d}'+latex(m782[1][1])+'+\int_{'+latex(m782[2][2])+'}^{'+latex(m782[2][3])+'}'+latex(m782[2][0])+'\;\mathrm{d}'+latex(m782[2][1])+'+\int_{'+latex(m782[3][2])+'}^{'+latex(m782[3][3])+'}'+latex(m782[3][0])+'\;\mathrm{d}'+latex(m782[3][1])+'+\int_{'+latex(m782[4][2])+'}^{'+latex(m782[4][3])+'}'+latex(m782[4][0])+'\;\mathrm{d}'+latex(m782[4][1])+'\\right)')





###############################################################################


"""
Complex profiles for a <= pi/2 
"""

# p-l-r for g1 profil. Calculated by AE in fig 7.2 minus AE in fig 7.3
p1 = 2*r*sin(s - 3*pi/2 + g1)*sin((g1 - pi/2 - s + a)/2) -  \
     2*r*sin((pi - a - 2*g1 + s)/4)*sin((pi - a +2*g1 - s)/4).simplify()

# p-l for g1 profiles
p2 = 2*r*sin(s/2)*sin(g1) - 2*r*sin((pi - a - 2*g1 + s)/4)*sin((pi - a + 2*g1 - s)/4)

# p-l for g1 profile. 
p3 = r*sin(g2) - 2*r*sin(g2/2 - a/4)*sin((pi - a + 2*g2 - s)/4).simplify()



#########################################################################################
# 7.10.1 animal: a <= pi.  Sensor: pi/2 <= s <= pi. Condition: a > s and a/2 <= s- pi/2 #
#########################################################################################


m7101 = [ [2*r*sin(s/2)*sin(g1),              g1, pi/2 - a/2 + s/2, pi/2            ],
          [p2,                                g1, s/2,              pi/2 - a/2 + s/2],
          [r*sin(a/2) + r*sin(s - pi/2 - g3), g3, 0,                s - pi/2        ],
          [r*sin(a/2),                        g3, s - pi/2,         s - pi/2 + a/2  ] ]

p7101 = pi**-1 * (integrate(m7101[0][0], m7101[0][1:]) + integrate(m7101[1][0], m7101[1][1:]) + integrate(m7101[2][0], m7101[2][1:]) + integrate(m7101[3][0], m7101[3][1:]) ).simplify()


rep7101 = {s:5*pi/8, a:6*pi/8} # Replacement values in range

# Define conditions for model
cond7101 = [a <= pi, pi/2 <= s, s <= pi, a/2 >= s/2, a/2 >= s - pi/2]
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

latexOutput.insert(7,'p7101 &= ' + latex(p7101))
longLatexOutput.insert(7,'p7101=\\frac{1}{\pi} \left(\int_{'+latex(m7101[0][2])+'}^{'+latex(m7101[0][3])+'}'+latex(m7101[0][0])+'\;\mathrm{d}'+latex(m7101[0][1])+'+\int_{'+latex(m7101[1][2])+'}^{'+latex(m7101[1][3])+'}'+latex(m7101[1][0])+'\;\mathrm{d}'+latex(m7101[1][1])+'+\int_{'+latex(m7101[2][2])+'}^{'+latex(m7101[2][3])+'}'+latex(m7101[2][0])+'\;\mathrm{d}'+latex(m7101[2][1])+'+\int_{'+latex(m7101[3][2])+'}^{'+latex(m7101[3][3])+'}'+latex(m7101[3][0])+'\;\mathrm{d}'+latex(m7101[3][1])+'\\right)')






#########################################################################################
# 7.10.2 animal: a <= pi.  Sensor: pi/2 <= s <= pi. Condition: a > s and a/2 <= s- pi/2 #
#########################################################################################

# THIS IS A PLACEHOLDER!!!!

m7102 = [ [2*r*sin(s/2)*sin(g1),              g1, pi/2 - a/2 + s/2, pi/2            ],
          [2*r*sin(s/2)*sin(g1) - p2,         g1, s/2,              pi/2 - a/2 + s/2],
          [2*r*sin(a/2),                      g3, 0,                s - pi/2 - a/2  ],
          [r*sin(a/2) + r*sin(s - pi/2 - g3), g3, s - pi/2 - a/2,   s - pi/2        ],
          [r*sin(a/2),                        g3, s - pi/2,         s - pi/2 - a/2  ] ]

p7102 = pi**-1 * (integrate(m7102[0][0], m7102[0][1:]) + integrate(m7102[1][0], m7102[1][1:]) + integrate(m7102[2][0], m7102[2][1:]) + integrate(m7102[3][0], m7102[3][1:]) + integrate(m7102[4][0], m7102[4][1:]) ).simplify()


rep7102 = {s:5*pi/8, a:6*pi/8} # Replacement values in range

# Define conditions for model
cond7102 = [a <= pi, pi/2 <= s, s <= pi, a/2 >= s/2, a/2 <= s - pi/2]
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

latexOutput.insert(7,'p7102 &= ' + latex(p7102))
longLatexOutput.insert(7,'p7102=\\frac{1}{\pi} \left(\int_{'+latex(m7102[0][2])+'}^{'+latex(m7102[0][3])+'}'+latex(m7102[0][0])+'\;\mathrm{d}'+latex(m7102[0][1])+'+\int_{'+latex(m7102[1][2])+'}^{'+latex(m7102[1][3])+'}'+latex(m7102[1][0])+'\;\mathrm{d}'+latex(m7102[1][1])+'+\int_{'+latex(m7102[2][2])+'}^{'+latex(m7102[2][3])+'}'+latex(m7102[2][0])+'\;\mathrm{d}'+latex(m7102[2][1])+'+\int_{'+latex(m7102[3][2])+'}^{'+latex(m7102[3][3])+'}'+latex(m7102[3][0])+'\;\mathrm{d}'+latex(m7102[3][1])+'+\int_{'+latex(m7102[4][2])+'}^{'+latex(m7102[4][3])+'}'+latex(m7102[4][0])+'\;\mathrm{d}'+latex(m7102[4][1])+'\\right)')






###################################################################################
# 7.9.1 animal: a <= pi.  Sensor: s <= pi/2. Condition: a <= pi - 2s & a <= s     #
###################################################################################

m791 = [ [p1, g1, pi/2 - s/2 + a/2, pi/2            ],
         [p2, g1, pi/2 - s/2,       pi/2 - s/2 + a/2],
         [p3, g2, s,                s + a/2         ] ]

p791 = pi**-1 * (integrate(m791[0][0], m791[0][1:]) + integrate(m791[1][0], m791[1][1:]) + integrate(m791[2][0], m791[2][1:])).simplify().trigsimp()


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

latexOutput.insert(9,'p791 &= ' + latex(p791))
longLatexOutput.insert(9,'p791=\\frac{1}{\pi} \left(\int_{'+latex(m791[0][2])+'}^{'+latex(m791[0][3])+'}p_1\;\mathrm{d}'+latex(m791[0][1])+'+\int_{'+latex(m791[1][2])+'}^{'+latex(m791[1][3])+'}p_2\;\mathrm{d}'+latex(m791[1][1])+'+\int_{'+latex(m791[2][2])+'}^{'+latex(m791[2][3])+'}p_3\;\mathrm{d}'+latex(m791[2][1])+'\\right)')


#######################################################################################
# 7.9.2 animal: a <= pi.  Sensor: s <= pi/2. Condition: a <= pi - 2s & s <= a <= 2s   #
#######################################################################################

print 'Stuff here! Bounds wrong in new bit'

m792 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 + s/2 - a/2, pi/2            ],
         [p2,                   g1, pi/2 - s/2,       pi/2 + s/2 - a/2],
         [p3,                   g2, s,                s + a/2         ],
         [r*sin(a/2),           g3, 0,                a/2 + s - pi/2   ] ]

p792 = pi**-1 * (integrate(m792[0][0], m792[0][1:]) + integrate(m792[1][0], m792[1][1:]) + integrate(m792[2][0], m792[2][1:]) + integrate(m792[3][0], m792[3][1:])).simplify().trigsimp()


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

latexOutput.insert(10,'p792 &= ' + latex(p792))
longLatexOutput.insert(10,'p792=\\frac{1}{\pi} \left(\int_{'+latex(m792[0][2])+'}^{'+latex(m792[0][3])+'}'+latex(m792[0][0])+'\;\mathrm{d}'+latex(m792[0][1])+'+\int_{'+latex(m792[1][2])+'}^{'+latex(m792[1][3])+'}p_2\;\mathrm{d}'+latex(m792[1][1])+'+\int_{'+latex(m792[2][2])+'}^{'+latex(m792[2][3])+'}p_3\;\mathrm{d}'+latex(m792[2][1])+'\\right)')


##################################################################################
# 7.9.3 animal: a <= pi.  Sensor: s <= pi/2. Condition: a <= pi - 2s &  2s <= a      #
##################################################################################


m793 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2    ],
         [r*sin(g1),            g1, s,          a/2     ],
         [p3,                   g2, a/2,        s + a/2 ],
         [r*sin(a/2),           g3, 0,          a/2 + s -pi/2   ] ]

p793 = pi**-1 * (integrate(m793[0][0], m793[0][1:]) + integrate(m793[1][0], m793[1][1:]) + integrate(m793[2][0], m793[2][1:]) + integrate(m793[3][0], m793[3][1:])).simplify().trigsimp()


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

# LaTeX output

latexOutput.insert(11,'p793 &= ' + latex(p793))
longLatexOutput.insert(11,'p793=\\frac{1}{\pi} \left(\int_{'+latex(m793[0][2])+'}^{'+latex(m793[0][3])+'}p_1\;\mathrm{d}'+latex(m793[0][1])+'+\int_{'+latex(m793[1][2])+'}^{'+latex(m793[1][3])+'p_2\;\mathrm{d}'+latex(m793[1][1])+'+\int_{'+latex(m793[2][2])+'}^{'+latex(m793[2][3])+'}p_3\;\mathrm{d}'+latex(m793[2][1])+'\\right)')

##################################################################################
# 7.9.4 animal: a <= pi.  Sensor: s <= pi/2. Condition: a > pi - 2s &  a <= s       #
##################################################################################

m794 = [ [p1, g1, pi/2 - s/2 + a/2, pi/2            ],
         [p2, g1, pi/2 - s/2,       pi/2 - s/2 + a/2],
         [p3, g2, s,                pi/2            ] ]

p794 = pi**-1 * (integrate(m794[0][0], m794[0][1:]) + integrate(m794[1][0], m794[1][1:]) + integrate(m794[2][0], m794[2][1:])).simplify().trigsimp()

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

latexOutput.insert(12,'p794 &= ' + latex(p794))
longLatexOutput.insert(12,'p794=\\frac{1}{\pi} \left(\int_{'+latex(m794[0][2])+'}^{'+latex(m794[0][3])+'}p_1\;\mathrm{d}'+latex(m794[0][1])+'+\int_{'+latex(m794[1][2])+'}^{'+latex(m794[1][3])+'}p_2\;\mathrm{d}'+latex(m794[1][1])+'+\int_{'+latex(m794[2][2])+'}^{'+latex(m794[2][3])+'}p_3\;\mathrm{d}'+latex(m794[2][1])+'\\right)')

##################################################################################
# 7.9.5 animal: a <= pi.  Sensor: s <= pi/2. Condition: a > pi - 2s &  s <= a <= 2s  #
##################################################################################


m795 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 + s/2 - a/2, pi/2            ],
         [p2,                   g1, pi/2 - s/2,       pi/2 + s/2 - a/2],
         [p3,                   g2, s,                pi/2        ],
         [r*sin(a/2),           g3, 0,                a/2 + s -pi/2   ] ]

p795 = pi**-1 * (integrate(m795[0][0], m795[0][1:]) + integrate(m795[1][0], m795[1][1:]) + integrate(m795[2][0], m795[2][1:]) + integrate(m795[3][0], m795[3][1:])).simplify().trigsimp()


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

latexOutput.insert(13,'p795 &= ' + latex(p795))
longLatexOutput.insert(13,'p795=\\frac{1}{\pi} \left(\int_{'+latex(m795[0][2])+'}^{'+latex(m795[0][3])+'}'+latex(m795[0][0])+'\;\mathrm{d}'+latex(m795[0][1])+'+\int_{'+latex(m795[1][2])+'}^{'+latex(m795[1][3])+'}p_2\;\mathrm{d}'+latex(m795[1][1])+'+\int_{'+latex(m795[2][2])+'}^{'+latex(m795[2][3])+'}p_3\;\mathrm{d}'+latex(m795[2][1])+'\\right)')





##################################################################################
# 7.9.6 animal: a <= pi.  Sensor: s <= pi/2. Condition: a > pi - 2s &  a > 2s      #
##################################################################################



m796 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2            ],
         [r*sin(g2),            g2, s,          a/2             ],
         [p3,                   g2, a/2,        pi/2            ],
         [r*sin(a/2),           g3, 0,          a/2 + s -pi/2   ] ]

p796 = pi**-1 * (integrate(m796[0][0], m796[0][1:]) + integrate(m796[1][0], m796[1][1:]) + integrate(m796[2][0], m796[2][1:]) + integrate(m796[3][0], m796[3][1:])).simplify()


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

latexOutput.insert(14,'p796 &= ' + latex(p796))
longLatexOutput.insert(14,'p796=\\frac{1}{\pi} \left(\int_{'+latex(m796[0][2])+'}^{'+latex(m796[0][3])+'}'+latex(m796[0][0])+'\;\mathrm{d}'+latex(m796[0][1])+'+\int_{'+latex(m796[1][2])+'}^{'+latex(m796[1][3])+'}'+latex(m796[1][0])+'\;\mathrm{d}'+latex(m796[1][1])+'+\int_{'+latex(m796[2][2])+'}^{'+latex(m796[2][3])+'}p_3\;\mathrm{d}'+latex(m796[2][1])+'\\right)')


#############################################
## 7.6 Still not sure how many models etc. ##
#############################################

# Transitions I think I've found

m76Trans = [ pi - a/2 + s/2, s/2 + pi/2, s/2 - a/2 + pi, 5*pi/2 - a/2 - s/2, 5*pi/2 - s/2, 3*pi - s/2 - a/2]

solve(m76Trans[1] - m76Trans[0], a)
solve(m76Trans[3] - m76Trans[1], a)
solve(m76Trans[3] - m76Trans[1], s)



#############################################################################
# 7.6.1 animal: a <= pi.  Sensor: s > pi. Condition: a <= s - pi              #
#############################################################################

m761 = [ [2*r*sin(a/2),               g4, pi/2,             pi/2 - a/2 + s/2  ],
         [r*(sin(a/2) + sin(g4-s/2)), g4, pi/2 - a/2 + s/2, s/2 + pi/2        ],
         [r*sin(a/2),                 g4, s/2+pi/2,         s/2 + pi - a/2    ]  ]

p761 = pi**-1 * (integrate(m761[0][0], m761[0][1:]) + integrate(m761[1][0], m761[1][1:]) + integrate(m761[2][0], m761[2][1:])).simplify()


rep761 = {s:3*pi/2, a:pi/3} # Replacement values in range


# Define conditions for model
cond761 = [a <= pi, s >= pi/2, a <= s - pi]
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

latexOutput.insert(15,'p761 &= ' + latex(p761))
longLatexOutput.insert(15,'p761=\\frac{1}{\pi} \left(2\int_{'+latex(m761[0][2])+'}^{'+latex(m761[0][3])+'}'+latex(m761[0][0])+'\;\mathrm{d}'+latex(m761[0][1])+'+2\int_{'+latex(m761[1][2])+'}^{'+latex(m761[1][3])+'}'+latex(m761[1][0])+'\;\mathrm{d}'+latex(m761[1][1])+'\\right)')



#########################################################################
# 7.6.2 animal: a <= pi.  Sensor: s > pi. Condition: a > s - pi        #
#########################################################################


m762 = [ [ 2*r*sin(a/2),                      g4, pi/2,         s/2       ],
         [ r*sin(a/2) - r*cos(g4),            g4, s/2,          pi        ], 
         [ r*sin(a/2),                        g4, pi,           2*pi - s/2],
         [ 2*r*(sin(a/2) + sin(s/2)*sin(g4)), g4, 2*pi - s/2,   3*pi/2    ]  ]

p762 = pi**-1 * (2*integrate(m762[0][0], m762[0][1:]) + integrate(m762[1][0], m762[1][1:]) + integrate(m762[2][0], m762[2][1:]) + 2*integrate(m762[3][0], m762[3][1:])).simplify()


rep762 = {s:4*pi/3, a:2*pi/3} # Replacement values in range

# Define conditions for model
cond762 = [a <= pi, s >= pi, a >= s - pi]
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

latexOutput.insert(16,'p762 &= ' + latex(p762))
longLatexOutput.insert(16,'p762 &= \\frac{1}{\pi} \left(\int_{'+latex(m762[0][2])+'}^{'+latex(m762[0][3])+'}'+latex(m762[0][0])+'\;\mathrm{d}'+latex(m762[0][1])+'+\int_{'+latex(m762[1][2])+'}^{'+latex(m762[1][3])+'}'+latex(m762[1][0])+'\;\mathrm{d}'+latex(m762[1][1])+'+\int_{'+latex(m762[2][2])+'}^{'+latex(m762[2][3])+'}'+latex(m762[2][0])+'\;\mathrm{d}'+latex(m762[2][1])+'+\int_{'+latex(m762[3][2])+'}^{'+latex(m762[3][3])+'}'+latex(m762[3][0])+'\;\mathrm{d}'+latex(m762[3][1])+'\\right)')




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

allComps = [
['gas', 'p73', {s:2*pi}],
['gas', 'p74', {a:pi}],
['p73', 'gas', {s:2*pi}],
['p73', 'p75', {s:pi, a:2*pi}],
['p73', 'p782',{s:pi}],
['p73', 'p762',{a:pi}],
['p75', 'p73', {s:pi, a:2*pi}],
['p75','p711',{s:pi/2,a:2*pi}],
['p75','p782',{a:2*pi}],
['p773','p711',{a:2*pi}],
['p773','p772',{a:2*pi-s}],
['p773','p772',{s:2*pi-a}],
['p773','p782',{s:pi/2}],
['p772','p772',{a:2*pi-s}],
['p772','p772',{s:2*pi-a}],
['p772','p771',{s:pi-a/2}],
['p772','p771',{a:2*pi-2*s}],
['p772','p782',{s:pi/2}],
['p771','p772',{s:2*pi-2*a}],
['p771','p772',{a:2*pi-2*s}],
['p771','p796',{a:pi}],
['p782','p773',{s:pi/2}],
['p782','p781',{a:2*pi-s}],
['p782','p781',{s:2*pi-a}],
['p782','p75',{a:2*pi}],
['p782','p73',{s:pi}],
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
['p7102','p762',{s:pi}],
['p762','p761',{a:s-pi}],
['p762','p761',{s:a+pi}],
['p762','p7102',{a:pi}],
['p762','p73',{a:pi}],
['p761','p762',{a:s-pi}],
['p761','p762',{s:a+pi}],
['p761','p74',{s:2*pi}],
['p74','p761',{s:2*pi}],
['p74','gas',{a:pi, s:2*pi}]
]

'''
for i in range(len(allComps)):
        if eval(allComps[i][0]).subs(allComps[i][2]).together() == eval(allComps[i][1]).subs(allComps[i][2]).together():
                print(allComps[i][0]+ ' and ' +allComps[i][1]+': OK')
        else:
                print(allComps[i][0]+ ' and ' +allComps[i][1]+': Incorrect')
'''

# Run through all the comparisons. Need simplify(). Even together() gives some false negatives.

checkFile = open('/home/tim/Dropbox/phd/Analysis/REM-chapter/checksFile.tex','w')

checkFile.write('All checks evaluated.\nTim Lucas - ' + str(datetime.now()) + '\n')
for i in range(len(allComps)):
        if eval(allComps[i][0]).subs(allComps[i][2]).simplify() == eval(allComps[i][1]).subs(allComps[i][2]).simplify():
                checkFile.write(allComps[i][0]+ ' and ' +allComps[i][1]+': OK\n')
        else:
                checkFile.write(allComps[i][0]+ ' and ' +allComps[i][1]+': Incorrect\n')

checkFile.close()


# And print to terminal
for i in range(len(allComps)):
        if not eval(allComps[i][0]).subs(allComps[i][2]).simplify() == eval(allComps[i][1]).subs(allComps[i][2]).simplify():
                print allComps[i][0] + ' and ' + allComps[i][1]+': Incorrect\n'

######################
### Write output   ###
######################

# write out full latex model solutions and model statements

latexFile = open('/home/tim/Dropbox/phd/Analysis/REM-chapter/ModelSolutions.tex', 'w')

latexFile.write('% LaTeX output. Solutions of all REM models.\n' + '%Tim Lucas - ' + str(datetime.now()) + '\n')
for i in range(len(latexOutput)):
        latexFile.write( '\\[' + latexOutput[i] + '\\]\n')

latexFile.close()




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











xRange = np.arange(pi,2*pi, 0.01)
yRange = [p73.subs({r:1, s:i}).n() for i in xRange]
plot73 = pl.plot(xRange, yRange)
pl.savefig('/home/tim/Dropbox/phd/Analysis/REM-chapter/imgs/p73Profile.pdf')
pl.close()










