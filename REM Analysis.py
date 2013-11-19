"""
Systematic analysis of REM models
Tim Lucas 
01/10/13
"""


from sympy import *
import numpy as np
import matplotlib.pyplot as pl


# Use LaTeX printing
%load_ext sympyprinting 
# Make LaTeX output white. Because I use a dark theme
init_printing(forecolor="White") 


# Load symbols used for symbolic maths
s, a, r, g1, g2, g3, g4 = symbols('theta_s theta_a r gamma_1 gamma_2 gamma_3 gamma_4')
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

# Confirm replacements
if not pi <= s.subs(rep73):
	print('rep73 incorrect')

# is average profile in range 0r-2r?
if not 0 < p73.subs(dict(rep73, **r1)) < 2:
	print('Total p73 not in 0, 2r')

# Are the individuals integrals >0r
for i in range(len(m73)):
	if not integrate(m73[i][0], m73[i][1:]).subs(dict(rep73, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p73 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m73)):
	if not 0 < (integrate(m73[i][0], m73[i][1:])/(m73[i][3]-m73[i][2])).subs(dict(rep73, **r1)) <= 2:
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

# LaTeX output

latexOutput.insert(0,'p73 = ' + latex(p73))



########################################################
# 7.5. animal: a = 2*pi.   sensor:  pi/2 < s < pi      #
########################################################

m75 = [ [2*r*sin(s/2)*sin(g1), g1, s/2, pi/2],
        [r - r*cos(g3 - s),    g3, 0, s - pi/2],
        [r,                    g3, s - pi/2, pi/2] ]

p75 = pi**-1 * (2*integrate(m75[0][0], m75[0][1:]) + 2*integrate(m75[1][0], m75[1][1:]) + integrate(m75[2][0], m75[2][1:])).trigsimp().simplify()


# Replacement values in range
rep75 = {s:3*pi/4} 

# Confirm replacements
if not pi/2 <= s.subs(rep75) <= pi:
	print('rep75 incorrect')

# is average profile in range 0r-2r?
if not 0 < p75.subs(dict(rep75, **r1)) < 2:
	print('Total p75 not in 0, 2r')

# Are the individuals integrals >0r
for i in range(len(m75)):
	if not integrate(m75[i][0], m75[i][1:]).subs(dict(rep75, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p75 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m75)):
	if not 0 < (integrate(m75[i][0], m75[i][1:])/(m75[i][3]-m75[i][2])).subs(dict(rep75, **r1)) <= 2:
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

# LaTeX output

latexOutput.insert(1,'p75 = ' + latex(p75))



"""
Multiple models in 7.7.
animal: a>pi.  Sensor: s < pi/2
"""


###################################################################
# 7.7.1 animal: a>pi.  Sensor: s < pi/2. Condition: a/2 < pi - s  #
###################################################################

m771 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s,          pi/2],
         [r*sin(g2),            g2, s,          pi - s/2],
         [r,                    g3, 0,          s] ]


p771 = pi**-1 * (integrate(m771[0][0], m771[0][1:]) + integrate(m771[1][0], m771[1][1:]) + integrate(m771[2][0], m771[2][1:]) + integrate(m771[3][0], m771[3][1:])).simplify().trigsimp()


rep771 = {s:pi/9, a:10*pi/9} # Replacement values in range

# Confirm replacements
if not (s.subs(rep771) <= pi/2 and a.subs(rep771) > pi and a.subs(rep771)/2 < pi - s.subs(rep771)):
	print('rep771 incorrect')

# is average profile in range 0r-2r?
if not 0 < p771.subs(dict(rep771, **r1)) < 2:
	print('Total p771 not in 0, 2r')

# Are the individuals integrals >0r
for i in range(len(m771)):
	if not integrate(m771[i][0], m771[i][1:]).subs(dict(rep771, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p771 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m771)):
	if not 0 < (integrate(m771[i][0], m771[i][1:])/(m771[i][3]-m771[i][2])).subs(dict(rep771, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p771 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m771)):
	if not (m771[i][3]-m771[i][2]).subs(rep771) > 0:
		print('Bounds ' + str(i+1) + ' in p771 has lower bounds bigger than upper bounds')	


# LaTeX output

latexOutput.insert(2,'p771 = ' + latex(p771))



###############################################################################
# 7.7.2 animal: a>pi.  Sensor: s < pi/2. Condition:  pi - s < a/2 < pi - s/2  #
###############################################################################

m772 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s,          pi/2],
         [r*cos(g1 - s/2),      g1, pi/2 - s/2, a/2 + s/2 - pi/2],
         [r,                    g3, 0,          s] ]

p772 = pi**-1 * (integrate(m772[0][0], m772[0][1:]) + 2*integrate(m772[1][0], m772[1][1:]) + integrate(m772[2][0], m772[2][1:]) + integrate(m772[3][0], m772[3][1:])).simplify().trigsimp()


rep772 = {s:3*pi/8, a:3*pi/2} # Replacement values in range

# Confirm replacements
if not (a.subs(rep772) > pi and s.subs(rep772) < pi/2 and pi - s.subs(rep772) < a.subs(rep772)/2 < pi - s.subs(rep772)/2):
	print('rep772 incorrect')

# is average profile in range 0r-2r?
if not 0 < p772.subs(dict(rep772, **r1)) < 2:
	print('Total p772 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m772)):
	if not integrate(m772[i][0], m772[i][1:]).subs(dict(rep772, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p772 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m772)):
	if not 0 < (integrate(m772[i][0], m772[i][1:])/(m772[i][3]-m772[i][2])).subs(dict(rep772, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p772 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m772)):
	if not (m772[i][3]-m772[i][2]).subs(rep772) > 0:
		print('Bounds ' + str(i+1) + ' in p772 has lower bounds bigger than upper bounds')	

# LaTeX output

latexOutput.insert(3,'p772 = ' + latex(p772))


###############################################################################
# 7.7.3 animal: a>pi.  Sensor: s < pi/2. Condition: pi - s/2 < a/2            #
###############################################################################

m773 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s,          pi/2],
         [r*cos(g1 - s/2),      g1, pi/2 - s/2, 3*pi/2 - a/2 + s/2],
         [r,                    g3, 0,          s   ] ]

p773 = pi**-1 * (integrate(m773[0][0], m773[0][1:]) + 2*integrate(m773[1][0], m773[1][1:]) + integrate(m773[2][0], m773[2][1:]) + integrate(m773[3][0], m773[3][1:])).simplify().trigsimp()


rep773 = {s:3*pi/8, a:29*pi/16} # Replacement values in range

# Confirm replacements
if not (a.subs(rep773) > pi and s.subs(rep773) <= pi/2 and pi - s.subs(rep773)/2 <= a.subs(rep773)/2):
	print('rep773 incorrect')

# is average profile in range 0r-2r?
if not 0 < p773.subs(dict(rep773, **r1)) < 2:
	print('Total p773 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m773)):
	if not integrate(m773[i][0], m773[i][1:]).subs(dict(rep773, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p773 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m773)):
	if not 0 < (integrate(m773[i][0], m773[i][1:])/(m773[i][3]-m773[i][2])).subs(dict(rep773, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p773 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m773)):
	if not (m773[i][3]-m773[i][2]).subs(rep773) > 0:
		print('Bounds ' + str(i+1) + ' in p773 has lower bounds bigger than upper bounds')	

# LaTeX output

latexOutput.insert(4,'p773 = ' + latex(p773))





"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Two models in 7.8.
animal: pi < a .  Sensor: pi/2 < s < pi
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


#################################################################################
# 7.8.1 animal:  a > pi.  Sensor: pi/2 < s < pi. Condition: a/2 < pi - s/2      #
#################################################################################

m781 = [ [2*r*sin(a/2),         g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s/2,        pi/2],
         [r*cos(g1 - s/2),      g1, pi/2 - s/2, a/2 + s/2 -pi/2],
         [r,                    g3, 0,          s   ] ]

p781 = pi**-1 * (integrate(m781[0][0], m781[0][1:]) + integrate(m781[1][0], m781[1][1:]) + integrate(m781[2][0], m781[2][1:]) + integrate(m781[3][0], m781[3][1:])).simplify()


rep781 = {s:3*pi/4, a:15*pi/8} # Replacement values in range

# Confirm replacements
if not (a.subs(rep781) > pi and s.subs(rep781) <= pi/2 and a.subs(rep781)/2 <= pi - s.subs(rep781)/2):
	print('rep781 incorrect')

# is average profile in range 0r-2r?
if not 0 < p781.subs(dict(rep781, **r1)) < 2:
	print('Total p781 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m781)):
	if not integrate(m781[i][0], m781[i][1:]).subs(dict(rep781, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p781 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m781)):
	if not 0 < (integrate(m781[i][0], m781[i][1:])/(m781[i][3]-m781[i][2])).subs(dict(rep781, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p781 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m781)):
	if not (m781[i][3]-m781[i][2]).subs(rep781) > 0:
		print('Bounds ' + str(i+1) + ' in p781 has lower bounds bigger than upper bounds')	

# LaTeX output

latexOutput.insert(5,'p781 = ' + latex(p781))



#################################################################################
# 7.8.2 animal: a > pi.  Sensor: pi/2 < s < pi. Condition: a/2 > pi - s/2  #
#################################################################################


m782 = [ [2*r*sin(a/2),         g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s/2,        pi/2],
         [r*cos(g1 - s/2),      g1, pi/2 - s/2, 3*pi/2 + s/2 - a/2],
         [2*r*sin(s/2)*sin(g1), g1, pi/2,       3*pi/2 + s/2 - a/2],
         [r,                    g3, 0,          s   ] ]

p782 = pi**-1 * (integrate(m782[0][0], m782[0][1:]) + integrate(m782[1][0], m782[1][1:]) + integrate(m782[2][0], m782[2][1:]) + integrate(m782[3][0], m782[3][1:]) + integrate(m782[4][0], m782[4][1:])).simplify()


rep782 = {s:3*pi/4, a:15*pi/8} # Replacement values in range

# Confirm replacements
if not (a.subs(rep782) > pi and s.subs(rep782) <= pi/2 and a.subs(rep782)/2 >= pi - s.subs(rep782)/2):
	print('rep782 incorrect')

# is average profile in range 0r-2r?
if not 0 < p782.subs(dict(rep782, **r1)) < 2:
	print('Total p782 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m782)):
	if not integrate(m782[i][0], m782[i][1:]).subs(dict(rep782, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p782 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m782)):
	if not 0 < (integrate(m782[i][0], m782[i][1:])/(m782[i][3]-m782[i][2])).subs(dict(rep782, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p782 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m782)):
	if not (m782[i][3]-m782[i][2]).subs(rep782) > 0:
		print('Bounds ' + str(i+1) + ' in p782 has lower bounds bigger than upper bounds')	

# LaTeX output

latexOutput.insert(6,'p782 = ' + latex(p782))





###############################################################################


"""
Complex profiles for a < pi/2 
"""

# p-l-r for g1 profil. Calculated by AE in fig 7.2 minus AE in fig 7.3
p1 = 2*r*sin(s - 3*pi/2 + g1)*sin((g1 - pi/2 - s + a)/2) -  \
     2*r*sin((pi - a - 2*g1 + s)/4)*sin((pi - a +2*g1 - s)/4).simplify()

# p-l for g1 profiles
p2 = 2*r*sin(s/2)*sin(g1) - 2*r*sin((pi - a - 2*g1 + s)/4)*sin((pi - a + 2*g1 - s)/4)

# p-l for g1 profile. 
p3 = r*sin(g2) - 2*r*sin(g2/2 - a/4)*sin((pi - a + 2*g2 - s)/4).simplify()



###############################################################################
# 7.10.1 animal: a < pi.  Sensor: pi/2 < s < pi. Condition: a > s             #
###############################################################################




m7101 = [ [2*r*sin(s/2)*sin(g1),  g1, pi/2 - a/2 + s/2, pi/2           ],
          [p2,                    g1, pi/2 - s/2,       pi/2 - a/2 +s/2],
          [r*sin(a/2),            g2, 0,                s - pi/2 + a/2 ] ]

p7101 = pi**-1 * (integrate(m7101[0][0], m7101[0][1:]) + integrate(m7101[1][0], m7101[1][1:]) + integrate(m7101[2][0], m7101[2][1:]) ).simplify()


rep7101 = {s:5*pi/8, a:7*pi/8} # Replacement values in range

# Confirm replacements
if not (a.subs(rep7101) < pi and pi/2 <= s.subs(rep7101) <= pi and a.subs(rep7101)/2 >= s.subs(rep7101)/2):
	print('rep7101 incorrect')

# is average profile in range 0r-2r?
if not 0 < p7101.subs(dict(rep7101, **r1)) < 2:
	print('Total p7101 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m7101)):
	if not integrate(m7101[i][0], m7101[i][1:]).subs(dict(rep7101, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p7101 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m7101)):
	if not 0 < (integrate(m7101[i][0], m7101[i][1:])/(m7101[i][3]-m7101[i][2])).subs(dict(rep7101, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p7101 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m7101)):
	if not (m7101[i][3]-m7101[i][2]).subs(rep7101) > 0:
		print('Bounds ' + str(i+1) + ' in p7101 has lower bounds bigger than upper bounds')	

# LaTeX output

latexOutput.insert(7,'p7101 = ' + latex(p7101))



###############################################################################
# 7.10.2 animal: a < pi.  Sensor: pi/2 < s < pi. Condition: a < s             #
###############################################################################

m7102 = [ [p1,         g1, pi/2 + a/2 - s/2, pi/2            ],
          [p2,         g1, pi/2 - s/2,       pi/2 + a/2 - s/2],
          [r*sin(a/2), g2, 0,                s - pi/2 + a/2  ] ]

p7102 = pi**-1 * (integrate(m7102[0][0], m7102[0][1:]) + integrate(m7102[1][0], m7102[1][1:]) + integrate(m7102[2][0], m7102[2][1:])).simplify()


rep7102 = {s:6*pi/8, a:4*pi/8} # Replacement values in range

# Confirm replacements
if not (a.subs(rep7102) < pi and pi/2 <= s.subs(rep7102) <= pi and a.subs(rep7102)/2 <= s.subs(rep7102)/2):
	print('rep7102 incorrect')

# is average profile in range 0r-2r?
if not 0 < p7102.subs(dict(rep7102, **r1)) < 2:
	print('Total p7102 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m7102)):
	if not integrate(m7102[i][0], m7102[i][1:]).subs(dict(rep7102, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p7102 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m7102)):
	if not 0 < (integrate(m7102[i][0], m7102[i][1:])/(m7102[i][3]-m7102[i][2])).subs(dict(rep7102, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p7102 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m7102)):
	if not (m7102[i][3]-m7102[i][2]).subs(rep7102) > 0:
		print('Bounds ' + str(i+1) + ' in p7102 has lower bounds bigger than upper bounds')	

# LaTeX output

latexOutput.insert(8,'p7102 = ' + latex(p7102))






###############################################################################
# 7.9.1 animal: a < pi.  Sensor: s < pi/2. Condition: a < pi - 2s & a < s     #
###############################################################################

m791 = [ [p1, g1, pi/2 - s/2 + a/2, pi/2            ],
         [p2, g1, pi/2 - s/2,       pi/2 - s/2 + a/2],
         [p3, g2, s,                s + a/2         ] ]

p791 = pi**-1 * (integrate(m791[0][0], m791[0][1:]) + integrate(m791[1][0], m791[1][1:]) + integrate(m791[2][0], m791[2][1:])).simplify().trigsimp()


rep791 = {s:2*pi/8, a:pi/8} # Replacement values in range

# Confirm replacements
if not (a.subs(rep791) < pi and s.subs(rep791) <= pi/2 and a.subs(rep791) <= pi - 2*s.subs(rep791) and a.subs(rep791) <= s.subs(rep791)):
	print('rep791 incorrect')

# is average profile in range 0r-2r?
if not 0 < p791.subs(dict(rep791, **r1)) < 2:
	print('Total p791 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m791)):
	if not integrate(m791[i][0], m791[i][1:]).subs(dict(rep791, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p791 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m791)):
	if not 0 < (integrate(m791[i][0], m791[i][1:])/(m791[i][3]-m791[i][2])).subs(dict(rep791, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p791 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m791)):
	if not (m791[i][3]-m791[i][2]).subs(rep791) > 0:
		print('Bounds ' + str(i+1) + ' in p791 has lower bounds bigger than upper bounds')	

# LaTeX output

latexOutput.insert(9,'p791 = ' + latex(p791))

##################################################################################
# 7.9.2 animal: a < pi.  Sensor: s < pi/2. Condition: a < pi - 2s & s < a < 2s   #
##################################################################################

m792 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2 + a/2, pi/2            ],
         [p2,                   g1, pi/2 - s/2,       pi/2 - s/2 + a/2],
         [p3,                   g2, s,                s + a/2         ] ]

p792 = pi**-1 * (integrate(m792[0][0], m792[0][1:]) + integrate(m792[1][0], m792[1][1:]) + integrate(m792[2][0], m792[2][1:])).simplify().trigsimp()


rep792 = {s:2*pi/8, a:pi/2-0.1} # Replacement values in range

# Confirm replacements
if not (a.subs(rep792) < pi and s.subs(rep792) <= pi/2 and a.subs(rep792) <= pi - 2*s.subs(rep792) and s.subs(rep792) <= a.subs(rep792) <= 2*s.subs(rep792)):
	print('rep792 incorrect')

# is average profile in range 0r-2r?
if not 0 < p792.subs(dict(rep792, **r1)) < 2:
	print('Total p792 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m792)):
	if not integrate(m792[i][0], m792[i][1:]).subs(dict(rep792, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p792 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m792)):
	if not 0 < (integrate(m792[i][0], m792[i][1:])/(m792[i][3]-m792[i][2])).subs(dict(rep792, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p792 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m792)):
	if not (m792[i][3]-m792[i][2]).subs(rep792) > 0:
		print('Bounds ' + str(i+1) + ' in p792 has lower bounds bigger than upper bounds')	

# LaTeX output

latexOutput.insert(10,'p792 = ' + latex(p792))

##################################################################################
# 7.9.3 animal: a < pi.  Sensor: s < pi/2. Condition: a < pi - 2s &  2s < a      #
##################################################################################



p793 = (1/pi) * ( integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) )       \
                + integrate( r*sin(g2), (g2, s, a/2) ) \
                + integrate( p3, (g2, a/2, s + a/2) ) ).simplify().trigsimp()

rep793 = {s:1*pi/8, a:pi/2} # Replacement values in range


m793 = [ [p1, g1, pi/2 - s/2 + a/2, pi/2            ],
         [p2, g1, pi/2 - s/2,       pi/2 - s/2 + a/2],
         [p3, g2, s,                s + a/2         ] ]

p793 = pi**-1 * (integrate(m793[0][0], m793[0][1:]) + integrate(m793[1][0], m793[1][1:]) + integrate(m793[2][0], m793[2][1:])).simplify().trigsimp()


rep793 = {s:2*pi/8, a:pi/8} # Replacement values in range

# Confirm replacements
if not (a.subs(rep793) < pi and s.subs(rep793) <= pi/2 and a.subs(rep793) <= pi - 2*s.subs(rep793) and a.subs(rep793) < s.subs(rep793)):
	print('rep793 incorrect')

# is average profile in range 0r-2r?
if not 0 < p793.subs(dict(rep793, **r1)) < 2:
	print('Total p793 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m793)):
	if not integrate(m793[i][0], m793[i][1:]).subs(dict(rep793, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p793 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m793)):
	if not 0 < (integrate(m793[i][0], m793[i][1:])/(m793[i][3]-m793[i][2])).subs(dict(rep793, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p793 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m793)):
	if not (m793[i][3]-m793[i][2]).subs(rep793) > 0:
		print('Bounds ' + str(i+1) + ' in p793 has lower bounds bigger than upper bounds')	

# LaTeX output

latexOutput.insert(11,'p793 = ' + latex(p793))


##################################################################################
# 7.9.4 animal: a < pi.  Sensor: s < pi/2. Condition: a > pi - 2s &  a < s       #
##################################################################################

m794 = [ [p1, g1, pi/2 - s/2 + a/2, pi/2            ],
         [p2, g1, pi/2 - s/2,       pi/2 - s/2 + a/2],
         [p3, g2, s,                pi/2            ] ]

p794 = pi**-1 * (integrate(m794[0][0], m794[0][1:]) + integrate(m794[1][0], m794[1][1:]) + integrate(m794[2][0], m794[2][1:])).simplify().trigsimp()

rep794 = {s:pi/2-0.1, a:pi/4} # Replacement values in range

# Confirm replacements
if not (a.subs(rep794) < pi and s.subs(rep794) <= pi/2 and a.subs(rep794) >= pi - 2*s.subs(rep794) and a.subs(rep794) < s.subs(rep794)):
	print('rep794 incorrect')

# is average profile in range 0r-2r?
if not 0 < p794.subs(dict(rep794, **r1)) < 2:
	print('Total p794 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m794)):
	if not integrate(m794[i][0], m794[i][1:]).subs(dict(rep794, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p794 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m794)):
	if not 0 < (integrate(m794[i][0], m794[i][1:])/(m794[i][3]-m794[i][2])).subs(dict(rep794, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p794 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m794)):
	if not (m794[i][3]-m794[i][2]).subs(rep794) > 0:
		print('Bounds ' + str(i+1) + ' in p794 has lower bounds bigger than upper bounds')	

# LaTeX output

latexOutput.insert(12,'p794 = ' + latex(p794))


##################################################################################
# 7.9.5 animal: a < pi.  Sensor: s < pi/2. Condition: a > pi - 2s &  s < a < 2s  #
##################################################################################


m795 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 + s/2 - a/2, pi/2            ],
         [p2,                   g1, pi/2 - s/2,       pi/2 + s/2 - a/2],
         [p3,                   g2, s/2,               pi/2            ] ]

p795 = pi**-1 * (integrate(m795[0][0], m795[0][1:]) + integrate(m795[1][0], m795[1][1:]) + integrate(m795[2][0], m795[2][1:])).simplify().trigsimp()


rep795 = {s:pi/2-0.1, a:pi/2} # Replacement values in range

# Confirm replacements
if not (a.subs(rep795) <= pi and s.subs(rep795) <= pi/2 and a.subs(rep795) >= pi - 2*s.subs(rep795) and s.subs(rep795) <= a.subs(rep795) <= 2*s.subs(rep795)):
	print('rep795 incorrect')

# is average profile in range 0r-2r?
if not 0 < p795.subs(dict(rep795, **r1)) < 2:
	print('Total p795 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m795)):
	if not integrate(m795[i][0], m795[i][1:]).subs(dict(rep795, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p795 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m795)):
	if not 0 < (integrate(m795[i][0], m795[i][1:])/(m795[i][3]-m795[i][2])).subs(dict(rep795, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p795 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m795)):
	if not (m795[i][3]-m795[i][2]).subs(rep795) > 0:
		print('Bounds ' + str(i+1) + ' in p795 has lower bounds bigger than upper bounds')	

# LaTeX output

latexOutput.insert(13,'p795 = ' + latex(p795))






##################################################################################
# 7.9.6 animal: a < pi.  Sensor: s < pi/2. Condition: a > pi - 2s &  a < 2s      #
##################################################################################



m796 = [ [2*r*sin(s/2)*sin(g1), g1, pi/2 - s/2, pi/2],
         [r*sin(g2),            g2, s,          a/2 ],
         [p3,                   g3, a/2,        pi/2] ]

p796 = pi**-1 * (integrate(m796[0][0], m796[0][1:]) + integrate(m796[1][0], m796[1][1:]) + integrate(m796[2][0], m796[2][1:])).simplify()


rep796 = {s:pi/4, a:3*pi/4} # Replacement values in range

# Confirm replacements
if not (a.subs(rep796) < pi and s.subs(rep796) <= pi/2 and a.subs(rep796) <= pi - 2*s.subs(rep796) and a.subs(rep796) < 2*s.subs(rep796)):
	print('rep796 incorrect')

# is average profile in range 0r-2r?
if not 0 < p796.subs(dict(rep796, **r1)) < 2:
	print('Total p796 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m796)):
	if not integrate(m796[i][0], m796[i][1:]).subs(dict(rep796, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p796 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m796)):
	if not 0 < (integrate(m796[i][0], m796[i][1:])/(m796[i][3]-m796[i][2])).subs(dict(rep796, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p796 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m796)):
	if not (m796[i][3]-m796[i][2]).subs(rep796) > 0:
		print('Bounds ' + str(i+1) + ' in p796 has lower bounds bigger than upper bounds')	

# LaTeX output

latexOutput.insert(14,'p796 = ' + latex(p796))




#############################################################################
# 7.6.1 animal: a < pi.  Sensor: s > pi. Condition: a < s - pi              #
#############################################################################

m761 = [ [2*r*sin(a/2), g4, pi/2, pi        ],
         [r*sin(a/2),   g4, pi,   2*pi - s/2] ]

p761 = pi**-1 * (integrate(m761[0][0], m761[0][1:]) + integrate(m761[1][0], m761[1][1:])).simplify()


rep761 = {s:3*pi/2, a:pi/3} # Replacement values in range

# Confirm replacements
if not (a.subs(rep761) < pi and s.subs(rep761) >= pi/2 and a.subs(rep761) <= s.subs(rep761) - pi):
	print('rep761 incorrect')

# is average profile in range 0r-2r?
if not 0 < p761.subs(dict(rep761, **r1)) < 2:
	print('Total p761 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m761)):
	if not integrate(m761[i][0], m761[i][1:]).subs(dict(rep761, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p761 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m761)):
	if not 0 < (integrate(m761[i][0], m761[i][1:])/(m761[i][3]-m761[i][2])).subs(dict(rep761, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p761 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m761)):
	if not (m761[i][3]-m761[i][2]).subs(rep761) > 0:
		print('Bounds ' + str(i+1) + ' in p761 has lower bounds bigger than upper bounds')	

# LaTeX output

latexOutput.insert(15,'p761 = ' + latex(p761))




#########################################################################
# 7.6.2 animal: a < pi.  Sensor: s < pi/2. Condition: a > s - pi        #
#########################################################################


array762 = [ [ 2*r*sin(a/2),                      g4, pi/2,         s/2       ],
	     [ r*sin(a/2) - r*cos(g4),            g4, s/2,          pi        ], 
	     [ r*sin(a/2),                        g4, pi,           2*pi - s/2],
	     [ 2*r*(sin(a/2) + sin(s/2)*sin(g4)), g4, 2*pi - s/2,   3*pi/2    ]  ]

p762 = pi**-1 * (2*integrate(m762[0][0], m762[0][1:]) + integrate(m762[1][0], m762[1][1:]) + integrate(m762[2][0], m762[2][1:]) + 2*integrate(m762[3][0], m762[3][1:])).simplify()


rep762 = {s:3*pi/2, a:pi/3} # Replacement values in range

# Confirm replacements
if not (a.subs(rep762) < pi and s.subs(rep762) <= pi/2 and a.subs(rep762) >= s.subs(rep762) - pi):
	print('rep762 incorrect')

# is average profile in range 0r-2r?
if not 0 < p762.subs(dict(rep762, **r1)) < 2:
	print('Total p762 not in 0, 2r')

# Are the individual integrals >0r
for i in range(len(m762)):
	if not integrate(m762[i][0], m762[i][1:]).subs(dict(rep762, **r1)) > 0:
		print('Integral ' + str(i+1) + ' in p762 is negative')

# Are the individual averaged integrals between 0 and  2r
for i in range(len(m762)):
	if not 0 < (integrate(m762[i][0], m762[i][1:])/(m762[i][3]-m762[i][2])).subs(dict(rep762, **r1)) <= 2:
		print('Integral ' + str(i+1) + ' in p762 has averaged integral outside 0<p<2r')
		
# Are the bounds the correct way around
for i in range(len(m762)):
	if not (m762[i][3]-m762[i][2]).subs(rep762) > 0:
		print('Bounds ' + str(i+1) + ' in p762 has lower bounds bigger than upper bounds')	

# LaTeX output

latexOutput.insert(16,'p762 = ' + latex(p762))







