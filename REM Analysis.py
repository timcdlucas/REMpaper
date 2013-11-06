"""
Systematic analysis of REM models
Tim Lucas 
01/10/13
"""

from sympy import *
from numpy import arange
import os


%load_ext sympyprinting


# Set working directory
os.chdir('/home/tim/Dropbox/phd/Analysis/REM chapter')


# Load symbols used for symbolic maths
s, a, r, g1, g2, g3, g4 = symbols('theta_s theta_a r gamma_1 gamma_2 gamma_3 gamma_4')


# A list that will recieve the latex ouput of each model

latexOutput = []

###################################################
# 7.3. animal: a > pi.  sensor: s > pi            #
###################################################

q73 = (1/pi) * ( 2*integrate(2*r, (g4, pi/2, s/2)) +             \
		 2*integrate(r + r*cos(g4 - s/2), (g4, s/2, pi)) \
               ).simplify()
 

rep73 = {s:3*pi/2} # Replacement values in range
# Confirm replacements
s.subs(rep73) > pi 

# is average profile in range 0r-2r?
q73.subs(rep73).n()

# Are the individual integrals >0r
integrate(2*r , (g4, pi/2, s/2)).subs(rep73).n()
integrate(r+r*cos(g4-s/2), (g4, s/2, pi)).subs(rep73)

# Plot function

xRange = arange(pi,2*pi, 0.01)
yRange = [q73.subs({r:1, s:i}).n() for i in xRange]
plot73 = plot(xRange, yRange)
savefig('imgs/q73Profile.pdf')

# LaTeX output

latexOutput.insert(0,'q75 = ' + latex(q73))

# Start testing against other functions

q73.subs({s:pi}).together()



######################################################
# 7.5. animal: a > pi.   sensor:  pi/2 < s < pi      #
######################################################

q75 = (1/pi) * ( 2*integrate( 2*r*sin(s/2)*sin(g1) , (g1, s/2, pi/2) )   \
               + 2*integrate( r - r*cos(g3 - s), (g3, 0, s - pi/2) )     \
               +   integrate( r, (g3, s - pi/2, pi/2) )                  \
               ).trigsimp().simplify()

q75 = q75.subs({2*sin(s/2)*cos(s/2):sin(s)}) # double angle formulae not implemented. There is a script.


rep75 = {s:3*pi/4} # Replacement values in range
# Confirm replacements
pi/2 < s.subs(rep75) < pi

# is average profile in range 0r-2r?
q75.subs(rep75).n() 

# Are the individuals integrals >0r
integrate( 2*r*sin(s/2)*sin(g1) , (g1, s/2, pi/2) ).subs(rep75).n()
integrate( r - r*cos(g3 - s), (g3, 0, s - pi/2) ).subs(rep75).n()
integrate( r, (g3, s - pi/2, pi/2) ).subs(rep75).n()


# Plot function

xRange = arange(pi/2,pi, 0.01)
yRange = [q75.subs({r:1, s:i}).n() for i in xRange]
plot73 = plot(xRange, yRange)
savefig('imgs/q75Profile.pdf')

# LaTeX output

latexOutput.insert(1,'q75 = ' + latex(q75))



"""
Multiple models in 7.7.
animal: a>pi.  Sensor: s < pi/2
"""


###################################################################
# 7.7.1 animal: a>pi.  Sensor: s < pi/2. Condition: a/2 < pi - s  #
###################################################################

q771 = (1/pi) * ( integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) ) \
                + integrate( r*sin(g2), (g2, s, pi/2) )                     \
                + integrate( r*sin(g2), (g2, s, pi - s/2) )                 \
                + s*r ).simplify().trigsimp()

rep771 = {s:pi/9, a:10*pi/9} # Replacement values in range
# Check replacements values
a.subs(rep771) > pi and s.subs(rep771) < pi/2 and a.subs(rep771)/2 < pi - s.subs(rep771)


# is average profile in range 0r-2r?
q771.subs(rep771).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) ).subs(rep771).n()
integrate( r*sin(g2), (g2, s, pi/2) ).subs(rep771).n()
integrate( r*sin(g2), (g2, s, pi - s/2) ).subs(rep771).n()

# Check sizes of integral bounds
pi/2 - s.subs(rep771)/2 < pi/2            and \
s.subs(rep771) < pi/2                     and \
s.subs(rep771) < pi - s.subs(rep771)/2

# LaTeX output

latexOutput.insert(2,'q771 = ' + latex(q771))



###############################################################################
# 7.7.2 animal: a>pi.  Sensor: s < pi/2. Condition:  pi - s < a/2 < pi - s/2  #
###############################################################################

q772 = (1/pi) * ( integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) )        \
              + 2*integrate( r*sin(g2), (g2, s, pi/2) )                            \
              +   integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2, a/2 + s/2 - pi/2) ) \
              +   s*r ).simplify().trigsimp()

rep772 = {s:3*pi/8, a:3*pi/2} # Replacement values in range

# Confirm our replacement rules are correct.
a.subs(rep772)>pi and s.subs(rep772)<pi/2 and pi-s.subs(rep772) < a.subs(rep772)/2 < pi - s.subs(rep772)/2



# is average profile in range 0r-2r?
q772.subs(rep772).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) ).subs(rep772).n()
integrate( r*sin(g2), (g2, s, pi/2) ).subs(rep772).n()
integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2, a/2 + s/2 - pi/2) ).subs(rep772).n()

# Check sizes of integral bounds
pi/2 - s.subs(rep772)/2 < pi/2            and \
s.subs(rep772) < pi/2                     and \
pi/2 - s.subs(rep772) < s.subs(rep772)/2 + a.subs(rep772) - pi/2

# LaTeX output

latexOutput.insert(2,'q772 = ' + latex(q772))


###############################################################################
# 7.7.3 animal: a>pi.  Sensor: s < pi/2. Condition: pi - s/2 < a/2            #
###############################################################################

q773 = (1/pi) * ( integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) )          \
              + 2*integrate( r*sin(g2), (g2, s, pi/2) )                              \
              +   integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2, 3*pi/2 - a/2 + s/2) ) \
              +   s*r ).simplify().trigsimp()

rep773 = {s:3*pi/8, a:29*pi/16} # Replacement values in range
# Confirm replacements
a.subs(rep773) > pi and s.subs(rep773) < pi/2 and  pi - s.subs(rep773)/2 < a.subs(rep773)/2 


# is average profile in range 0r-2r?
q773.subs(rep773).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) ).subs(rep773).n()
integrate( r*sin(g2), (g2, s, pi/2) ) .subs(rep773).n()
integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2, 3*pi/2 - a/2 + s/2) ).subs(rep773).n()

# Check sizes of integral bounds TODO IS FALSE
pi/2 - s.subs(rep773)/2 < pi/2                                         and \
s.subs(rep773) < pi/2                                                  and \
pi/2 - s.subs(rep773) < 3* pi/2 + s.subs(rep773)/2 - a.subs(rep773)    and \
pi/2 <  3* pi/2 + s.subs(rep773)/2 - a.subs(rep773) 


# LaTeX output

latexOutput.insert(3,'q773 = ' + latex(q773))








"""
Two models in 7.8.
animal: pi < a .  Sensor: pi/2 < s < pi
"""


#################################################################################
# 7.8.1 animal:  a > pi.  Sensor: pi/2 < s < pi. Condition: a/2 < pi - s/2      #
#################################################################################

q781 = (1/pi) * ( integrate( 2*r*sin(a/2), (g1, pi/2 - s/2, pi/2) )                \
                + integrate( r*sin(g2), (g2, s/2, pi/2) )                          \
                + integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2, a/2 + s/2 -pi/2) )  \
                + s*r/2 ).simplify().trigsimp()

rep781 = {s:3*pi/4, a:15*pi/8} # Replacement values in range
# Check replacements values
a.subs(rep781) > pi and pi/2 < s.subs(rep781) < pi and a.subs(rep781)/2 < pi - s.subs(rep781)/2


# is average profile in range 0r-2r?
q781.subs(rep781).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(a/2), (g1, pi/2 - s/2, pi/2) ).subs(rep781).n()
integrate( r*sin(g2), (g2, s/2, pi/2) ).subs(rep781).n()
integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2, a/2 + s/2 -pi/2) ).subs(rep781).n()

# Check size of integral bounds

# LaTeX output

latexOutput.insert(3,'q781 = ' + latex(q781))



#################################################################################
# 7.8.2 animal: a > pi.  Sensor: pi/2 < s < pi. Condition: a/2 > pi - s/2  #
#################################################################################

q782 = (1/pi) * ( integrate( 2*r*sin(a/2), (g1, pi/2 - s/2, pi/2) )                   \
                + integrate( r*sin(g2), (g2, s/2, pi/2) )                             \
                + integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2, 3*pi/2 + s/2 - a/2) )  \
                + integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2, 3*pi/2 + s/2 - a/2) )   \
                + s*r/2 ).simplify().trigsimp()

rep782 = {s:3*pi/4, a:15*pi/8} # TODO Replacement values in range
# Check replacements values
a.subs(rep782) > pi and pi/2 < s.subs(rep782) < pi and a.subs(rep782)/2 < pi - s.subs(rep782)/2


# is average profile in range 0r-2r?
q782.subs(rep782).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(a/2), (g1, pi/2 - s/2, pi/2) ) .subs(rep782).n()
integrate( r*sin(g2), (g2, s/2, pi/2) ).subs(rep782).n()
integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2,  3*pi/2 + s/2 - a/2) ).subs(rep782).n()
integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2, 3*pi/2 + s/2 - a/2) ).subs(rep782).n()

# Check size of integral bounds

# LaTeX output

latexOutput.insert(4,'q782 = ' + latex(q782))





"""
Complex profiles for a < pi/2 
"""

# p-l-r for g1 profil. Calculated by AE in fig 7.2 minus AE in fig 7.3
p1 = 2*r*sin(s - 3*pi/2 + g1)*sin((g1 - pi/2 - s + a)/2) -  \
     2*r*sin((pi - a - 2*g1 + s)/4)*sin((pi - a +2*g1 - s)/4).simplify()

# p-l for g1 profiles
p2 = 2*r*sin(s/2)*sin(g1) - 2*r*sin((pi - a - 2*g1 + s)/4)*sin((pi - a + 2*g1 - s)/4)

# p-l for g1 profile. 
p3 = r*sin(g2) - 2*r*sin(g2/2 - a/4)*sin((pi - a + 2*g2 - s)/4).simpligy()














