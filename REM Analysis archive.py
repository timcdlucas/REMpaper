"""
Systematic analysis of REM models
Tim Lucas 
01/10/13
"""


from sympy import *
from numpy import arange
import os
import pylab as pl


%load_ext sympyprinting


# Set working directory
os.chdir('/home/tim/Dropbox/phd/Analysis/REM-chapter')


# Load symbols used for symbolic maths
s, a, r, g1, g2, g3, g4 = symbols('theta_s theta_a r gamma_1 gamma_2 gamma_3 gamma_4')


# A list that will recieve the latex ouput of each model

latexOutput = []

###################################################
# 7.3. animal: a > pi.  sensor: s > pi            #
###################################################

p73 = (1/pi) * ( 2*integrate(2*r, (g4, pi/2, s/2)) +             \
		 2*integrate(r + r*cos(g4 - s/2), (g4, s/2, pi)) \
                 ).simplify()
 

rep73 = {s:3*pi/2} # Replacement values in range
# Confirm replacements
s.subs(rep73) > pi 

# is average profile in range 0r-2r?
p73.subs(rep73).n()

# Are the individual integrals >0r
integrate(2*r , (g4, pi/2, s/2)).subs(rep73).n()
integrate(r+r*cos(g4-s/2), (g4, s/2, pi)).subs(rep73).n()

# Plot function

xRange = arange(pi,2*pi, 0.01)
yRange = [p73.subs({r:1, s:i}).n() for i in xRange]
plot73 = plot(xRange, yRange)
savefig('imgs/p73Profile.pdf')

# LaTeX output

latexOutput.insert(0,'p75 = ' + latex(p73))

# Start testing against other functions

p73.subs({s:pi}).together()



########################################################
# 7.5. animal: a = 2*pi.   sensor:  pi/2 < s < pi      #
########################################################

p75 = (1/pi) * ( 2*integrate( 2*r*sin(s/2)*sin(g1) , (g1, s/2, pi/2) )   \
               + 2*integrate( r - r*cos(g3 - s), (g3, 0, s - pi/2) )     \
               +   integrate( r, (g3, s - pi/2, pi/2) )                  \
                 ).trigsimp().simplify()



rep75 = {s:3*pi/4} # Replacement values in range
# Confirm replacements
pi/2 < s.subs(rep75) < pi

# is average profile in range 0r-2r?
p75.subs(rep75).n() 

# Are the individuals integrals >0r
integrate( 2*r*sin(s/2)*sin(g1) , (g1, s/2, pi/2) ).subs(rep75).n()
integrate( r - r*cos(g3 - s), (g3, 0, s - pi/2) ).subs(rep75).n()
integrate( r, (g3, s - pi/2, pi/2) ).subs(rep75).n()


# Plot function

xRange = arange(pi/2,pi, 0.01)
yRange = [p75.subs({r:1, s:i}).n() for i in xRange]
plot73 = plot(xRange, yRange)
savefig('imgs/p75Profile.pdf')

# LaTeX output

latexOutput.insert(1,'p75 = ' + latex(p75))



"""
Multiple models in 7.7.
animal: a>pi.  Sensor: s < pi/2
"""


###################################################################
# 7.7.1 animal: a>pi.  Sensor: s < pi/2. Condition: a/2 < pi - s  #
###################################################################

p771 = (1/pi) * ( integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) ) \
                + integrate( r*sin(g2), (g2, s, pi/2) )                     \
                + integrate( r*sin(g2), (g2, s, pi - s/2) )                 \
                + s*r ).simplify().trigsimp()

rep771 = {s:pi/9, a:10*pi/9} # Replacement values in range
# Check replacements values
a.subs(rep771) > pi and s.subs(rep771) < pi/2 and a.subs(rep771)/2 < pi - s.subs(rep771)


# is average profile in range 0r-2r?
p771.subs(rep771).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) ).subs(rep771).n()
integrate( r*sin(g2), (g2, s, pi/2) ).subs(rep771).n()
integrate( r*sin(g2), (g2, s, pi - s/2) ).subs(rep771).n()

# Check sizes of integral bounds
pi/2 - s.subs(rep771)/2 < pi/2            and \
s.subs(rep771) < pi/2                     and \
s.subs(rep771) < pi - s.subs(rep771)/2

# LaTeX output

latexOutput.insert(2,'p771 = ' + latex(p771))



###############################################################################
# 7.7.2 animal: a>pi.  Sensor: s < pi/2. Condition:  pi - s < a/2 < pi - s/2  #
###############################################################################

p772 = (1/pi) * ( integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) )        \
              + 2*integrate( r*sin(g2), (g2, s, pi/2) )                            \
              +   integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2, a/2 + s/2 - pi/2) ) \
              +   s*r ).simplify().trigsimp()

rep772 = {s:3*pi/8, a:3*pi/2} # Replacement values in range

# Confirm our replacement rules are correct.
a.subs(rep772)>pi and s.subs(rep772)<pi/2 and pi-s.subs(rep772) < a.subs(rep772)/2 < pi - s.subs(rep772)/2



# is average profile in range 0r-2r?
p772.subs(rep772).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) ).subs(rep772).n()
integrate( r*sin(g2), (g2, s, pi/2) ).subs(rep772).n()
integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2, a/2 + s/2 - pi/2) ).subs(rep772).n()

# Check sizes of integral bounds
pi/2 - s.subs(rep772)/2 < pi/2            and \
s.subs(rep772) < pi/2                     and \
pi/2 - s.subs(rep772) < s.subs(rep772)/2 + a.subs(rep772) - pi/2

# LaTeX output

latexOutput.insert(2,'p772 = ' + latex(p772))


###############################################################################
# 7.7.3 animal: a>pi.  Sensor: s < pi/2. Condition: pi - s/2 < a/2            #
###############################################################################

p773 = (1/pi) * ( integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) )          \
              + 2*integrate( r*sin(g2), (g2, s, pi/2) )                              \
              +   integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2, 3*pi/2 - a/2 + s/2) ) \
              +   s*r ).simplify().trigsimp()

rep773 = {s:3*pi/8, a:29*pi/16} # Replacement values in range
# Confirm replacements
a.subs(rep773) > pi and s.subs(rep773) < pi/2 and  pi - s.subs(rep773)/2 < a.subs(rep773)/2 


# is average profile in range 0r-2r?
p773.subs(rep773).n()

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

latexOutput.insert(3,'p773 = ' + latex(p773))








"""
Two models in 7.8.
animal: pi < a .  Sensor: pi/2 < s < pi
"""


#################################################################################
# 7.8.1 animal:  a > pi.  Sensor: pi/2 < s < pi. Condition: a/2 < pi - s/2      #
#################################################################################

p781 = (1/pi) * ( integrate( 2*r*sin(a/2), (g1, pi/2 - s/2, pi/2) )                \
                + integrate( r*sin(g2), (g2, s/2, pi/2) )                          \
                + integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2, a/2 + s/2 -pi/2) )  \
                + s*r/2 ).simplify().trigsimp()

rep781 = {s:3*pi/4, a:15*pi/8} # Replacement values in range
# Check replacements values
a.subs(rep781) > pi and pi/2 < s.subs(rep781) < pi and a.subs(rep781)/2 < pi - s.subs(rep781)/2


# is average profile in range 0r-2r?
p781.subs(rep781).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(a/2), (g1, pi/2 - s/2, pi/2) ).subs(rep781).n()
integrate( r*sin(g2), (g2, s/2, pi/2) ).subs(rep781).n()
integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2, a/2 + s/2 -pi/2) ).subs(rep781).n()

# Check size of integral bounds

# LaTeX output

latexOutput.insert(3,'p781 = ' + latex(p781))



#################################################################################
# 7.8.2 animal: a > pi.  Sensor: pi/2 < s < pi. Condition: a/2 > pi - s/2  #
#################################################################################

p782 = (1/pi) * ( integrate( 2*r*sin(a/2), (g1, pi/2 - s/2, pi/2) )                   \
                + integrate( r*sin(g2), (g2, s/2, pi/2) )                             \
                + integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2, 3*pi/2 + s/2 - a/2) )  \
                + integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2, 3*pi/2 + s/2 - a/2) )   \
                + s*r/2 ).simplify().trigsimp()

rep782 = {s:3*pi/4, a:15*pi/8} # TODO Replacement values in range
# Check replacements values
a.subs(rep782) > pi and pi/2 < s.subs(rep782) < pi and a.subs(rep782)/2 < pi - s.subs(rep782)/2


# is average profile in range 0r-2r?
p782.subs(rep782).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(a/2), (g1, pi/2 - s/2, pi/2) ) .subs(rep782).n()
integrate( r*sin(g2), (g2, s/2, pi/2) ).subs(rep782).n()
integrate( r*cos(g1 - s/2), (g1, pi/2 - s/2,  3*pi/2 + s/2 - a/2) ).subs(rep782).n()
integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2, 3*pi/2 + s/2 - a/2) ).subs(rep782).n()

# Check size of integral bounds

# LaTeX output

latexOutput.insert(4,'p782 = ' + latex(p782))





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

p7101 = (1/pi) * ( integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - a/2 + s/2, pi/2) )     \
                 + integrate( p2, (g1, pi/2 - s/2, pi/2 -a/2 +s/2) )                   \
                 + integrate( r*sin(a/2), (g2, 0, s - pi/2 + a/2) ) ).simplify().trigsimp()

rep7101 = {s:5*pi/8, a:7*pi/8} # Replacement values in range
# Confirm replacements
a.subs(rep7101) < pi and pi/2 < s.subs(rep7101) < pi and a.subs(rep7101) > s.subs(rep7101)


# is average profile in range 0r-2r?
p7101.subs(rep7101).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - a/2 + s/2, pi/2) ).subs(rep7101).n()
integrate( p2, (g1, pi/2 - s/2, pi/2 -a/2 +s/2) ).subs(rep7101).n()
integrate( r*sin(a/2), (g2, 0, s - pi/2 + a/2) ) ).subs(rep7101).n()

# Check sizes of integral bounds 
pi/2 - a.subs(rep7101)/2 + s.subs(rep7101)/2 < pi/2 and \
pi/2 - s.subs(rep7101)/2 <  pi/2 - a.subs(rep7101)/2 + s.subs(rep7101)/2 and \
0 < s.subs(rep7101) - pi/2 + a.subs(rep(7101)/2


# LaTeX output

latexOutput.insert(5,'p7101 = ' + latex(p7101))





###############################################################################
# 7.10.2 animal: a < pi.  Sensor: pi/2 < s < pi. Condition: a < s             #
###############################################################################

p7102 = (1/pi) * ( integrate( p1, (g1, pi/2 - a/2 + s/2, pi/2) )     \
                 + integrate( p2, (g1, pi/2 - s/2, pi/2 -a/2 +s/2) )                   \
                 + integrate( r*sin(a/2), (g2, 0, s - pi/2 + a/2) ) ).simplify().trigsimp()

rep7102 = {s:6*pi/8, a:4*pi/8} # Replacement values in range
# Confirm replacements
a.subs(rep7102) < pi and pi/2 < s.subs(rep7102) < pi and a.subs(rep7102) < s.subs(rep7102)


# is average profile in range 0r-2r?
p7102.subs(rep7102).n()

# Are the individuals integrals >0r
integrate( p, (g1, pi/2 - a/2 + s/2, pi/2) ).subs(rep7102).n()
integrate( p2, (g1, pi/2 - s/2, pi/2 -a/2 +s/2) ).subs(rep7102).n()
integrate( r*sin(a/2), (g2, 0, s - pi/2 + a/2) ) ).subs(rep7102).n()

# Check sizes of integral bounds 
pi/2 - a.subs(rep7102)/2 + s.subs(rep7102)/2 < pi/2 and \
pi/2 - s.subs(rep7102)/2 <  pi/2 - a.subs(rep7102)/2 + s.subs(rep7102)/2 and \
0 < s.subs(rep7102) - pi/2 + a.subs(rep7102)/2

# LaTeX output

latexOutput.insert(6,'p7102 = ' + latex(p7102))








###############################################################################
# 7.9.1 animal: a < pi.  Sensor: s < pi/2. Condition: a < pi - 2s & a < s     #
###############################################################################

p791 = (1/pi) * ( integrate( p1, (g1, pi/2 - s/2 + a/2, pi/2) )       \
                + integrate( p2, (g1, pi/2 - s/2, pi/2 - s/2 + a/2) ) \
                + integrate( p3, (g2, s, s + a/2) ) ).simplify().trigsimp()

rep791 = {s:2*pi/8, a:pi/8} # Replacement values in range
# Confirm replacements
a.subs(rep791) < pi and  s.subs(rep791) < pi/2 and a.subs(rep791) < pi - 2*s.subs(rep791) and a.subs(rep791) < s.subs(rep791)


# is average profile in range 0r-2r?
p791.subs(rep791).n()

# Are the individuals integrals >0r
integrate( p1, (g1, pi/2 - s/2 + a/2, pi/2) ).subs(rep791).n()
integrate( p2, (g1, pi/2 - s/2, pi/2 - s/2 + a/2) ).subs(rep791).n()
integrate( p3, (g2, s, s + a/2) ).subs(rep791).n()

# Check sizes of integral bounds 
pi/2 - s.subs(rep791)/2 + a.subs(rep791)/2 < pi/2 and \
pi/2 - s.subs(rep791)/2 < pi/2 - s.subs(rep791)/2 + a.subs(rep791)/2 and \
s.subs(rep791) < s.subs(rep791) + a.subs(rep791)/2

# LaTeX output

latexOutput.insert(7,'p791 = ' + latex(p791))


##################################################################################
# 7.9.2 animal: a < pi.  Sensor: s < pi/2. Condition: a < pi - 2s & s < a < 2s   #
##################################################################################

p792 = (1/pi) * ( integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2 + a/2, pi/2) )       \
                + integrate( p2, (g1, pi/2 - s/2, pi/2 - s/2 + a/2) ) \
                + integrate( p3, (g2, s, s + a/2) ) ).simplify().trigsimp()

rep792 = {s:2*pi/8, a:pi/2-0.1} # Replacement values in range
# Confirm replacements
a.subs(rep792) < pi and  s.subs(rep792) < pi/2 and a.subs(rep791) < pi - 2*s.subs(rep791) and s.subs(rep792) < a.subs(rep792) < 2*s.subs(rep792)


# is average profile in range 0r-2r?
p792.subs(rep792).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2 + a/2, pi/2) ).subs(rep792).n()
integrate( p2, (g1, pi/2 - s/2, pi/2 - s/2 + a/2) ).subs(rep792).n()
integrate( p3, (g2, s, s + a/2) ).subs(rep792).n()

# Check sizes of integral bounds 
pi/2 - s.subs(rep792)/2 + a.subs(rep792)/2 < pi/2 and \
pi/2 - s.subs(rep792)/2 < pi/2 - s.subs(rep792)/2 + a.subs(rep792)/2 and \
s.subs(rep792) < s.subs(rep792) + a.subs(rep792)/2

# LaTeX output

latexOutput.insert(8,'p792 = ' + latex(p792))




##################################################################################
# 7.9.3 animal: a < pi.  Sensor: s < pi/2. Condition: a < pi - 2s &  2s < a      #
##################################################################################



p793 = (1/pi) * ( integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) )       \
                + integrate( r*sin(g2), (g2, s, a/2) ) \
                + integrate( p3, (g2, a/2, s + a/2) ) ).simplify().trigsimp()

rep793 = {s:1*pi/8, a:pi/2} # Replacement values in range
# Confirm replacements
a.subs(rep793) < pi and  s.subs(rep793) < pi/2 and a.subs(rep791) < pi - 2*s.subs(rep791) and 2*s.subs(rep793) < a.subs(rep793)


# is average profile in range 0r-2r?
p793.subs(rep793).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) ).subs(rep793).n()
integrate( r*sin(g2), (g2, s, a/2) ).subs(rep793).n()
integrate( p3, (g2, a/2, s + a/2) ) ).n()

# Check sizes of integral bounds 
 pi/2 - s.subs(rep793)/2 < pi/2 and \
s.subs(rep793) < a.subs(rep793)/2 and \
a.subs(rep793)/2 < s.subs(rep793) + a.subs(rep793)/2

# LaTeX output

latexOutput.insert(9,'p793 = ' + latex(p793))





##################################################################################
# 7.9.4 animal: a < pi.  Sensor: s < pi/2. Condition: a > pi - 2s &  a < s       #
##################################################################################



p794 = (1/pi) * ( integrate( p1, (g1, pi/2 - s/2 + a/2, pi/2) )       \
                + integrate( p2, (g1, pi/2 - s/2, pi/2 - s/2 + a/2) ) \
                + integrate( p3, (g2, s, pi/2) ) ).simplify().trigsimp()

rep794 = {s:pi/2-0.1, a:pi/4} # Replacement values in range
# Confirm replacements
a.subs(rep794) < pi and  s.subs(rep794) < pi/2 and a.subs(rep791) > pi - 2*s.subs(rep791) and a.subs(rep794) < s.subs(rep794)


# is average profile in range 0r-2r?
p794.subs(rep794).n()

# Are the individuals integrals >0r
integrate( p1, (g1, pi/2 - s/2 + a/2, pi/2) ).subs(rep794).n()
integrate( p2, (g1, pi/2 - s/2, pi/2 - s/2 + a/2) ).subs(rep794).n()
integrate( p3, (g2, s, pi/2) ).n()

# Check sizes of integral bounds 
pi/2 - s.subs(rep794)/2 + a.subs(rep794)/2 < pi/2 and \
pi/2 - s.subs(rep794)/2 < pi/2 - s.subs(rep794)/2 + a.subs(rep794)/2 and \
s.subs(rep794) < pi/2

# LaTeX output

latexOutput.insert(10,'p794 = ' + latex(p794))





##################################################################################
# 7.9.5 animal: a < pi.  Sensor: s < pi/2. Condition: a > pi - 2s &  s < a < 2s  #
##################################################################################



p795 = (1/pi) * ( integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 + s/2 - a/2, pi/2) )       \
                + integrate( p2, (g1, pi/2 - s/2,  pi/2 + s/2 - a/2) ) \
                + integrate( p3, (g2, s/2, pi/2) ) ).simplify().trigsimp()

rep795 = {s:1*pi/8, a:pi/2} # Replacement values in range
# Confirm replacements
a.subs(rep795) < pi and  s.subs(rep795) < pi/2 and a.subs(rep791) > pi - 2*s.subs(rep791) and s.subs(rep795) < a.subs(rep795) < 2*s.subs(rep795) 


# is average profile in range 0r-2r?
p795.subs(rep795).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 + s/2 - a/2, pi/2) ).subs(rep795).n()
integrate( p2, (g1, pi/2 - s/2,  pi/2 + s/2 - a/2) ).subs(rep795).n()
integrate( p3, (g2, s/2, pi/2) ) ).n()

# Check sizes of integral bounds 
pi/2 + s.subs(rep795)/2 - a.subs(rep795)/2 < pi/2 and \
pi/2 - ssubs(rep795)/2 <a  pi/2 + ssubs(rep795)/2 - asubs(rep795)/2 and \
s.subs(rep795)/2 < pi/2

# LaTeX output

latexOutput.insert(11,'p795 = ' + latex(p795))








##################################################################################
# 7.9.6 animal: a < pi.  Sensor: s < pi/2. Condition: a > pi - 2s &  a < 2s      #
##################################################################################



p796 = (1/pi) * ( integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) )       \
                + integrate( r*sin(g2), (g2, s, a/2) ) \
                + integrate( p3, (g3, a/2, pi/2) ) ).simplify().trigsimp()

rep796 = {s:pi/4, a:3*pi/4} # Replacement values in range
# Confirm replacements
a.subs(rep796) < pi and  s.subs(rep796) < pi/2 and a.subs(rep791) > pi - 2*s.subs(rep791) and  a.subs(rep796) < 2*s.subs(rep796)


# is average profile in range 0r-2r?
p796.subs(rep796).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(s/2)*sin(g1), (g1, pi/2 - s/2, pi/2) ).subs(rep796).n()
integrate( r*sin(g2), (g2, s, a/2) ).subs(rep796).n()
integrate( p3, (g3, a/2, pi/2) ) ).subs(rep796).n()

# Check sizes of integral bounds 
pi/2 - s.subs(rep796)/2 < pi/2 and \
s.subs(rep796) < a.subs(rep796)/2 and \
a.subs(rep796)/2 < pi/2

# LaTeX output

latexOutput.insert(12,'p796 = ' + latex(p796))





#############################################################################
# 7.6.1 animal: a < pi.  Sensor: s > pi. Condition: a < s - pi              #
#############################################################################



p761 = (1/pi) * ( integrate( 2*r*sin(a/2), (g4, pi/2, pi) )       \
                + integrate( r*sin(a/2), (g4, pi, 2*pi - s/2) ) ).simplify().trigsimp()

rep761 = {s:3*pi/2, a:pi/3} # Replacement values in range
# Confirm replacements
a.subs(rep761) < pi and s.subs(rep761) > pi and a.subs(rep761) < s.subs(rep761) - pi


# is average profile in range 0r-2r?
p761.subs(rep761).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(a/2), (g4, pi/2, pi) ).subs(rep761).n()
integrate( r*sin(a/2), (g4, pi, 2*pi - s/2) ) ).subs(rep761).n()


# Check sizes of integral bounds 
pi/2 < pi and \
pi < 2*pi - s.subs(rep761)/2

# LaTeX output

latexOutput.insert(13,'p761 = ' + latex(p761))






#########################################################################
# 7.6.2 animal: a < pi.  Sensor: s < pi/2. Condition: a > s - pi        #
#########################################################################




p762 = (1/pi) * ( integrate( 2*r*sin(a/2), (g4, pi/2, s/2) )       \
                + integrate( r*sin(a/2) - r*cos(g4), (g4, s/2, pi) ) \
                + integrate( r*sin(a/2), (g4, pi, 2*pi - s/2) ) \
                + integrate( 2*r*(sin(a/2) + sin(d/2)*sin(g4) ), (g4, 2*pi - s/2, 3*pi/2) ) ).simplify().trigsimp()

rep762 = {s:3*pi/2, a:2*pi/3} # Replacement values in range
# Confirm replacements
a.subs(rep762) < pi and s.subs(rep762) > pi and a.subs(rep762) > s.subs(rep762) - pi


# is average profile in range 0r-2r?
p762.subs(rep762).n()

# Are the individuals integrals >0r
integrate( 2*r*sin(a/2), (g4, pi/2, s/2) ) .subs(rep762).n()
integrate( r*sin(a/2) - r*cos(g4), (g4, s/2, pi) ).subs(rep762).n()
integrate( r*sin(a/2), (g4, pi, 2*pi - s/2) ).subs(rep762).n()
integrate( 2*r*(sin(a/2) + sin(d/2)*sin(g4) ), (g4, 2*pi - s/2, 3*pi/2) ).subs(rep762).n()


# Check sizes of integral bounds 
pi/2 < s.subs(rep762)/2 and \
s.subs(rep762)/2 < pi and \
pi < 2*pi - s.subs(rep762)/2 and \
2*pi - s.subs(rep762)/2 < 3*pi/2

# LaTeX output

latexOutput.insert(14,'p762 = ' + latex(p762))




