import pylab
import numpy as np
from numpy import *
import math
from matplotlib import *
from pylab import *
import subprocess
from subprocess import call

# Getting the data from a file
# m = genfromtxt("starcov04.txt")
m = genfromtxt("bezier_curves_regressions.dat")
# Name of the output file
# title = "site-cov04a.eps"
title = "regressions_two_sites_2.eps"
# Boundingbox - deffines the size of the figure in points
BBLL = array([-10, 0])  # lower left corner
BBUR = array([433, 140])  # upper right corner
# Opening the file and writing the eps header
file = open(title, "w")
file.write("%!PS-Adobe-3.0 EPSF-1.2\n")
file.write("%%BoundingBox: {} {} {} {}\n\n".format(BBLL[0], BBLL[1], BBUR[0], BBUR[1]))
file.write("/Times-Roman findfont 6 scalefont setfont\n\n")
file.write("1 setlinejoin\n\n")
# Positioning the origin on the page
file.write("25 20 translate\n\n")
# Drawing the line representing the sites

# Drawing the tic marks for sites and nucleotides within each site
for i in range(40):
    locx = 10 * i
    add = 2
    if i + 1 < 10:
        add = 4
    file.write("{} -16 moveto\n".format(locx + add))
    file.write("({}) show\n".format(i + 1))
    file.write("{} -1 moveto\n".format(locx))
    file.write("0 -5 rlineto\n")
    file.write("stroke\n")
    file.write("gsave\n")
    file.write("/Times-Roman findfont 3.5 scalefont setfont\n")
    for j in range(4):
        locx = 10 * i + 2 * (j + 1)
        if j == 0:
            file.write("0 1 1 0 setcmykcolor\n") #following color scheme from paper
            file.write("{} -10 moveto\n".format(locx - 1.5))
            file.write("(A) show\n")
        if j == 1:
            file.write("1 1 0 0 setcmykcolor\n")
            file.write("{} -8 moveto\n".format(locx - 1.3))
            file.write("(C) show\n")
        if j == 2:
            file.write("1 0 1 0 setcmykcolor\n")
            file.write("{} -10 moveto\n".format(locx - 1.4))
            file.write("(G) show\n")
        if j == 3:
            file.write("0 0 1 0.3 setcmykcolor\n")
            file.write("{} -8 moveto\n".format(locx - 1.2))
            file.write("(U) show\n") #change to U for RNA
        file.write("{} -1 moveto\n".format(locx))
        file.write("0 -3 rlineto\n")
        file.write("stroke\n")

    file.write("grestore\n")

file.write("0 -1 moveto\n")
file.write("400 -1 lineto\n")
file.write("stroke\n\n")

file.write("-30 -9 moveto\n")
file.write("(Nucleotide) show\n")

file.write("-20 -16 moveto\n")
file.write("(Site) show\n")

file.write("gsave\n")
file.write("/Times-Roman findfont 10 scalefont setfont\n")
file.write("-20 100 moveto\n")
file.write("(Regressions of phenotypes onto site pairs) show\n")
file.write("-20 90 moveto\n")
file.write("(at different sites) show\n")

file.write("1 1 0 0 setcmykcolor\n")
file.write("-20 80 moveto\n")
file.write("30 0 rlineto\n")
file.write("stroke\n")
file.write("15 77.5 moveto\n")
file.write("(Positive) show\n")

file.write("0 0 1 .3 setcmykcolor\n")
file.write("-20 70 moveto\n")
file.write("30 0 rlineto\n")
file.write("stroke\n")
file.write("15 67 moveto\n")
file.write("(Negative) show\n")

file.write("grestore\n")

# Drawing the curved lines connecting nucleotides at sites.
# There were 32 pairs in this case
# 103 pairs for case with abs(cov)>.04
# change range when doing abs(cov)>.05
for i in range(14):
    llocx = 10. * m[i][0] + 2 * (m[i][2] + 1)
    llocy = 0
    rlocx = 10. * m[i][1] + 2 * (m[i][3] + 1)
    rlocy = 0
    d = rlocx - llocx
    file.write("gsave\n")
    #file.write("{} 0 translate\n\n".format(llocx))
    if m[i][4] > 0:
        file.write("1 1 0 0 setcmykcolor\n")
        file.write("{} 0 translate\n\n".format(llocx))
        file.write("0 0 moveto\n")
        file.write("{} {}\n".format(d / 3, d / 3))
        file.write("{} {}\n".format(2 * d / 3, d / 3))
        file.write("{} 0 curveto\n".format(d))
        file.write("stroke\n\n")
        file.write("grestore\n\n")
    if m[i][4] < 0:
        file.write("0 0 1 .3 setcmykcolor\n")
        # The following lines define Bezier curves connecting the sites/nucleotides
        file.write("{} 0 translate\n\n".format(llocx))
        file.write("0 0 moveto\n")
        file.write("{} {}\n".format(d / 3, d / 3))
        file.write("{} {}\n".format(2 * d / 3, d / 3))
        file.write("{} 0 curveto\n".format(d))
        file.write("stroke\n\n")
        file.write("grestore\n\n")

file.write("showpage\n")
file.close()
# This line converts the eps file to a pdf file.
# It might work only on Linux
call(["epstopdf", title])  # need to figure out how to do this part in windows, works in Linux
