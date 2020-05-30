####THIS is an attempt to make a command line tool using click.

from numpy import *
import numpy as np
#import sr
import time
from inspect import getsourcefile
import os
import click

#filepath = os.path.join(os.path.dirname(os.path.abspath(getsourcefile(lambda: 0))), filename)


start_time = time.time()

@click.command(help='program to compute orthogonal polynomials up to 3rd order')
@click.option('--n', default=1, help='Population size or number of samples')
@click.option('--dm', default=4, help='dimension of vector, e.g., this is =4 when input is DNA/RNA')
@click.option('--sites', default=3, help='number of sites in a sequence')
#@click.option()
@click.argument('filename', type=click.File('rb'))
def orthogonal_polynomial(filename, phenotype, sites, dm, n):
    """Program to compute orthogonal polynomials up to 3rd order"""
    with open(filename) as f:
        seq = f.readlines()
    global i
    F = Fest = Fon1 = Fon2i1 = Fon12 = genfromtxt(phenotype)  # file containing trait values that will be mapped to sequence, # vectors that must be the same size as F
    for i in range(n):
        Fest[i] = 0
        Fon1[i] = 0
        Fon2i1[i] = 0
        Fon12[i] = 0

    # ----Initializing various terms that we will use.--------------
    # 3 sites, each a dm dim vector, in n individuals
    # nOTE: For application to Amino Acid sequences, increase
    # the size of the arrays accordingly.
    import sr
    phi = array([[[0.0 for k in range(dm)] for i in range(n)] for j in range(sites)])  # general enough for all sites
    mean = array([[0.0 for z in range(dm)] for i in range(sites)])
    var = array([[0.0 for z in range(dm)] for i in range(sites)])
    phi2 = array([[[[[0.0 for k in range(dm)] for i in range(dm)] for j in range(n)] for l in range(sites)] for m in range(sites)])
    phi2m = array([[[[0.0 for k in range(dm)] for i in range(dm)] for l in range(sites)] for m in range(sites)])
    phi3 = array([[[[0.0 for k in range(dm)] for i in range(dm)] for j in range(dm)] for l in range(n)])
    phi3m = array([[[0.0 for k in range(dm)] for i in range(dm)] for j in range(dm)])
    P = array([[[0.0 for z in range(dm)] for j in range(n)] for i in range(sites)])
    # CM = array([[[0.0 for k in range(dm)] for i in range(sites)] for m in range(0, M)])
    cov = array([[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(sites)] for k in range(sites)])
    # ------------Converting letters to vectors---------------
    # phi[individual][site][state]. phi[i][j] = vector for site j in individual i.
    for i in range(0, n):  # individual
        for j in range(sites):
            if seq[i][j] == 'A':
                phi[j][i][0] = 1.0
            if seq[i][j] == 'C':
                phi[j][i][1] = 1.0
            if seq[i][j] == 'G':
                phi[j][i][2] = 1.0
            if seq[i][j] == 'T':
                phi[j][i][3] = 1.0
    # keep in alpha order
    # ---------------------------------First order terms --------------------------------
    # calculate mean vectors #stays same
    for i in range(0, n):
        for j in range(sites):
            mean[j] += phi[j][i] / n
    naming = os.path.basename(f.name)
    np.save(naming + str('_mean'), mean)
    print("mean") #to show progress, can do something much more efficient/elegant

    for j in range(sites):  # site
        for i in range(0, n):  # indiv
            P[j][i] = phi[j][i] - mean[j]
    # var[site][nucleotide]
    for k in range(sites):
        for i in range(0, dm):  # nucleotide
            for j in range(0, n):  # individual
                var[k][i] += ((P[k][j][i]) ** 2) / n

    # # Covariances between nucleotides at sites i and j
    # # this is a matrix
    # # the cov matrix for the two sites is just the mean, across all individuals, of the outer product of P1 and P2
    # #P2 is site 2 with means subtracted out
    for j in range(sites):
        for k in range(sites):
            for i in range(n):
                cov[j][k] += sr.outer_general(P[j][i], P[k][i]) / n
    print("cov")
    naming = os.path.basename(f.name)
    np.save(naming + str('_cov'), cov)

    Pa = array([[[0.0 for z in range(dm)] for j in range(n)] for i in range(sites)])
    reg11 = array([[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(sites)] for k in range(sites)])
    # # Regression of site k on site l
    # # reg11[2][1] is a matrix (hence two more indices) containing the regressions of
    # # each element of the site 2 vector on each element of the site 1 vector
    for k in range(sites):
        for l in range(sites):
            for i in range(0, dm):
                for j in range(0, dm):
                    if var[l][j] > 0.0000000000001:
                        reg11[k][l][i][j] = cov[k][l][i][j] / var[l][j]
                    else:
                        reg11[k][l][i][j] = 0
    r1on2 = reg11[0][1]
    r2on1 = reg11[1][0]

    # # # # # First order terms with zeros except for the value that is
    # # # # # present (i.e. orthogonalized within each vector)
    for k in range(sites):  # site
        for i in range(0, n):  # indiv
            Pa[k][i] = sr.inner_general(sr.outer_general(phi[k][i], phi[k][i]), P[k][i])
            PakiDim = np.array(Pa[k][i])
            # srOuter11Dim = np.array(sr.outer_general(phi[k][i],phi[k][i]), P[k][i])
            # srinnerGenDim = np.array(sr.outer_general(phi[k][i],phi[k][i]), P[k][i])
    # # #print("PakiDim =" + str(PakiDim.ndim))
    #np.save("PakiDim_star_6sites.npy", PakiDim)
    # PakiDim = np.load("PakiDim_star.npy")
    # # # # # Site j orthogonalized wrt site k
    # # # # P1i1[j][k][i] =  first order phi of site j independent of site k for individual i
    P1i1 = array([[[[0.0 for z in range(dm)] for i in range(n)] for j in range(sites)] for k in range(sites)])
    for i in range(0, n):  # Individuals
        for j in range(sites):
            for k in range(sites):
                if k != j:
                    P1i1[j][k][i] = P[j][i] - sr.inner_general(reg11[j][k], Pa[k][i])
    # P1i1 = np.load("P1i1_star_22.npy")
    P2i1 = P1i1[1][0]
    #np.save("P1i1_star_6sites.npy", P1i1)
    # # # # Variance in P2i1
    varP1i1 = array([[[0.0 for z in range(dm)] for i in range(sites)] for j in range(sites)])
    for i in range(0, n):  # individuals
        for j in range(sites):
            for k in range(sites):
                if k != j:
                    varP1i1[j][k] += (P1i1[j][k][i] ** 2) / n

    # # # # cov11i1[j][k][l] = cov between site j and (site k independent of l)
    cov11i1 = array(
        [[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(sites)] for k in range(sites)] for l in
         range(sites)])
    for j in range(sites):
        for k in range(sites):
            for l in range(sites):
                for i in range(n):
                    cov11i1[j][k][l] += sr.outer_general(P[j][i], P1i1[k][l][i]) / n

    # # # # # regression of site j on (site k independent of l)
    reg11i1 = array([[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(sites)] for k in range(sites)] for l in range(sites)])
    for k in range(sites):
        for l in range(sites):
            for m in range(sites):
                for i in range(0, dm):
                    for j in range(0, dm):
                        if varP1i1[l][m][j] > 0.0000000000001:
                            reg11i1[k][l][m][i][j] = cov11i1[k][l][m][i][j] / varP1i1[l][m][j]
                        else:
                            reg11i1[k][l][m][i][j] = 0
    # #dump cov11i1 = none
    cov11i1 = none
    # # # # # same as P1i1, except with all elements = 0 except the one present.
    Pa1i1 = array([[[[0.0 for z in range(dm)] for i in range(n)] for j in range(sites)] for k in range(sites)])
    #
    for i in range(0, n):  # indiv
        for j in range(sites):
            for k in range(sites):
                if k != j:
                    Pa1i1[j][k][i] = sr.inner_general(sr.outer_general(phi[j][i], phi[j][i]), P1i1[j][k][i])
    #np.save("Pa1i1_star_6sites.npy", Pa1i1)
    # Pa1i1 = np.load("Pa1i1_star.npy")
    # # # # # P1D[j][i] = first order poly of site j independent of all other sites, for individual i
    P1D = array([[[0.0 for z in range(dm)] for i in range(n)] for j in range(sites)])
    # #
    for i in range(0, n):  # indiv
        for j in range(sites):
            for k in range(sites):
                if k != j:
                    for l in range(sites):
                        if l != k & l != j:
                            P1D[j][i] = P[j][i] - sr.inner_general(reg11i1[j][k][l], Pa1i1[k][l][i]) - sr.inner_general(reg11[j][l], Pa[l][i])
    #np.save("P1D_star_6sites.npy", P1D)
    # P1D = np.load("P1D_star.npy")
    # reg11i1 = none
    # # # #variance in P1D
    varP1D = array([[0.0 for z in range(dm)] for i in range(sites)])
    #
    for k in range(sites):
        for i in range(0, dm):  # nucleotide
            for j in range(0, n):  # individual
                varP1D[k][i] += ((P1D[k][j][i]) ** 2) / n
    #np.save("varP1D_star_6sites.npy", varP1D)
    # varP1D = np.load("varP1D_star.npy")
    Pa2i1 = Pa1i1[1][0]
    varP2i1 = varP1i1[1][0]
    # # # #-------------------------------------------------------------
    # # # # ------------------------Second Order Terms -----------------
    # # # #-------------------------------------------------------------
    # # # #recursive polynomial builder
    # # #
    # # # def recursive_poly_builder(ord_one, ord_two, d):
    # # #     if ord_one - ord_two == 0:
    # # #         return 0
    # # #     else:
    # # #         return [recursive_poly_builder(ord_one - 1, ord_two, d) for i in range(d)]
    # # #         # uses the concept of list comprehension]
    # # #
    # # # # Second order phenotypes.
    for k in range(n):
        for i in range(sites):
            for j in range(sites):
                if j != i:
                    phi2[i][j][k] = sr.outer(phi[i][k], phi[j][k])
                    phi2m[i][j] += phi2[i][j][k] / n
    phi12 = phi2[0][1]
    phi12m = phi2m[0][1]
    # # # # Q12 contains the 2'nd order phenotypes with the means subtracted out.
    Q2 = array([[[[[0.0 for k in range(dm)] for i in range(dm)] for j in range(n)] for l in range(sites)] for m in range(sites)])
    for i in range(n):  # indiv
        for j in range(sites):
            for k in range(sites):
                if k != j:
                    Q2[j][k][i] = phi2[j][k][i] - phi2m[j][k]
        # Q12[i] = phi12[i] - phi12m
    Q12 = Q2[0][1]

    # # # # Covariance between elements of the 2'nd order phenotype matrix and
    # # # # the 1'st order phenotype.
    cov2w1 = array([[[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)] for k in range(sites)] for l in range(sites)] for m in range(sites)])
    for i in range(n):
        for j in range(sites):
            for k in range(sites):
                if k != j:
                    for l in range(sites):
                        cov2w1[j][k][l] += sr.outer_general(Q2[j][k][i], P[l][i]) / n

    # # # # Covariance of second order phenotype matrices with first order phenotypes.
    cov2w1a = array([[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)] for k in range(sites)] for l in range(sites)])
    cov2w1b = array([[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)] for k in range(sites)] for l in range(sites)])
    #cov2w1c = array([[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)] for k in range(sites)] for l in range(sites)])

    for i in range(n):
        for j in range(sites):
            for k in range(sites):
                if k != j:
                    cov2w1a[j][k] += sr.outer_general(Q2[j][k][i], P[0][i]) / n
                    cov2w1b[j][k] += sr.outer_general(Q2[j][k][i], P1i1[1][0][i]) / n
                    #cov2w1c[j][k] += sr.outer_general(Q2[j][k][i], P1D[2][i]) / n #will need this when doing third order
    # # #print(cov2w1a[0][1]-cov2w1test[0][0][1])

    # # # # regressions of second order phenotype matrices on first order phenotypes.
    r2on1a = array([[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)] for k in range(sites)] for l in range(sites)])
    r2on1b = array([[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)] for k in range(sites)] for l in range(sites)])
    #r2on1c = array([[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)] for k in range(sites)] for l in range(sites)]) #need this when doing third order

    for i in range(sites):
        for j in range(sites):
            if j != i:
                for k in range(dm):
                    for l in range(dm):
                        for m in range(dm):
                            if var[0][m] > 0.0000000001:
                                r2on1a[i][j][k][l][m] = cov2w1a[i][j][k][l][m] / var[0][m]
                            else:
                                r2on1a[i][j][k][l][m] = 0
                            if varP1i1[1][0][m] > 0.0000000001:
                                r2on1b[i][j][k][l][m] = cov2w1b[i][j][k][l][m] / varP1i1[1][0][m]
                            else:
                                r2on1b[i][j][k][l][m] = 0
                            # if varP1D[2][m] > 0.0000000001:
                            #     r2on1c[i][j][k][l][m] = cov2w1c[i][j][k][l][m] / varP1D[2][m]
                            # else:
                            #     r2on1c[i][j][k][l][m] = 0
    # cov2w1a = none
    # cov2w1b = none
    # cov2w1c = none
    # # # Second order polynomials
    P2 = array([[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(n)] for k in range(sites)] for l in range(sites)])
    P1Da = array([[[0.0 for z in range(dm)] for i in range(n)] for j in range(sites)])
    for i in range(sites):
        for j in range(sites):
            if j != i:
                for k in range(n):
                    P2[i][j][k] = Q2[i][j][k] - sr.inner_general(r2on1a[i][j], Pa[0][k]) - sr.inner_general(r2on1b[i][j], Pa1i1[1][0][k]) #- sr.inner_general(r2on1c[i][j], P1Da[2][k])
    r12on1 = r2on1a[0][1]
    r12on2i1 = r2on1b[0][1]
    PP12 = P2[0][1]
    #np.save("P2_star_6sites.npy", P2)
    # # # #print(var2[0][1])
    # # # Second order terms with zeros except for the value that is
    # # # present (i.e. orthogonalized within each matrix)
    P2a = array([[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(n)] for k in range(sites)] for l in
                 range(sites)])
    for h in range(0, n):  # indiv
        for i in range(sites):
            for j in range(sites):
                if j != i:
                    P2a[i][j][h] = sr.inner_general(sr.outer_general(phi2[i][j][h], phi2[i][j][h]), P2[i][j][h])
    #np.save("P2a_star_6sites.npy", P2a)
    PPa12 = P2a[0][1]
    # # # # Covariances between second order phenotypes
    cov2w2 = array([[[[[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)] for k in range(dm)] for l in
                       range(sites)] for m in range(sites)] for n in range(sites)] for p in range(sites)])
    for i in range(n):
        for j in range(sites):
            for k in range(sites):
                if k != j:
                    for l in range(sites):
                        for m in range(sites):
                            if m != l:
                                cov2w2[j][k][l][m] += sr.outer_general(P2[j][k][i], P2[l][m][i]) / n
    #np.save("cov2w2_star_6sites.npy", cov2w2)
    # # Variances of second order phenotypes
    var2 = array([[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(sites)] for k in range(sites)])
    for i in range(sites):
        for j in range(sites):
            if j != i:
                for k in range(dm):
                    for l in range(dm):
                        var2[i][j][k][l] = cov2w2[i][j][i][j][k][l][k][l]
    var12 = var2[0][1]
    #np.save("var2_star_6sites.npy", var2)
    # # # Variances of second order phenotypes
    # # for i in range(sites):
    # #     for j in range(sites):
    # #         if j >> i:
    # #             for k in range(n):
    # #                 for l in range(dm):
    # #                     for m in range(dm):
    # #                         var2[i][j][l][m] += (P2[i][j][k][l][m]**2)/n - (P2m[i][j][l][m]**2)/n
    #
    # # # # regressions of second order phenotypes on one another
    reg2on2 = array([[[[[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)] for k in range(dm)] for l in
                        range(sites)] for m in range(sites)] for n in range(sites)] for p in range(sites)])
    for i in range(sites):
        for j in range(sites):
            if j != i:
                for k in range(sites):
                    for l in range(sites):
                        if l != k:
                            for m in range(dm):
                                for n in range(dm):
                                    for o in range(dm):
                                        for p in range(dm):
                                            if var2[k][l][o][p] > 0.0000000001:
                                                reg2on2[i][j][k][l][m][n][o][p] = cov2w2[i][j][k][l][m][n][o][p] / \
                                                                                  var2[k][l][o][p]
                                            else:
                                                reg2on2[i][j][k][l][m][n][o][p] = 0
    # var2 = none
    # cov2w2 = none
    # # # # Second order phenotypes independent of one another
    P2i2 = array([[[[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(n)] for k in range(sites)] for l in
                    range(sites)] for m in range(sites)] for n in range(sites)])
    P2i2a = array([[[[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(n)] for k in range(sites)] for l in
                     range(sites)] for m in range(sites)] for n in range(sites)])
    for i in range(n):
        for j in range(sites):
            for k in range(sites):
                if k != j:
                    for l in range(sites):
                        for m in range(sites):
                            if m != l:
                                P2i2[j][k][l][m][i] = P2[j][k][i] - sr.inner_general(reg2on2[j][k][l][m], P2a[l][m][i])
                                P2i2a[j][k][l][m][i] = sr.inner_general(sr.outer_general(phi2[j][k][i], phi2[j][k][i]), P2i2[j][k][l][m][i])

    # # # # cov of 2'nd order phi with another independent of the third
    cov2w2i2 = array([[[[[[[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)] for k in range(dm)] for l
                           in range(sites)] for m in range(sites)] for n in range(sites)] for p in range(sites)] for q
                       in range(sites)] for r in range(sites)])
    var2i2 = array([[[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(sites)] for k in range(sites)] for l
                     in range(sites)] for m in range(sites)])
    for i in range(n):
        for j in range(sites):
            for k in range(sites):
                if k != j:
                    for l in range(sites):
                        for m in range(sites):
                            if m != l:
                                for n in range(sites):
                                    for o in range(sites):
                                        if o != n:
                                            cov2w2i2[j][k][l][m][n][o] += sr.outer_general(P2[j][k][i],
                                                                                           P2i2[l][m][n][o][i]) / n
                                for p in range(dm):
                                    for q in range(dm):
                                        var2i2[j][k][l][m][p][q] += (P2i2[j][k][l][m][i][p][q] ** 2) / n
    #np.save("cov2w2i2_star_6sites.npy", cov2w2i2)
    #np.save("var2i2_star_6sites.npy", var2i2)
    # # # # regressions corresponding to the above
    reg2on2i2 = array([[[[[[[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)] for k in range(dm)] for
                            l in range(sites)] for m in range(sites)] for n in range(sites)] for p in range(sites)] for
                        q in range(sites)] for r in range(sites)])
    for i in range(sites):
        for j in range(sites):
            if j != i:
                for k in range(sites):
                    for l in range(sites):
                        if l != k:
                            for k1 in range(sites):
                                for l1 in range(sites):
                                    if l1 != k1:
                                        for m in range(dm):
                                            for n in range(dm):
                                                for o in range(dm):
                                                    for p in range(dm):
                                                        if var2i2[k][l][k1][l1][o][p] > 0.0000000001:
                                                            reg2on2i2[i][j][k][l][k1][l1][m][n][o][p] = \
                                                            cov2w2i2[i][j][k][l][k1][l1][m][n][o][p] / \
                                                            var2i2[k][l][k1][l1][o][p]
                                                        else:
                                                            reg2on2i2[i][j][k][l][k1][l1][m][n][o][p] = 0
    # cov2w2i2 = none
    # var2i2 = none
    #np.save("reg2on2i2_star_6sites.npy", reg2on2i2)
    # # # # 2'nd order phi independent of all others
    P2D = array([[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(n)] for k in range(sites)] for l in
                 range(sites)])
    P2Da = array([[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(n)] for k in range(sites)] for l in
                  range(sites)])
    for i in range(1):
        for j in range(sites):
            for k in range(sites):
                if k != j:
                    for l in range(sites):
                        for m in range(sites):
                            if m != l and (l, m) != (j, k) and (l, m) != (k, j):
                                for n in range(sites):
                                    for o in range(sites):
                                        if o != n and (n, o) != (l, m) and (n, o) != (m, l) and (n, o) != (j, k) and (
                                        n, o) != (k, j):
                                            P2D[j][k][i] = P2[j][k][i] - sr.inner_general(reg2on2i2[j][k][l][m][n][o],P2i2a[l][m][n][o][i]) - sr.inner_general(reg2on2[j][k][n][o], P2a[n][o][i])
                                            P2Da[j][k][i] = sr.inner_general(sr.outer_general(phi2[j][k][i], phi2[j][k][i]), P2D[j][k][i])
    # reg2on2 = none
    # reg2on2i2 = none

    # P2Da = none
    # # # #variance in P2D
    var2D = array([[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(sites)] for k in range(sites)])
    for i in range(n):
        for j in range(sites):
            for k in range(sites):
                if k != j:
                    for l in range(dm):
                        for m in range(dm):
                            var2D[j][k][l][m] += (P2D[j][k][i][l][m] ** 2) / n
    # ------------Projecting the trait values into the space of orthogonal ----
    # -------------polynomials --------------------------------------------------------------------
    # initializing arrays
    cov1FP = array([[0.0 for z in range(dm)] for i in range(sites)])
    covFP = array([[0.0 for z in range(dm)] for i in range(sites)])
    covFPP = array([[0.0 for z in range(dm)] for i in range(dm)])
    # rFon1 = array([0.0 for z in range(dm)])
    # rFon2 = array([0.0 for z in range(dm)])
    rFon2i1 = array([0.0 for z in range(dm)])
    rFon12 = array([[0.0 for z in range(dm)] for i in range(dm)])
    covFw1 = array([[0.0 for z in range(dm)] for i in range(sites)])
    covFw1i1 = array([[[0.0 for z in range(dm)] for i in range(sites)] for j in range(sites)])
    covFw1D = array([[0.0 for z in range(dm)] for i in range(sites)])
    rFon1 = array([[0.0 for z in range(dm)] for i in range(sites)])
    rFon1i1 = array([[[0.0 for z in range(dm)] for i in range(sites)] for j in range(sites)])
    rFon1D = array([[0.0 for z in range(dm)] for i in range(sites)])
    covFw2 = array([[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(sites)] for k in range(sites)])
    covFw2i2 = array([[[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(sites)] for k in range(sites)] for
                       l in range(sites)] for m in range(sites)])
    covFw2D = array([[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(sites)] for k in range(sites)])
    rFon2 = array([[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(sites)] for k in range(sites)])
    rFon2i2 = array([[[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(sites)] for k in range(sites)] for
                      l in range(sites)] for m in range(sites)])
    rFon2D = array([[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(sites)] for k in range(sites)])
    covFw3 = array([[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)])
    rFon3 = array([[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)])
    # Calculating the mean trait value
    Fm = 0
    for i in range(0, n):  # individuals
        Fm += F[i] / n
    # Covariances of the trait with each element of the 1'st order vectors.
    # We can use the 'dot' operator here to get the inner product of a first and a second rank tensor (a vector and a matrix).
    covFP[0] = dot(F, P[0]) / n  # for site 1
    cov1FP[1] = dot(F, P[1]) / n
    covFP[1] = dot(F, P2i1) / n  # for site 2 independent of 1
    # print(covFP[0])
    # print(cov1FP[0])
    # print(covFP[1])
    for i in range(n):
        for j in range(sites):
            for k in range(dm):
                covFw1[j][k] += F[i] * P[j][i][k] / n
                covFw1D[j][k] += F[i] * P1D[j][i][k] / n
            for l in range(sites):
                if l != j:
                    for m in range(dm):
                        covFw1i1[j][l][m] += F[i] * P1i1[j][l][i][m] / n
    for i in range(n):
        for j in range(sites):
            for k in range(sites):
                if k != j:
                    for l in range(dm):
                        for m in range(dm):
                            covFw2[j][k][l][m] += F[i] * P2[j][k][i][l][m] / n
                            covFw2D[j][k][l][m] += F[i] * P2D[j][k][i][l][m] / n
            for n in range(sites):
                for o in range(sites):
                    if n != o:
                        for p in range(dm):
                            for q in range(dm):
                                covFw2i2[j][k][n][o][p][q] += F[i] * P2i2[j][k][n][o][i][p][q] / n
    # for i in range(n):
    #     for j in range(dm):
    #         for k in range(dm):
    #             for l in range(dm):
    #                 covFw3[j][k][l] += F[i] * P3[i][j][k][l] / n

    # Covariance of the trait with each element of the second order phenotype.
    # note: We can nOT use the 'dot' operator here because we do not have a matrix times a vector.
    for i in range(0, n):  # indiv
        covFPP += (F[i] * PP12[i] / n)
    #np.save("covFPP_star_6sites.npy", covFPP)
    # print("covFPP ="+str(covFPP))
    # Regressions of the trait on each element of the first order
    # phenotype vectors.
    for j in range(sites):
        for i in range(0, dm):  # nucleotides
            if var[j][i] > 0.0000000001:
                rFon1[j][i] = covFw1[j][i] / var[j][i]
            else:
                rFon1[j][i] = 0
            if varP1D[j][i] > 0.0000000001:
                rFon1D[j][i] = covFw1D[j][i] / varP1D[j][i]
            else:
                rFon1D[j][i] = 0

    # plt.hist(list(rFon1), color = none, edgecolor = 'black', bins = int(100/5))
    for i in range(sites):
        for j in range(sites):
            if j != i:
                for k in range(dm):
                    if varP1i1[i][j][k] > 0.00000000001:
                        rFon1i1[i][j][k] = covFw1i1[i][j][k] / varP1i1[i][j][k]
                    else:
                        rFon1i1[i][j][k] = 0
    # Regressions of the trait on each element of the second order
    # phenotype matrices.
    for i in range(sites):
        for j in range(sites):
            if j != i:
                for k in range(dm):
                    for l in range(dm):
                        if var2[i][j][k][l] > 0.00000000001:
                            rFon2[i][j][k][l] = covFw2[i][j][k][l] / var2[i][j][k][l]
                        else:
                            rFon2[i][j][k][l] = 0
                        if var2D[i][j][k][l] > 0.00000000001:
                            rFon2D[i][j][k][l] = covFw2D[i][j][k][l] / var2D[i][j][k][l]
                        else:
                            rFon2D[i][j][k][l] = 0
    #we need rFon2D when doing up to 3rd order and just rFon2 when doing up to 2nd order
    np.save(naming + str('_rFon2'), rFon2)  # rFon2D is needed as we're working with 3 sites
    print("computed rFon2")
    np.save(naming + str('_rFon2D'), rFon2D)  # rFon2D is needed as we're working with 3 sites
    print("computed rFon2D")

    for i in range(sites):
        for j in range(sites):
            if j != i:
                for k in range(sites):
                    for l in range(sites):
                        if l != k:
                            for m in range(dm):
                                for n in range(dm):
                                    if var2i2[i][j][k][l][m][n] > 0.0000000001:
                                        rFon2i2[i][j][k][l][m][n] = covFw2i2[i][j][k][l][m][n] / var2i2[i][j][k][l][m][n]
                                    else:
                                        rFon2i2[i][j][k][l][m][n] = 0
    # for i in range(dm):
    #     for j in range(dm):
    #         for k in range(dm):
    #             if var3[i][j][k] > 0.0000000001:
    #                 rFon3[i][j][k] = covFw3[i][j][k] / var3[i][j][k]
    #             else:
    #                 rFon3[i][j][k] = 0

    # print("rFon2i1"+str(rFon2i1))
    #
    # # Regressions of the trait on each element of the second order phenotype matrix.
    for i in range(0, dm):  # nucleotide 1
        for j in range(0, dm):  # nucleotide 2
            if var12[i][j] > 0.0000000001:
                rFon12[i][j] = covFPP[i][j] / var12[i][j]
            else:
                rFon12[i] = 0

    # print("rFon12"+str(rFon12))
    # Contribution of site 1 for each individual.
    # This is the regression of the trait on site 1 times (inner product)
    # the individual's site 1 vector that has been orthogonalized within
    # the vector.
    for i in range(0, n):
        # test_val = sr.inner_general(rFon1[0],Pa[0][i])
        Fon1[i] = sr.inner_general(rFon1[0], Pa[0][i])

    # Contribution of site 2 independent of 1 for each individual.
    for i in range(0, n):
        Fon2i1[i] = sr.inner_general(rFon1i1[1][0], Pa1i1[1][0][i])

    # # Contribution of the second order phenotype for each individual.
    for i in range(0, n):
        Fon12[i] = sr.inner_general(rFon2[0][1], P2a[0][1][i])

    # contribution of third order phenotype for each individual......
    # for i in range(0, n):
    #     Fon3[i] = sr.inner_general(rFon3[0], P3a)
    # Ignoring very small values that would be due to roundoff error.
    # Change or delete this for a large data set.
    for i in range(0, n):
        if fabs(Fon1[i]) < 0.0000000000001:
            Fon1[i] = 0
        if fabs(Fon2i1[i]) < 0.0000000000001:
            Fon2i1[i] = 0
        if fabs(Fon12[i]) < 0.0000000000001:
            Fon12[i] = 0
    # -----------------------Listing the main results------------------
    print('Regression of trait on site 1')
    print(rFon1)
    print('Regression on 1st order polynomial - orthogonalized within - rFon1D')
    print(rFon1D)
    print('Regression of trait on site 2')
    print(rFon2)
    print('Regression on 2nd order polynomial - orthogonalized within - rFon2D')
    print(rFon2D)
    print('Regression of trait on site 2 independent of 1')
    print(rFon2i1)
    print('Regression on (site 1)x(site 2), independent of first order')
    print(rFon12)
    #print('Regression of trait on site 3')
    #print(rFon3)
    print("--- %s seconds ---" % (time.time() - start_time))

if __name__ == '__main__':
    orthogonal_polynomial()
