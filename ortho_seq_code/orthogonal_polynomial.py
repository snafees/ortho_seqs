import numpy as np
from numpy.linalg import *
import time
import os
import ortho_seq_code.sr as sr
from ortho_seq_code.constants_orthoseqs import *
import click
import itertools


def create_dir_if_not_exists(out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    else:
        print("Path {} already exists, might be overwriting data".format(out_dir))


def orthogonal_polynomial(
    filename,
    pheno_file,
    molecule,
    sites,
    dm,
    pop_size,
    poly_order,
    precomputed,
    out_dir,
):
    """Program to compute orthogonal polynomials up to 2nd order"""
    create_dir_if_not_exists(out_dir)
    start_time = time.time()
    with open(filename) as f:
        seq = f.readlines()
    global i
    # file containing trait values that will be mapped to sequence
    # vectors that must be the same size as F
    with open(pheno_file) as f2:
        phenotype = f2.readlines()
    print(phenotype)
    F = np.genfromtxt(phenotype)  # this needs to stay this way!
    Fest = np.genfromtxt(phenotype)  # this needs to stay this way!
    Fon1 = np.genfromtxt(phenotype)  # this needs to stay this way!
    Fon2i1 = np.genfromtxt(phenotype)  # this needs to stay this way!
    Fon12 = np.genfromtxt(phenotype)  # this needs to stay this way!
    print(Fest)
    Fest = [0] * pop_size
    Fon1 = [0] * pop_size
    Fon2i1 = [0] * pop_size
    Fon12 = [0] * pop_size
    np.set_printoptions(precision=10)
    range_sites = range(sites)
    range_popsize = range(pop_size)
    range_dm = range(dm)
    # ----Initializing various terms that we will use.--------------
    # 3 sites, each a dm dim vector, in n individuals
    # nOTE: For application to Amino Acid sequences, increase
    # the size of the arrays accordingly.
    phi = np.zeros((sites, pop_size, dm))
    mean = np.zeros((sites, dm))
    var = np.zeros((sites, dm))
    phi2 = np.zeros((sites, sites, pop_size, dm, dm))
    phi2m = np.zeros((sites, sites, dm, dm))
    P = np.zeros((sites, pop_size, dm))
    cov = np.zeros((sites, sites, dm, dm))
    # ------------Converting letters to vectors---------------
    # phi[individual][site][state]. phi[i][j] = vector for site j in individual i.
    alphabets = DM_ALPHABETS[dm]
    for dna_alphabet_index in range(len(alphabets)):
        for i in range_popsize:
            for j in range_sites:
                if seq[i][j] == alphabets[dna_alphabet_index]:
                    phi[j][i][dna_alphabet_index] = 1.0

    if dm == 3 and molecule == "protein_pnp":
        iterator = itertools.product(
            range(len(PROTEIN_ALPHABETS_POLAR)),
            range(len(PROTEIN_ALPHABETS_NONPOLAR)),
            range_popsize,
            range_sites,
        )
        for (
            protein_alphabet_index_polar,
            protein_alphabet_index_nonpolar,
            i,
            j,
        ) in iterator:
            if seq[i][j] == PROTEIN_ALPHABETS_POLAR[protein_alphabet_index_polar]:
                phi[j][i][0] = 1.0
            if seq[i][j] == PROTEIN_ALPHABETS_NONPOLAR[protein_alphabet_index_nonpolar]:
                phi[j][i][1] = 1.0
            if seq[i][j] == "n":
                phi[j][i][2] = 1.0

    naming = os.path.basename(f.name)
    if precomputed:
        precomputed_array = np.load(os.path.join(out_dir, naming + ".npz"))
        mean = precomputed_array[naming + "_mean"]
        P = precomputed_array[naming + "_P"]
        var = precomputed_array[naming + "_var"]
        cov = precomputed_array[naming + "_cov"]
        reg11 = precomputed_array[naming + "_reg11"]
        Pa = precomputed_array[naming + "_Pa"]
        P1i1 = precomputed_array[naming + "_P1i1"]
        P2i1 = P1i1[1][0]
        varP1i1 = precomputed_array[naming + "_varP1i1"]
        cov11i1 = precomputed_array[naming + "_cov11i1"]
        reg11i1 = precomputed_array[naming + "_reg11i1"]
        Pa1i1 = precomputed_array[naming + "_Pa1i1"]
        P1D = precomputed_array[naming + "_P1D"]
        varP1D = precomputed_array[naming + "_varP1D"]
        phi2 = precomputed_array[naming + "_phi2"]
        phi2m = precomputed_array[naming + "_phi2m"]
        Q2 = precomputed_array[naming + "_Q2"]
        cov2w1 = precomputed_array[naming + "_cov2w1"]
        cov2w1a = precomputed_array[naming + "_cov2w1a"]
        cov2w1b = precomputed_array[naming + "_cov2w1b"]
        r2on1a = precomputed_array[naming + "_r2on1a"]
        r2on1b = precomputed_array[naming + "_r2on1b"]
        P2 = precomputed_array[naming + "_P2"]
        PP12 = P2[0][1]
        P2a = precomputed_array[naming + "_P2a"]
        cov2w2 = precomputed_array[naming + "_cov2w2"]
        var2 = precomputed_array[naming + "_var2"]
        var12 = var2[0][1]
        reg2on2 = precomputed_array[naming + "_reg2on2"]
        P2i2 = precomputed_array[naming + "_P2i2"]
        P2i2a = precomputed_array[naming + "_P2i2a"]
        cov2w2i2 = precomputed_array[naming + "_cov2w2i2"]
        var2i2 = precomputed_array[naming + "_var2i2"]
        reg2on2i2 = precomputed_array[naming + "_reg2on2i2"]
        P2D = precomputed_array[naming + "_P2D"]
        P2Da = precomputed_array[naming + "_P2Da"]
        var2D = precomputed_array[naming + "_var2D"]
    else:
        # keep in alpha order
        # ---------------------------------First order terms ----------------------
        # calculate mean vectors first
        arrays_save = {}
        if poly_order == "first" and not precomputed:
            for i, j in itertools.product(range_popsize, range_sites):
                mean[j] += phi[j][i] / pop_size
            arrays_save[naming + "_mean"] = mean
            #  to show progress, can do something much more efficient/elegant
            print("computed mean")
            for j, i in itertools.product(range_sites, range_popsize):  # site, indiv
                P[j][i] = phi[j][i] - mean[j]
            arrays_save[naming + "_P"] = P

            # var[site][nucleotide]
            for k, i, j in itertools.product(
                range_sites, range_dm, range_popsize
            ):  # nucleotide, indiv
                var[k][i] += ((P[k][j][i]) ** 2) / pop_size
            print("computed variance")
            arrays_save[naming + "_var"] = var

            # Covariances between nucleotides at sites i and j
            # this is a matrix
            # the cov matrix for the two sites is just the mean,
            # across all individuals, of the outer product of P1 and P2
            # #P2 is site 2 with means subtracted out
            for j, k, i in itertools.product(range_sites, range_sites, range_popsize):
                cov[j][k] += sr.outer_general(P[j][i], P[k][i]) / pop_size
            print("computed covariance")
            arrays_save[naming + "_cov"] = cov

            Pa = np.zeros((sites, pop_size, dm))
            reg11 = np.zeros((sites, sites, dm, dm))
            # Regression of site k on site l
            # reg11[2][1] is a matrix (hence two more indices)
            # containing the regressions of
            # each element of the site 2 vector on each element of the site 1 vector
            for k, l, i, j in itertools.product(
                range_sites, range_sites, range_dm, range_dm
            ):
                if var[l][j] > 0.0000000000001:
                    reg11[k][l][i][j] = cov[k][l][i][j] / var[l][j]
                else:
                    reg11[k][l][i][j] = 0
            print("computed reg11")
            arrays_save[naming + "_reg11"] = reg11

            # # # # # First order terms with zeros except for the value that is
            # # # # # present (i.e. orthogonalized within each vector)
            for k, i in itertools.product(range_sites, range_popsize):  # site, indiv
                Pa[k][i] = sr.inner_general(
                    sr.outer_general(phi[k][i], phi[k][i]), P[k][i]
                )
            print("computed Pa: first order orthogonalized within each vector")
            arrays_save[naming + "_Pa"] = Pa

            # Site j orthogonalized wrt site k
            # P1i1[j][k][i] =  first order phi of site j independent of
            # site k for individual i
            P1i1 = np.zeros((sites, sites, pop_size, dm))
            for i, j, k in itertools.product(range_popsize, range_sites, range_sites):
                if k != j:
                    P1i1[j][k][i] = P[j][i] - sr.inner_general(reg11[j][k], Pa[k][i])
            print("computed P1i1")
            arrays_save[naming + "_P1i1"] = P1i1
            P2i1 = P1i1[1][0]

            # # # # Variance in P2i1
            varP1i1 = np.zeros((sites, sites, dm))
            for i, j, k in itertools.product(range_popsize, range_sites, range_sites):
                if k != j:
                    varP1i1[j][k] += (P1i1[j][k][i] ** 2) / pop_size
            print("computed varP1i1")
            arrays_save[naming + "_varP1i1"] = varP1i1

            # # # # cov11i1[j][k][l] = cov between site j and (site k independent of l)
            cov11i1 = np.zeros((sites, sites, sites, dm, dm))
            for j, k, l, i in itertools.product(
                range_sites, range_sites, range_sites, range_popsize
            ):
                cov11i1[j][k][l] += sr.outer_general(P[j][i], P1i1[k][l][i]) / pop_size
            print("computed cov11i1")
            arrays_save[naming + "_cov11i1"] = cov11i1
            # # # # # regression of site j on (site k independent of l)
            reg11i1 = np.zeros((sites, sites, sites, dm, dm))
            for k, l, m, i, j in itertools.product(
                range_sites, range_sites, range_sites, range_dm, range_dm
            ):
                if varP1i1[l][m][j] > 0.0000000000001:
                    reg11i1[k][l][m][i][j] = cov11i1[k][l][m][i][j] / varP1i1[l][m][j]
                else:
                    reg11i1[k][l][m][i][j] = 0
            print("computed reg11i1")
            arrays_save[naming + "_reg11i1"] = reg11i1

            # # # # # same as P1i1, except with all elements = 0 except the one present
            Pa1i1 = np.zeros((sites, sites, pop_size, dm))
            for i, j, k in itertools.product(range_popsize, range_sites, range_sites):
                if k != j:
                    Pa1i1[j][k][i] = sr.inner_general(
                        sr.outer_general(phi[j][i], phi[j][i]), P1i1[j][k][i]
                    )
            print("computed Pa1i1")
            arrays_save[naming + "_Pa1i1"] = Pa1i1
            # # # # # P1D[j][i] = first order poly of site j independent of all other
            # sites, for individual i
            P1D = np.zeros((sites, pop_size, dm))
            for i, j, k in itertools.product(range_popsize, range_sites, range_sites):
                if k != j:
                    for l in range_sites:
                        if l != k & l != j:
                            P1D[j][i] = (
                                P[j][i]
                                - sr.inner_general(reg11i1[j][k][l], Pa1i1[k][l][i])
                                - sr.inner_general(reg11[j][l], Pa[l][i])
                            )
            print("computed P1D")
            arrays_save[naming + "_P1D"] = P1D

            # # # #variance in P1D
            varP1D = np.zeros((sites, dm))
            for k, i, j in itertools.product(range_sites, range_dm, range_popsize):
                varP1D[k][i] += ((P1D[k][j][i]) ** 2) / pop_size
            print("computed varP1D")
            arrays_save[naming + "_varP1D"] = varP1D
        # Pa2i1 = Pa1i1[1][0]
        # varP2i1 = varP1i1[1][0]
        # # # #-------------------------------------------------------------
        # # # # ------------------------Second Order Terms -----------------
        # # # #-------------------------------------------------------------

        # # # # Second order phenotypes.
        if poly_order == "second" and not precomputed:
            # this is also written under first order
            # calculate mean vectors
            for i, j in itertools.product(range_popsize, range_sites):
                mean[j] += phi[j][i] / pop_size
            arrays_save[naming + "_mean"] = mean
            #  to show progress, can do something much more efficient/elegant
            print("computed mean")

            for j, i in itertools.product(range_sites, range_popsize):
                P[j][i] = phi[j][i] - mean[j]
            arrays_save[naming + "_P"] = P

            # var[site][nucleotide]
            for k, i, j in itertools.product(range_sites, range_dm, range_popsize):
                var[k][i] += ((P[k][j][i]) ** 2) / pop_size
            print("computed variance")
            arrays_save[naming + "_var"] = var

            # Covariances between nucleotides at sites i and j
            # this is a matrix
            # the cov matrix for the two sites is just the mean,
            # across all individuals, of the outer product of P1 and P2
            # #P2 is site 2 with means subtracted out
            for j, k, i in itertools.product(range_sites, range_sites, range_popsize):
                cov[j][k] += sr.outer_general(P[j][i], P[k][i]) / pop_size
            print("computed covariance")
            arrays_save[naming + "_cov"] = cov

            Pa = np.zeros((sites, pop_size, dm))
            reg11 = np.zeros((sites, sites, dm, dm))
            # Regression of site k on site l
            # reg11[2][1] is a matrix (hence two more indices)
            # containing the regressions of
            # each element of the site 2 vector on each element of the site 1 vector
            for k, l, i, j in itertools.product(
                range_sites, range_sites, range_dm, range_dm
            ):
                if var[l][j] > 0.0000000000001:
                    reg11[k][l][i][j] = cov[k][l][i][j] / var[l][j]
                else:
                    reg11[k][l][i][j] = 0
            print("computed reg11")
            arrays_save[naming + "_reg11"] = reg11

            # # # # # First order terms with zeros except for the value that is
            # # # # # present (i.e. orthogonalized within each vector)
            for k, i in itertools.product(range_sites, range_popsize):  # indiv
                Pa[k][i] = sr.inner_general(
                    sr.outer_general(phi[k][i], phi[k][i]), P[k][i]
                )
            print("computed Pa: first order orthogonalized within each vector")
            arrays_save[naming + "_Pa"] = Pa

            P1i1 = np.zeros((sites, sites, pop_size, dm))
            for i, j, k in itertools.product(range_popsize, range_sites, range_sites):
                if k != j:
                    P1i1[j][k][i] = P[j][i] - sr.inner_general(reg11[j][k], Pa[k][i])
            P2i1 = P1i1[1][0]
            print("computed P1i1")
            arrays_save[naming + "_P1i1"] = P1i1

            # # # # Variance in P2i1
            varP1i1 = np.zeros((sites, sites, dm))
            for i, j, k in itertools.product(range_popsize, range_sites, range_sites):
                if k != j:
                    varP1i1[j][k] += (P1i1[j][k][i] ** 2) / pop_size
            print("computed varP1i1")
            arrays_save[naming + "_varP1i1"] = varP1i1
            # # # # cov11i1[j][k][l] = cov between site j and (site k independent of l)
            cov11i1 = np.zeros((sites, sites, sites, dm, dm))
            for j, k, l, i in itertools.product(
                range_sites, range_sites, range_sites, range_popsize
            ):
                cov11i1[j][k][l] += sr.outer_general(P[j][i], P1i1[k][l][i]) / pop_size
            print("computed cov11i1")
            arrays_save[naming + "_cov11i1"] = cov11i1
            # # # # # regression of site j on (site k independent of l)
            reg11i1 = np.zeros((sites, sites, sites, dm, dm))
            for k, l, m, i, j in itertools.product(
                range_sites, range_sites, range_sites, range_dm, range_dm
            ):
                if varP1i1[l][m][j] > 0.0000000000001:
                    reg11i1[k][l][m][i][j] = cov11i1[k][l][m][i][j] / varP1i1[l][m][j]
                else:
                    reg11i1[k][l][m][i][j] = 0
            print("computed reg11i1")
            arrays_save[naming + "_reg11i1"] = reg11i1
            # # # # # same as P1i1, except with all elements = 0 except the one present
            Pa1i1 = np.zeros((sites, sites, pop_size, dm))
            for i, j, k in itertools.product(range_popsize, range_sites, range_sites):
                if k != j:
                    Pa1i1[j][k][i] = sr.inner_general(
                        sr.outer_general(phi[j][i], phi[j][i]), P1i1[j][k][i]
                    )
            print("computed Pa1i1")
            arrays_save[naming + "_Pa1i1"] = Pa1i1
            # # # # # P1D[j][i] = first order poly of site j independent of all other
            # sites, for individual i
            P1D = np.zeros((sites, pop_size, dm))
            for i, j, k in itertools.product(
                range_popsize, range_sites, range_sites
            ):  # indiv
                if k != j:
                    for l in range_sites:
                        if l != k & l != j:
                            P1D[j][i] = (
                                P[j][i]
                                - sr.inner_general(reg11i1[j][k][l], Pa1i1[k][l][i])
                                - sr.inner_general(reg11[j][l], Pa[l][i])
                            )
            print("computed P1D")
            arrays_save[naming + "_P1D"] = P1D
            # # # #variance in P1D
            varP1D = np.zeros((sites, dm))
            for k, i, j in itertools.product(range_sites, range_dm, range_popsize):
                varP1D[k][i] += ((P1D[k][j][i]) ** 2) / pop_size
            print("computed varP1D")
            arrays_save[naming + "_varP1D"] = varP1D
            #  end of that part from first order

            for k, i, j in itertools.product(range_popsize, range_sites, range_sites):
                if j != i:
                    phi2[i][j][k] = sr.outer_general(phi[i][k], phi[j][k])
                    phi2m[i][j] += phi2[i][j][k] / pop_size
            # phi12 = phi2[0][1]
            # phi12m = phi2m[0][1]
            print("computed phi2 and phi2m")
            arrays_save[naming + "_phi2"] = phi2
            arrays_save[naming + "_phi2m"] = phi2m

            # Q12 contains the 2'nd order phenotypes with the means subtracted out
            Q2 = np.zeros((sites, sites, pop_size, dm, dm))
            for i, j, k in itertools.product(range_popsize, range_sites, range_sites):
                if k != j:
                    Q2[j][k][i] = phi2[j][k][i] - phi2m[j][k]
            # Q12[i] = phi12[i] - phi12m
            # Q12 = Q2[0][1]
            print("computed Q2")
            arrays_save[naming + "_Q2"] = Q2
            # # # # Covariance between elements of the 2'nd order phenotype matrix and
            # # # # the 1'st order phenotype.
            cov2w1 = np.zeros((sites, sites, sites, dm, dm, dm))
            for i, j, k in itertools.product(range_popsize, range_sites, range_sites):
                if k != j:
                    for l in range_sites:
                        cov2w1[j][k][l] += (
                            sr.outer_general(Q2[j][k][i], P[l][i]) / pop_size
                        )
            print("computed cov2w1")
            arrays_save[naming + "_cov2w1"] = cov2w1
            # # # # Covariance of second order phenotype matrices with first
            # order phenotypes.
            cov2w1a = np.zeros((sites, sites, dm, dm, dm))
            cov2w1b = np.zeros((sites, sites, dm, dm, dm))
            # will need this when doing third order
            # cov2w1c = array([[[[[0.0 for z in range_dm] for i in range_dm]
            # for j in range_dm] for k in range_sites] for l in range_sites])

            for i, j, k in itertools.product(range_popsize, range_sites, range_sites):
                if k != j:
                    cov2w1a[j][k] += sr.outer_general(Q2[j][k][i], P[0][i]) / pop_size
                    cov2w1b[j][k] += (
                        sr.outer_general(Q2[j][k][i], P1i1[1][0][i]) / pop_size
                    )
                    # will need this when doing third order
                    # cov2w1c[j][k] += sr.outer_general(
                    # Q2[j][k][i], P1D[2][i]) / n
            print("computed cov2w1a,cov2w1b")
            arrays_save[naming + "_cov2w1a"] = cov2w1a
            arrays_save[naming + "_cov2w1b"] = cov2w1b

            # regressions of second order phenotype matrices on first order phenotypes.
            r2on1a = np.zeros((sites, sites, dm, dm, dm))
            r2on1b = np.zeros((sites, sites, dm, dm, dm))
            # will need this when doing third order
            # r2on1c = array([[[[[0.0 for z in range_dm] for i in range_dm] for j
            # in range_dm] for k in range_sites] for l in range_sites])

            for i, j in itertools.product(range_sites, range_sites):
                if j != i:
                    for k, l, m in itertools.product(range_dm, range_dm, range_dm):
                        if var[0][m] > 0.0000000001:
                            r2on1a[i][j][k][l][m] = cov2w1a[i][j][k][l][m] / var[0][m]
                        else:
                            r2on1a[i][j][k][l][m] = 0
                        if varP1i1[1][0][m] > 0.0000000001:
                            r2on1b[i][j][k][l][m] = (
                                cov2w1b[i][j][k][l][m] / varP1i1[1][0][m]
                            )
                        else:
                            r2on1b[i][j][k][l][m] = 0
                        # Need for third degree polynomial
                        # if varP1D[2][m] > 0.0000000001:
                        #     r2on1c[i][j][k][l][m] = \
                        # cov2w1c[i][j][k][l][m] / varP1D[2][m]
                        # else:
                        #     r2on1c[i][j][k][l][m] = 0
            print("computed r2on1a, r2on1b")
            arrays_save[naming + "_r2on1a"] = r2on1a
            arrays_save[naming + "_r2on1b"] = r2on1b

            # # # Second order polynomials
            P2 = np.zeros((sites, sites, pop_size, dm, dm))
            # P1Da = np.array(
            #   [[[0.0 for z in range_dm] for i in range_popsize]
            # for j in range_sites])
            for i, j in itertools.product(range_sites, range_sites):
                if j != i:
                    for k in range_popsize:
                        P2[i][j][k] = (
                            Q2[i][j][k]
                            - sr.inner_general(r2on1a[i][j], Pa[0][k])
                            - sr.inner_general(r2on1b[i][j], Pa1i1[1][0][k])
                        )
                        # - sr.inner_general(r2on1c[i][j], P1Da[2][k]) # noqa
            # r12on1 = r2on1a[0][1]
            # r12on2i1 = r2on1b[0][1]
            print("computed P2")
            PP12 = P2[0][1]
            arrays_save[naming + "_P2"] = P2
            # # # Second order terms with zeros except for the value that is
            # # # present (i.e. orthogonalized within each matrix)
            P2a = np.zeros((sites, sites, pop_size, dm, dm))
            for h, i, j in itertools.product(range_popsize, range_sites, range_sites):
                if j != i:
                    P2a[i][j][h] = sr.inner_general(
                        sr.outer_general(phi2[i][j][h], phi2[i][j][h]), P2[i][j][h]
                    )
            print("computed P2a")
            arrays_save[naming + "_P2a"] = P2a
            # PPa12 = P2a[0][1]
            # # # # Covariances between second order phenotypes
            cov2w2 = np.zeros((sites, sites, sites, sites, dm, dm, dm, dm))
            for i, j, k in itertools.product(range_popsize, range_sites, range_sites):
                if k != j:
                    for l, m in itertools.product(range_sites, range_sites):
                        if m != l:
                            cov2w2[j][k][l][m] += (
                                sr.outer_general(P2[j][k][i], P2[l][m][i]) / pop_size
                            )
            print("computed cov2w2")
            arrays_save[naming + "_cov2w2"] = cov2w2
            # # Variances of second order phenotypes
            var2 = np.zeros((sites, sites, dm, dm))
            for i, j in itertools.product(range_sites, range_sites):
                if j != i:
                    for k, l in itertools.product(range_dm, range_dm):
                        var2[i][j][k][l] = cov2w2[i][j][i][j][k][l][k][l]
            var12 = var2[0][1]
            print("computed var2")
            arrays_save[naming + "_var2"] = var2
            # #Variances of second order phenotypes (this is another way of computing variances, by using the polynomial itself)
            # # for i in range_sites:
            # #     for j in range_sites:
            # #         if j >> i:
            # #             for k in range_popsize:
            # #                 for l in range_dm:
            # #                     for m in range_dm:
            # #                         var2[i][j][l][m] += \
            # #   (P2[i][j][k][l][m]**2)/n - (P2m[i][j][l][m]**2)/n

            # # # # regressions of second order phenotypes on one another
            reg2on2 = np.zeros((sites, sites, sites, sites, dm, dm, dm, dm))
            for i, j in itertools.product(range_sites, range_sites):
                if j != i:
                    for k, l in itertools.product(range_sites, range_sites):
                        if l != k:
                            for m, n, o, p in itertools.product(
                                range_dm, range_dm, range_dm, range_dm
                            ):
                                if var2[k][l][o][p] > 0.0000000001:
                                    numerator = cov2w2[i][j][k][l][m][n][o][p]
                                    denominator = var2[k][l][o][p]
                                    reg2on2[i][j][k][l][m][n][o][p] = (
                                        numerator / denominator
                                    )
                                else:
                                    reg2on2[i][j][k][l][m][n][o][p] = 0
            print("computed reg2on2")
            arrays_save[naming + "_reg2on2"] = reg2on2
            # # # # Second order phenotypes independent of one another
            P2i2 = np.zeros((sites, sites, sites, sites, pop_size, dm, dm))
            P2i2a = np.zeros((sites, sites, sites, sites, pop_size, dm, dm))
            for i, j, k in itertools.product(range_popsize, range_sites, range_sites):
                if k != j:
                    for l, m in itertools.product(range_sites, range_sites):
                        if m != l:
                            P2i2[j][k][l][m][i] = P2[j][k][i] - sr.inner_general(
                                reg2on2[j][k][l][m], P2a[l][m][i]
                            )
                            P2i2a[j][k][l][m][i] = sr.inner_general(
                                sr.outer_general(phi2[j][k][i], phi2[j][k][i]),
                                P2i2[j][k][l][m][i],
                            )
            print("computed P2i2, P2i2a")
            arrays_save[naming + "_P2i2"] = P2i2
            arrays_save[naming + "_P2i2a"] = P2i2a
            # # # # cov of 2'nd order phi with another independent of the third
            cov2w2i2 = np.zeros(
                (sites, sites, sites, sites, sites, sites, dm, dm, dm, dm)
            )
            var2i2 = np.zeros((sites, sites, sites, sites, dm, dm))
            for i in range(pop_size):
                for j in range(sites):
                    for k in range(sites):
                        if k != j:
                            for l in range(sites):
                                for m in range(sites):
                                    if m != l:
                                        for n in range(sites):
                                            for o in range(sites):
                                                if o != n:
                                                    cov2w2i2[j][k][l][m][n][o] += (
                                                        sr.outer_general(
                                                            P2[j][k][i],
                                                            P2i2[l][m][n][o][i],
                                                        )
                                                        / pop_size
                                                    )
                                        for p in range(dm):
                                            for q in range(dm):
                                                var2i2[j][k][l][m][p][q] += (
                                                    P2i2[j][k][l][m][i][p][q] ** 2
                                                ) / pop_size
            print("computed cov2w2i2, var2i2")
            arrays_save[naming + "_cov2w2i2"] = cov2w2i2
            arrays_save[naming + "_var2i2"] = var2i2
            # # # # regressions corresponding to the above
            reg2on2i2 = np.zeros(
                (sites, sites, sites, sites, sites, sites, dm, dm, dm, dm)
            )
            for i, j in itertools.product(range_sites, range_sites):
                if j != i:
                    for k, l in itertools.product(range_sites, range_sites):
                        if l != k:
                            for k1, l1 in itertools.product(range_sites, range_sites):
                                if l1 != k1:
                                    for m, n, o, p in itertools.product(
                                        range_dm, range_dm, range_dm, range_dm
                                    ):
                                        if var2i2[k][l][k1][l1][o][p] > 0.0000000001:
                                            numerator = cov2w2i2[i][j][k][l][k1][l1][m][
                                                n
                                            ][o][p]
                                            denominator = var2i2[k][l][k1][l1][o][p]
                                            reg2on2i2[i][j][k][l][k1][l1][m][n][o][
                                                p
                                            ] = (numerator / denominator)
                                        else:
                                            reg2on2i2[i][j][k][l][k1][l1][m][n][o][
                                                p
                                            ] = 0
            print("computed reg2on2i2")
            arrays_save[naming + "_reg2on2i2"] = reg2on2i2
            # # # # 2'nd order phi independent of all others
            P2D = np.zeros((sites, sites, pop_size, dm, dm))
            P2Da = np.zeros((sites, sites, pop_size, dm, dm))
            for i, j, k in itertools.product(range(1), range_sites, range_sites):
                if k != j:
                    for l, m in itertools.product(range_sites, range_sites):
                        if m != l and (l, m) != (j, k) and (l, m) != (k, j):
                            for n, o in itertools.product(range_sites, range_sites):
                                if (
                                    o != n
                                    and (n, o) != (l, m)
                                    and (n, o) != (m, l)
                                    and (n, o) != (j, k)
                                    and (n, o) != (k, j)
                                ):
                                    P2D[j][k][i] = (
                                        P2[j][k][i]
                                        - sr.inner_general(
                                            reg2on2i2[j][k][l][m][n][o],
                                            P2i2a[l][m][n][o][i],
                                        )
                                        - sr.inner_general(
                                            reg2on2[j][k][n][o],
                                            P2a[n][o][i],
                                        )
                                    )
                                    P2Da[j][k][i] = sr.inner_general(
                                        sr.outer_general(phi2[j][k][i], phi2[j][k][i]),
                                        P2D[j][k][i],
                                    )
            print("computed P2D, P2Da")
            arrays_save[naming + "_P2D"] = P2D
            arrays_save[naming + "_P2Da"] = P2Da
            # # # #variance in P2D
            var2D = np.zeros((sites, sites, dm, dm))
            for i, j, k in itertools.product(range_popsize, range_sites, range_sites):
                if k != j:
                    for l, m in itertools.product(range_dm, range_dm):
                        var2D[j][k][l][m] += (P2D[j][k][i][l][m] ** 2) / pop_size
            print("computed var2D")
            arrays_save[naming + "_var2D"] = var2D
        output_npz_file = os.path.join(out_dir, naming + ".npz")
        print("Saving to {}".format(output_npz_file))

        np.savez_compressed(output_npz_file, **arrays_save)
    # ------------Projecting the trait values into the space of orthogonal ----
    # -------------polynomials ------------------------------------------------
    # initializing arrays
    cov1FP = np.zeros((sites, dm))
    covFP = np.zeros((sites, dm))
    covFPP = np.zeros((dm, dm))
    rFon2i1 = np.zeros((dm))
    rFon12 = np.zeros((dm, dm))
    covFw1 = np.zeros((sites, dm))
    covFw1i1 = np.zeros((sites, sites, dm))
    covFw1D = np.zeros((sites, dm))
    rFon1 = np.zeros((sites, dm))
    rFon1i1 = np.zeros((sites, sites, dm))
    rFon1D = np.zeros((sites, dm))
    covFw2 = np.zeros((sites, sites, dm, dm))
    covFw2i2 = np.zeros((sites, sites, sites, sites, dm, dm))
    covFw2D = np.zeros((sites, sites, dm, dm))
    rFon2 = np.zeros((sites, sites, dm, dm))
    rFon2i2 = np.zeros((sites, sites, sites, sites, dm, dm))
    rFon2D = np.zeros((sites, sites, dm, dm))

    # Need for third degree polynomial
    # covFw3 = np.array(
    #     [[[0.0 for z in range_dm] for i in range_dm] for j in range_dm])
    # rFon3 = np.array(
    #     [[[0.0 for z in range_dm] for i in range_dm] for j in range_dm])
    # Calculating the mean trait value

    naming_phenotype = os.path.basename(f2.name)
    cov_with_F_save = {}
    # calculating mean phenotype value
    Fm = sum([F[i] / pop_size for i in range_popsize])
    np.save(os.path.join(out_dir, naming_phenotype + str("_Fm")), Fm)
    # Covariances of the trait with each element of the 1'st order vectors.
    # We can use the 'dot' operator here to get the inner product of a
    # first and a second rank tensor (a vector and a matrix).
    if poly_order == "first":
        covFP[0] = np.dot(F, P[0]) / pop_size  # for site 1
        cov1FP[1] = np.dot(F, P[1]) / pop_size
        covFP[1] = np.dot(F, P2i1) / pop_size  # for site 2 independent of 1
        cov_with_F_save[naming_phenotype + "_covFP[0]"] = covFP[0]
        cov_with_F_save[naming_phenotype + "_cov1FP[1]"] = cov1FP[1]
        cov_with_F_save[naming_phenotype + "_covFP[1]"] = covFP[1]

        for i in range(pop_size):
            for j in range(sites):
                for k in range(dm):
                    covFw1[j][k] += F[i] * P[j][i][k] / pop_size
                    covFw1D[j][k] += F[i] * P1D[j][i][k] / pop_size
                for l in range(sites):
                    if l != j:
                        for m in range(dm):
                            covFw1i1[j][l][m] += F[i] * P1i1[j][l][i][m] / pop_size
        cov_with_F_save[naming_phenotype + "_covFw1i1"] = covFw1i1
    if poly_order == "second":
        # this part is from first order
        covFP[0] = np.dot(F, P[0]) / pop_size  # for site 1
        cov1FP[1] = np.dot(F, P[1]) / pop_size
        covFP[1] = np.dot(F, P2i1) / pop_size  # for site 2 independent of 1
        cov_with_F_save[naming_phenotype + "_covFP[0]"] = covFP[0]
        cov_with_F_save[naming_phenotype + "_cov1FP[1]"] = cov1FP[1]
        cov_with_F_save[naming_phenotype + "_covFP[1]"] = covFP[1]

        for i in range(pop_size):
            for j in range(sites):
                for k in range(dm):
                    covFw1[j][k] += F[i] * P[j][i][k] / pop_size
                    covFw1D[j][k] += F[i] * P1D[j][i][k] / pop_size
                for l in range(sites):
                    if l != j:
                        for m in range(dm):
                            covFw1i1[j][l][m] += F[i] * P1i1[j][l][i][m] / pop_size
        cov_with_F_save[naming_phenotype + "_covFw1i1"] = covFw1i1
        #  end of that part

        for i in range(pop_size):
            for j in range(sites):
                for k in range(sites):
                    if k != j:
                        for l in range(dm):
                            for m in range(dm):
                                covFw2[j][k][l][m] += (
                                    F[i] * P2[j][k][i][l][m] / pop_size
                                )
                                covFw2D[j][k][l][m] += (
                                    F[i] * P2D[j][k][i][l][m] / pop_size
                                )
                for n in range(sites):
                    for o in range(sites):
                        if n != o:
                            for p in range(dm):
                                for q in range(dm):
                                    covFw2i2[j][k][n][o][p][q] += (
                                        F[i] * P2i2[j][k][n][o][i][p][q] / pop_size
                                    )
        cov_with_F_save[naming_phenotype + "_covFw2"] = covFw2
        cov_with_F_save[naming_phenotype + "_covFw2D"] = covFw2D
        cov_with_F_save[naming_phenotype + "_covFw2i2"] = covFw2i2

        # Need for third degree polynomial
        # for i in range_popsize:
        #     for j in range_dm:
        #         for k in range_dm:
        #             for l in range_dm:
        #                 covFw3[j][k][l] += F[i] * P3[i][j][k][l] / pop_size

        # Covariance of the trait with each element of the second order phenotype.
        # note: We can nOT use the 'dot' operator here because
        # we do not have a matrix times a vector.
        covFPP = sum([F[i] * PP12[i] / pop_size for i in range_popsize])
        cov_with_F_save[naming_phenotype + "_covFPP"] = covFPP
    # Regressions of the trait on each element of the first order
    # phenotype vectors.
    if poly_order == "first":
        for j, i in itertools.product(range_sites, range_dm):
            if var[j][i] > 0.0000000001:
                rFon1[j][i] = covFw1[j][i] / var[j][i]
            else:
                rFon1[j][i] = 0
            if varP1D[j][i] > 0.0000000001:
                rFon1D[j][i] = covFw1D[j][i] / varP1D[j][i]
            else:
                rFon1D[j][i] = 0

        for i, j in itertools.product(range_sites, range_sites):
            if j != i:
                for k in range_dm:
                    if varP1i1[i][j][k] > 0.00000000001:
                        rFon1i1[i][j][k] = covFw1i1[i][j][k] / varP1i1[i][j][k]
                    else:
                        rFon1i1[i][j][k] = 0
        # Contribution of site 1 for each individual.
        # This is the regression of the trait on site 1 times (inner product)
        # the individual's site 1 vector that has been orthogonalized within
        # the vector.
        # test_val = sr.inner_general(rFon1[0],Pa[0][i])
        Fon1 = [sr.inner_general(rFon1[0], Pa[0][i]) for i in range_popsize]

        # Contribution of site 2 independent of 1 for each individual.
        Fon2i1 = [
            sr.inner_general(rFon1i1[1][0], Pa1i1[1][0][i]) for i in range_popsize
        ]

    # Regressions of the trait on each element of the second order
    # phenotype matrices.
    if poly_order == "second":
        # this part is from first order
        for j, i in itertools.product(range_sites, range_dm):
            if var[j][i] > 0.0000000001:
                rFon1[j][i] = covFw1[j][i] / var[j][i]
            else:
                rFon1[j][i] = 0
            if varP1D[j][i] > 0.0000000001:
                rFon1D[j][i] = covFw1D[j][i] / varP1D[j][i]
            else:
                rFon1D[j][i] = 0

        for i, j in itertools.product(range_sites, range_sites):
            if j != i:
                for k in range_dm:
                    if varP1i1[i][j][k] > 0.00000000001:
                        rFon1i1[i][j][k] = covFw1i1[i][j][k] / varP1i1[i][j][k]
                    else:
                        rFon1i1[i][j][k] = 0
        # Contribution of site 1 for each individual.
        # This is the regression of the trait on site 1 times (inner product)
        # the individual's site 1 vector that has been orthogonalized within
        # the vector.
        # test_val = sr.inner_general(rFon1[0],Pa[0][i])
        Fon1 = [sr.inner_general(rFon1[0], Pa[0][i]) for i in range_popsize]

        # Contribution of site 2 independent of 1 for each individual.
        Fon2i1 = [
            sr.inner_general(rFon1i1[1][0], Pa1i1[1][0][i]) for i in range_popsize
        ]
        # end of that section from first order

        for i, j in itertools.product(range_sites, range_sites):
            if j != i:
                for k, l in itertools.product(range_dm, range_dm):
                    if var2[i][j][k][l] > 0.00000000001:
                        rFon2[i][j][k][l] = covFw2[i][j][k][l] / var2[i][j][k][l]
                    else:
                        rFon2[i][j][k][l] = 0
                    if var2D[i][j][k][l] > 0.00000000001:
                        rFon2D[i][j][k][l] = covFw2D[i][j][k][l] / var2D[i][j][k][l]
                    else:
                        rFon2D[i][j][k][l] = 0

        for i, j in itertools.product(range_sites, range_sites):
            if j != i:
                for k, l in itertools.product(range_sites, range_sites):
                    if l != k:
                        for m, n in itertools.product(range_dm, range_dm):
                            if var2i2[i][j][k][l][m][n] > 0.0000000001:
                                numerator = covFw2i2[i][j][k][l][m][n]
                                denominator = var2i2[i][j][k][l][m][n]
                                rFon2i2[i][j][k][l][m][n] = numerator / denominator
                            else:
                                rFon2i2[i][j][k][l][m][n] = 0
        # Need when doing 3rd order
        # for i in range_dm:
        #     for j in range_dm:
        #         for k in range_dm:
        #             if var3[i][j][k] > 0.0000000001:
        #                 rFon3[i][j][k] = covFw3[i][j][k] / var3[i][j][k]
        #             else:
        #                 rFon3[i][j][k] = 0

        # print("rFon2i1"+str(rFon2i1))
        #
        # # Regressions of the trait on each element
        # of the second order phenotype matrix.
        # nucleotide1, nucleotide2
        for i, j in itertools.product(range_dm, range_dm):
            if var12[i][j] > 0.0000000001:
                rFon12[i][j] = covFPP[i][j] / var12[i][j]
            else:
                rFon12[i] = 0
        # # Contribution of the second order phenotype for each individual.
        Fon12 = [sr.inner_general(rFon2[0][1], P2a[0][1][i]) for i in range_popsize]
        Fon12 = [
            0 if np.fabs(Fon12[i]) < 0.0000000000001 else Fon12[i]
            for i in range_popsize
        ]

        # ----------Calculating the expected trait value for each individual
        # ----------given it's phenotype and the regressions calculated
        # -----------above (to check  whether or not everything works).

        Fest = [Fm + Fon1[i] + Fon2i1[i] + Fon12[i] for i in range_popsize]
        Fest = [
            0 if np.fabs(Fest[i]) < 0.0000000000001 else Fest[i] for i in range_popsize
        ]
    # contribution of third order phenotype for each individual......
    # for i in range_popsize:
    #     Fon3[i] = sr.inner_general(rFon3[0], P3a)
    # Ignoring very small values that would be due to roundoff error.
    # Change or delete this for a large data set.
    Fon1 = [0 if np.fabs(Fon1[i]) < 0.0000000000001 else Fon1[i] for i in range_popsize]
    Fon2i1 = [
        0 if np.fabs(Fon2i1[i]) < 0.0000000000001 else Fon2i1[i] for i in range_popsize
    ]

    output_npz_file = os.path.join(out_dir, naming_phenotype + "_covs_with_F.npz")
    print("Saving to {}".format(output_npz_file))
    np.savez_compressed(output_npz_file, **cov_with_F_save)

    # -----------------------Listing the main results------------------
    regression_results = {}
    print("Regression of trait on site 1")
    print(rFon1)
    print("Regression on 1st order polynomial - orthogonalized within - rFon1D")
    print(rFon1D)
    print("Regression of trait on site 2 independent of 1")
    print(rFon2i1)
    regression_results[naming_phenotype + "_rFon1"] = rFon1
    regression_results[naming_phenotype + "_rFon1D"] = rFon1D
    regression_results[naming_phenotype + "_rFon2i1"] = rFon2i1
    regression_results[naming_phenotype + "_Fon1"] = Fon1
    regression_results[naming_phenotype + "_Fon2i1"] = Fon2i1
    regression_results[naming_phenotype + "_Fest"] = Fest

    print("computed rFon1")
    # rFon1D is needed as we're working with 3 sites
    print("computed rFon1D")

    if poly_order == "second":
        print("Regression of trait on site 2")
        print(rFon2)
        print("Regression on 2nd order polynomial - orthogonalized within - rFon2D")
        print(rFon2D)
        print("Regression on (site 1)x(site 2), independent of first order")
        print(rFon12)
        # we need rFon2D when doing up to 3rd order and just rFon2
        # when doing up to 2nd order
        # rFon2D is needed as we're working with 3 sites
        regression_results[naming_phenotype + "_rFon2"] = rFon2
        regression_results[naming_phenotype + "_rFon2D"] = rFon2D
        regression_results[naming_phenotype + "_Fon12"] = Fon12
        regression_results[naming_phenotype + "_rFon12"] = rFon12
        print("computed rFon2")
        # rFon2D is needed as we're working with 3 sites
        print("computed rFon2D")

    regression_npz_file = os.path.join(out_dir, naming_phenotype + "_regressions.npz")
    print("Saving regression results to to {}".format(regression_npz_file))
    np.savez_compressed(regression_npz_file, **regression_results)

    print("Trait values estimated from regressions")
    print(Fest)
    print("--- %s seconds ---" % (time.time() - start_time))


@click.command(help="program to compute orthogonal polynomials up to 2nd order")  # noqa
@click.argument("filename", type=str)  # noqa
@click.option(
    "--pop_size", default=1, help="Population size or number of sequences"
)  # noqa
@click.option(
    "--dm",
    default=4,
    help="dimension of vector, e.g., this is =4 when input is DNA/RNA",
)  # noqa
@click.option(
    "--sites", default=2, help="number of sites in a sequence"
)  # starting off with two sites to run full second order
@click.option(
    "--molecule", default="DNA", help="can provide DNA or amino acid sequence"
)
@click.option(
    "--pheno_file", type=str, help="phenotype text file corresponding to sequence data"
)
@click.option(
    "--poly_order", default="first", help="can do first and second order so far"
)
@click.option(
    "--precomputed", default=False, help="if true, then saved results are used"
)
@click.option(
    "--out_dir", help="directory to save output/debug files to", type=str
)  # noqa
# @click.argument('pheno_file', type=click.File('rb'))
def cli(
    filename,
    pop_size,
    dm,
    sites,
    molecule,
    pheno_file,
    poly_order,
    precomputed,
    out_dir,
):
    orthogonal_polynomial(
        filename,
        pheno_file,
        molecule,
        sites,
        dm,
        pop_size,
        poly_order,
        precomputed,
        out_dir,
    )
