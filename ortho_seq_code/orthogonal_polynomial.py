import numpy as np
from numpy.linalg import *
import time
import os
import ortho_seq_code.sr as sr
import click


@click.command(help='program to compute orthogonal polynomials up to 2nd order') # noqa
@click.option('--pop_size', default=1, help='Population size or number of sequences') # noqa
@click.option('--dm', default=4, help='dimension of vector, e.g., this is =4 when input is DNA/RNA') # noqa
@click.option('--sites', default=2, help='number of sites in a sequence') #starting off with two sites to run full second order
@click.option('--molecule', default='DNA', help='can provide DNA or amino acid sequence')
@click.option('--pheno_file', type=str, help="phenotype text file corresponding to sequence data")
@click.option('--poly_order', default='first', help='can do first and second order so far')
@click.option('--precomputed', default='False', help='if true, then saved results are used')
@click.option('--out_dir', help="directory to save output/debug files to", type=str) # noqa
@click.argument('filename', type=str) # noqa
#@click.argument('pheno_file', type=click.File('rb'))
def orthogonal_polynomial(filename, pheno_file, molecule, sites, dm, pop_size, poly_order, precomputed, out_dir):
    """Program to compute orthogonal polynomials up to 2nd order"""
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    else:
        print(
            "Path {} already exists, might be overwriting data".format(
                out_dir))

    start_time = time.time()
    with open(filename) as f:
        seq = f.readlines()
    global i
    # file containing trait values that will be mapped to sequence
    # vectors that must be the same size as F
    with open(pheno_file) as f2:
         phenotype = f2.readlines()
    #f2 = open(pheno_file, "rb")
    print(phenotype)
    F = np.genfromtxt(phenotype)  # this needs to stay this way!
    Fest = np.genfromtxt(phenotype) # this needs to stay this way!
    Fon1 = np.genfromtxt(phenotype) # this needs to stay this way!
    Fon2i1 = np.genfromtxt(phenotype) # this needs to stay this way!
    Fon12 = np.genfromtxt(phenotype) # this needs to stay this way!
    print(Fest)
    for i in range(pop_size):
        Fest[i] = 0
        Fon1[i] = 0
        Fon2i1[i] = 0
        Fon12[i] = 0
    np.set_printoptions(precision=10)
    # ----Initializing various terms that we will use.--------------
    # 3 sites, each a dm dim vector, in n individuals
    # nOTE: For application to Amino Acid sequences, increase
    # the size of the arrays accordingly.
    phi = np.array(
        [[[0.0 for k in range(dm)]
            for i in range(pop_size)] for j in range(sites)])
    mean = np.array([[0.0 for z in range(dm)] for i in range(sites)])
    var = np.array([[0.0 for z in range(dm)] for i in range(sites)])
    phi2 = np.array(
        [[[[[0.0 for k in range(dm)] for i in range(dm)]
            for j in range(pop_size)] for l in range(sites)]
            for m in range(sites)])
    phi2m = np.array(
        [[[[0.0 for k in range(dm)] for i in range(dm)]
            for l in range(sites)]
            for m in range(sites)])
    P = np.array(
        [[[0.0 for z in range(dm)] for j in range(pop_size)]
            for i in range(sites)])
    cov = np.array(
        [[[[0.0 for z in range(dm)] for i in range(dm)]
            for j in range(sites)] for k in range(sites)])
    # ------------Converting letters to vectors---------------
    # phi[individual][site][state]. phi[i][j] = vector for site j
    # in individual i.
    DNA_ALPHABETS = ['A', 'C', 'G', 'T']
    if dm == 4:
        if molecule == 'DNA':
            for dna_alphabet_index in range(len(DNA_ALPHABETS)):
                for i in range(pop_size):
                    for j in range(sites):
                        if seq[i][j] == DNA_ALPHABETS[dna_alphabet_index]:
                            phi[j][i][dna_alphabet_index] = 1.0

    DNA_ALPHABETS_n = ['A', 'C', 'G', 'T', 'n']
    if dm == 5:
        if molecule == 'DNA_n':
            for dna_alphabet_index in range(len(DNA_ALPHABETS_n)):
                for i in range(pop_size):
                    for j in range(sites):
                        if seq[i][j] == DNA_ALPHABETS_n[dna_alphabet_index]:
                            phi[j][i][dna_alphabet_index] = 1.0


    PROTEIN_ALPHABETS = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G','H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S','T', 'W', 'Y', 'V']
    if dm == 20:
        if molecule == 'protein':
            for protein_alphabet_index in range(len(PROTEIN_ALPHABETS)):
                for i in range(pop_size):
                    for j in range(sites):
                        if seq[i][j] == PROTEIN_ALPHABETS[protein_alphabet_index]:
                            phi[j][i][protein_alphabet_index] = 1.0

    PROTEIN_ALPHABETS_n = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G','H', 'I',
    'L', 'K', 'M', 'F', 'P', 'S','T', 'W', 'Y', 'V', 'n']
    if dm == 21:
        if molecule == 'protein_n':
            for protein_alphabet_index in range(len(PROTEIN_ALPHABETS_n)):
                for i in range(pop_size):
                    for j in range(sites):
                        if seq[i][j] == PROTEIN_ALPHABETS_n[protein_alphabet_index]:
                            phi[j][i][protein_alphabet_index] = 1.0

    PROTEIN_ALPHABETS_polar =  ['R', 'N', 'D', 'C', 'E', 'Q','H',
   'K','S','T', 'Y']
    PROTEIN_ALPHABETS_nonpolar = ['A', 'G','I', 'L', 'M', 'F', 'P', 'W', 'V']
    if dm == 3:
       if molecule == 'protein_pnp':
           for protein_alphabet_index_polar in range(len(PROTEIN_ALPHABETS_polar)):
               for protein_alphabet_index_nonpolar in range(len(PROTEIN_ALPHABETS_nonpolar)):
                   for i in range(pop_size):
                       for j in range(sites):
                           if seq[i][j] == PROTEIN_ALPHABETS_polar[protein_alphabet_index_polar]:
                               phi[j][i][0] = 1.0
                           if seq[i][j] == PROTEIN_ALPHABETS_nonpolar[protein_alphabet_index_nonpolar]:
                               phi[j][i][1] = 1.0
                           if seq[i][j] == 'n':
                               phi[j][i][2] = 1.0
    # if dm == 3:
    #     if molecule == 'protein_pnp':
    #         for i in range(pop_size):  # individual
    #             for j in range(sites):
    #                 if seq[i][j] == 'R':
    #                     phi[j][i][0] = 1.0
    #                 if seq[i][j] == 'N':
    #                     phi[j][i][0] = 1.0
    #                 if seq[i][j] == 'D':
    #                     phi[j][i][0] = 1.0
    #                 if seq[i][j] == 'C':
    #                     phi[j][i][0] = 1.0
    #                 if seq[i][j] == 'E':
    #                     phi[j][i][0] = 1.0
    #                 if seq[i][j] == 'Q':
    #                     phi[j][i][0] = 1.0
    #                 if seq[i][j] == 'H':
    #                     phi[j][i][0] = 1.0
    #                 if seq[i][j] == 'K':
    #                     phi[j][i][0] = 1.0
    #                 if seq[i][j] == 'S':
    #                     phi[j][i][0] = 1.0
    #                 if seq[i][j] == 'T':
    #                     phi[j][i][0] = 1.0
    #                 if seq[i][j] == 'Y':
    #                     phi[j][i][0] = 1.0
    #                 if seq[i][j] == 'A':
    #                     phi[j][i][1] = 1.0
    #                 if seq[i][j] == 'G':
    #                     phi[j][i][1] = 1.0
    #                 if seq[i][j] == 'I':
    #                     phi[j][i][1] = 1.0
    #                 if seq[i][j] == 'L':
    #                     phi[j][i][1] = 1.0
    #                 if seq[i][j] == 'M':
    #                     phi[j][i][1] = 1.0
    #                 if seq[i][j] == 'F':
    #                     phi[j][i][1] = 1.0
    #                 if seq[i][j] == 'P':
    #                     phi[j][i][1] = 1.0
    #                 if seq[i][j] == 'W':
    #                     phi[j][i][1] = 1.0
    #                 if seq[i][j] == 'V':
    #                     phi[j][i][1] = 1.0
    #                 if seq[i][j] == 'n':
    #                       phi[j][i][2] = 1.0

    # keep in alpha order
    # ---------------------------------First order terms ----------------------
    # calculate mean vectors first
    naming = os.path.basename(f.name)
    if poly_order == 'first':
        if precomputed == 'True':
            mean = np.load(os.path.join(out_dir, naming + str('_mean.npy')))
        else:
            for i in range(pop_size):
                for j in range(sites):
                    mean[j] += phi[j][i] / pop_size
            np.save(os.path.join(out_dir, naming + str('_mean')), mean)
            #  to show progress, can do something much more efficient/elegant
            print("computed mean")

        if precomputed == 'True':
            P = np.load(os.path.join(out_dir, naming + str('_P.npy')))
        else:
            for j in range(sites):  # site
                for i in range(0, pop_size):  # indiv
                    P[j][i] = phi[j][i] - mean[j]
            naming = os.path.basename(f.name)
            np.save(os.path.join(out_dir, naming + str('_P')), P)

        # var[site][nucleotide]
        if precomputed == 'True':
            var = np.load(os.path.join(out_dir, naming + str('_var.npy')))
        else:
            for k in range(sites):
                for i in range(0, dm):  # nucleotide
                    for j in range(0, pop_size):  # individual
                        var[k][i] += ((P[k][j][i]) ** 2) / pop_size
            print("computed variance")
            naming = os.path.basename(f.name)
            np.save(os.path.join(out_dir, naming + str('_var')), var)

        # Covariances between nucleotides at sites i and j
        # this is a matrix
        # the cov matrix for the two sites is just the mean,
        # across all individuals, of the outer product of P1 and P2
        # #P2 is site 2 with means subtracted out
        if precomputed == 'True':
            cov = np.load(os.path.join(out_dir, naming + str('_cov.npy')))
        else:
            for j in range(sites):
                for k in range(sites):
                    for i in range(pop_size):
                        cov[j][k] += sr.outer_general(P[j][i], P[k][i]) / pop_size
            print("computed covariance")
            naming = os.path.basename(f.name)
            np.save(os.path.join(out_dir, naming + str('_cov')), cov)

        Pa = np.array(
            [[[0.0 for z in range(dm)]
                for j in range(pop_size)] for i in range(sites)])
        reg11 = np.array(
            [[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(sites)] for k in range(sites)])
        # Regression of site k on site l
        # reg11[2][1] is a matrix (hence two more indices)
        # containing the regressions of
        # each element of the site 2 vector on each element of the site 1 vector
        if precomputed == 'True':
            reg11 = np.load(os.path.join(out_dir, naming + str('_reg11.npy')))
        else:
            for k in range(sites):
                for l in range(sites):
                    for i in range(0, dm):
                        for j in range(0, dm):
                            if var[l][j] > 0.0000000000001:
                                reg11[k][l][i][j] = cov[k][l][i][j] / var[l][j]
                            else:
                                reg11[k][l][i][j] = 0
            print("computed reg11")
            naming = os.path.basename(f.name)
            np.save(os.path.join(out_dir, naming + str('_reg11')), reg11)

        # # # # # First order terms with zeros except for the value that is
        # # # # # present (i.e. orthogonalized within each vector)
        if precomputed == 'True':
            Pa = np.load(os.path.join(out_dir, naming + str('_Pa.npy')))
        else:
            for k in range(sites):  # site
                for i in range(pop_size):  # indiv
                    Pa[k][i] = sr.inner_general(
                        sr.outer_general(phi[k][i], phi[k][i]), P[k][i])
            print("computed Pa: first order orthogonalized within each vector")
            naming = os.path.basename(f.name)
            np.save(os.path.join(out_dir, naming + str('_Pa')), Pa)

        if precomputed == 'True':
            P1i1 = np.load(os.path.join(out_dir, naming + str('_P1i1.npy')))
        else:
        # Site j orthogonalized wrt site k
        # P1i1[j][k][i] =  first order phi of site j independent of
        # site k for individual i
            P1i1 = np.array(
                [[[[0.0 for z in range(dm)] for i in range(pop_size)]
                    for j in range(sites)]
                    for k in range(sites)])
            for i in range(0, pop_size):  # Individuals
                for j in range(sites):
                    for k in range(sites):
                        if k != j:
                            P1i1[j][k][i] = P[j][i] - sr.inner_general(
                                reg11[j][k], Pa[k][i])
            print("computed P1i1")
            naming = os.path.basename(f.name)
            np.save(os.path.join(out_dir, naming + str('_P1i1')), P1i1)
        P2i1 = P1i1[1][0]

        if precomputed == 'True':
            varP1i1 = np.load(os.path.join(out_dir, naming + str('_varP1i1.npy')))
        else:
        # # # # Variance in P2i1
            varP1i1 = np.array(
                [[[0.0 for z in range(dm)]
                    for i in range(sites)] for j in range(sites)])
            for i in range(0, pop_size):  # individuals
                for j in range(sites):
                    for k in range(sites):
                        if k != j:
                            varP1i1[j][k] += (P1i1[j][k][i] ** 2) / pop_size
            print("computed varP1i1")
            naming = os.path.basename(f.name)
            np.save(os.path.join(out_dir, naming + str('_varP1i1')), varP1i1)

        # # # # cov11i1[j][k][l] = cov between site j and (site k independent of l)
        if precomputed == 'True':
            cov11i1 = np.load(os.path.join(out_dir, naming + str('_cov11i1.npy')))
        else:
            cov11i1 = np.array(
                [[[[[0.0 for z in range(dm)] for i in range(dm)]
                    for j in range(sites)] for k in range(sites)] for l in
                 range(sites)])
            for j in range(sites):
                for k in range(sites):
                    for l in range(sites):
                        for i in range(pop_size):
                            cov11i1[j][k][l] += \
                                sr.outer_general(P[j][i], P1i1[k][l][i]) / pop_size
            print("computed cov11i1")
            naming = os.path.basename(f.name)
            np.save(os.path.join(out_dir, naming + str('_cov11i1')), cov11i1)

        if precomputed == 'True':
            reg11i1 = np.load(os.path.join(out_dir, naming + str('_reg11i1.npy')))
        else:
        # # # # # regression of site j on (site k independent of l)
            reg11i1 = np.array(
                [[[[[0.0 for z in range(dm)] for i in range(dm)]
                    for j in range(sites)] for k in range(sites)]
                    for l in range(sites)])
            for k in range(sites):
                for l in range(sites):
                    for m in range(sites):
                        for i in range(0, dm):
                            for j in range(0, dm):
                                if varP1i1[l][m][j] > 0.0000000000001:
                                    reg11i1[k][l][m][i][j] = \
                                        cov11i1[k][l][m][i][j] / varP1i1[l][m][j]
                                else:
                                    reg11i1[k][l][m][i][j] = 0
            print("computed reg11i1")
            naming = os.path.basename(f.name)
            np.save(os.path.join(out_dir, naming + str('_reg11i1')), reg11i1)

        if precomputed == 'True':
            Pa1i1 = np.load(os.path.join(out_dir, naming + str('_Pa1i1.npy')))
        else:
        # # # # # same as P1i1, except with all elements = 0 except the one present
            Pa1i1 = np.array(
                [[[[0.0 for z in range(dm)]
                    for i in range(pop_size)] for j in range(sites)]
                    for k in range(sites)])
            for i in range(0, pop_size):  # indiv
                for j in range(sites):
                    for k in range(sites):
                        if k != j:
                            Pa1i1[j][k][i] = sr.inner_general(
                                sr.outer_general(phi[j][i], phi[j][i]), P1i1[j][k][i])
            print("computed Pa1i1")
            naming = os.path.basename(f.name)
            np.save(os.path.join(out_dir, naming + str('_Pa1i1')), Pa1i1)
        # # # # # P1D[j][i] = first order poly of site j independent of all other
        # sites, for individual i
        if precomputed == 'True':
            P1D = np.load(os.path.join(out_dir, naming + str('_P1D.npy')))
        else:
            P1D = np.array(
                [[[0.0 for z in range(dm)]
                    for i in range(pop_size)] for j in range(sites)])
            for i in range(pop_size):  # indiv
                for j in range(sites):
                    for k in range(sites):
                        if k != j:
                            for l in range(sites):
                                if l != k & l != j:
                                    P1D[j][i] = P[j][i] - sr.inner_general(reg11i1[j][k][l], Pa1i1[k][l][i]) - \
                                                sr.inner_general(reg11[j][l], Pa[l][i])
            print("computed P1D")
            naming = os.path.basename(f.name)
            np.save(os.path.join(out_dir, naming + str('_P1D')), P1D)

        if precomputed == 'True':
            varP1D = np.load(os.path.join(out_dir, naming + str('_varP1D.npy')))
        else:
            # # # #variance in P1D
            varP1D = np.array([[0.0 for z in range(dm)] for i in range(sites)])
            for k in range(sites):
                for i in range(0, dm):  # nucleotide
                    for j in range(0, pop_size):  # individual
                        varP1D[k][i] += ((P1D[k][j][i]) ** 2) / pop_size
            print("computed varP1D")
            naming = os.path.basename(f.name)
            np.save(os.path.join(out_dir, naming + str('_varP1D')), varP1D)
    # Pa2i1 = Pa1i1[1][0]
    # varP2i1 = varP1i1[1][0]
    # # # #-------------------------------------------------------------
    # # # # ------------------------Second Order Terms -----------------
    # # # #-------------------------------------------------------------

    # # # # Second order phenotypes.
    if poly_order == 'second':
        #this is also written under first order
        # calculate mean vectors
        for i in range(pop_size):
            for j in range(sites):
                mean[j] += phi[j][i] / pop_size
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_mean')), mean)
        #  to show progress, can do something much more efficient/elegant
        print("computed mean")

        for j in range(sites):  # site
            for i in range(0, pop_size):  # indiv
                P[j][i] = phi[j][i] - mean[j]
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_P')), P)

        # var[site][nucleotide]
        for k in range(sites):
            for i in range(0, dm):  # nucleotide
                for j in range(0, pop_size):  # individual
                    var[k][i] += ((P[k][j][i]) ** 2) / pop_size
        print("computed variance")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_var')), var)

        # Covariances between nucleotides at sites i and j
        # this is a matrix
        # the cov matrix for the two sites is just the mean,
        # across all individuals, of the outer product of P1 and P2
        # #P2 is site 2 with means subtracted out
        for j in range(sites):
            for k in range(sites):
                for i in range(pop_size):
                    cov[j][k] += sr.outer_general(P[j][i], P[k][i]) / pop_size
        print("computed covariance")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_cov')), cov)

        Pa = np.array(
            [[[0.0 for z in range(dm)]
                for j in range(pop_size)] for i in range(sites)])
        reg11 = np.array(
            [[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(sites)] for k in range(sites)])
        # Regression of site k on site l
        # reg11[2][1] is a matrix (hence two more indices)
        # containing the regressions of
        # each element of the site 2 vector on each element of the site 1 vector
        for k in range(sites):
            for l in range(sites):
                for i in range(0, dm):
                    for j in range(0, dm):
                        if var[l][j] > 0.0000000000001:
                            reg11[k][l][i][j] = cov[k][l][i][j] / var[l][j]
                        else:
                            reg11[k][l][i][j] = 0
        print("computed reg11")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_reg11')), reg11)

        # # # # # First order terms with zeros except for the value that is
        # # # # # present (i.e. orthogonalized within each vector)
        for k in range(sites):  # site
            for i in range(pop_size):  # indiv
                Pa[k][i] = sr.inner_general(
                    sr.outer_general(phi[k][i], phi[k][i]), P[k][i])
        print("computed Pa: first order orthogonalized within each vector")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_Pa')), Pa)

        P1i1 = np.array(
            [[[[0.0 for z in range(dm)] for i in range(pop_size)]
                for j in range(sites)]
                for k in range(sites)])
        for i in range(0, pop_size):  # Individuals
            for j in range(sites):
                for k in range(sites):
                    if k != j:
                        P1i1[j][k][i] = P[j][i] - sr.inner_general(
                            reg11[j][k], Pa[k][i])
        P2i1 = P1i1[1][0]
        print("computed P1i1")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_P1i1')), P1i1)

        # # # # Variance in P2i1
        varP1i1 = np.array(
            [[[0.0 for z in range(dm)]
                for i in range(sites)] for j in range(sites)])
        for i in range(0, pop_size):  # individuals
            for j in range(sites):
                for k in range(sites):
                    if k != j:
                        varP1i1[j][k] += (P1i1[j][k][i] ** 2) / pop_size
        print("computed varP1i1")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_varP1i1')), varP1i1)
        # # # # cov11i1[j][k][l] = cov between site j and (site k independent of l)
        cov11i1 = np.array(
            [[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(sites)] for k in range(sites)] for l in
             range(sites)])
        for j in range(sites):
            for k in range(sites):
                for l in range(sites):
                    for i in range(pop_size):
                        cov11i1[j][k][l] += \
                            sr.outer_general(P[j][i], P1i1[k][l][i]) / pop_size
        print("computed cov11i1")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_cov11i1')), cov11i1)
        # # # # # regression of site j on (site k independent of l)
        reg11i1 = np.array(
            [[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(sites)] for k in range(sites)]
                for l in range(sites)])
        for k in range(sites):
            for l in range(sites):
                for m in range(sites):
                    for i in range(0, dm):
                        for j in range(0, dm):
                            if varP1i1[l][m][j] > 0.0000000000001:
                                reg11i1[k][l][m][i][j] = \
                                    cov11i1[k][l][m][i][j] / varP1i1[l][m][j]
                            else:
                                reg11i1[k][l][m][i][j] = 0
        print("computed reg11i1")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_reg11i1')), reg11i1)
        # # # # # same as P1i1, except with all elements = 0 except the one present
        Pa1i1 = np.array(
            [[[[0.0 for z in range(dm)]
                for i in range(pop_size)] for j in range(sites)]
                for k in range(sites)])
        for i in range(0, pop_size):  # indiv
            for j in range(sites):
                for k in range(sites):
                    if k != j:
                        Pa1i1[j][k][i] = sr.inner_general(
                            sr.outer_general(phi[j][i], phi[j][i]), P1i1[j][k][i])
        print("computed Pa1i1")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_Pa1i1')), Pa1i1)
        # # # # # P1D[j][i] = first order poly of site j independent of all other
        # sites, for individual i
        P1D = np.array(
            [[[0.0 for z in range(dm)]
                for i in range(pop_size)] for j in range(sites)])
        for i in range(pop_size):  # indiv
            for j in range(sites):
                for k in range(sites):
                    if k != j:
                        for l in range(sites):
                            if l != k & l != j:
                                P1D[j][i] = P[j][i] - sr.inner_general(reg11i1[j][k][l], Pa1i1[k][l][i]) - \
                                            sr.inner_general(reg11[j][l], Pa[l][i])
        print("computed P1D")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_P1D')), P1D)
        # # # #variance in P1D
        varP1D = np.array([[0.0 for z in range(dm)] for i in range(sites)])
        for k in range(sites):
            for i in range(0, dm):  # nucleotide
                for j in range(0, pop_size):  # individual
                    varP1D[k][i] += ((P1D[k][j][i]) ** 2) / pop_size
        print("computed varP1D")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_varP1D')), varP1D)
        ### end of that part from first order

        for k in range(pop_size):
            for i in range(sites):
                for j in range(sites):
                    if j != i:
                        phi2[i][j][k] = sr.outer_general(phi[i][k], phi[j][k])
                        phi2m[i][j] += phi2[i][j][k] / pop_size
        # phi12 = phi2[0][1]
        # phi12m = phi2m[0][1]
        print("computed phi2 and phi2m")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_phi2')), phi2)
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_phi2m')), phi2m)

        # Q12 contains the 2'nd order phenotypes with the means subtracted out
        Q2 = np.array(
            [[[[[0.0 for k in range(dm)] for i in range(dm)]
                for j in range(pop_size)] for l in range(sites)]
                for m in range(sites)])
        for i in range(pop_size):  # indiv
            for j in range(sites):
                for k in range(sites):
                    if k != j:
                        Q2[j][k][i] = phi2[j][k][i] - phi2m[j][k]
        # Q12[i] = phi12[i] - phi12m
        # Q12 = Q2[0][1]
        print("computed Q2")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_Q2')), Q2)
        # # # # Covariance between elements of the 2'nd order phenotype matrix and
        # # # # the 1'st order phenotype.
        cov2w1 = np.array(
            [[[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)]
                for k in range(sites)]
                for l in range(sites)] for m in range(sites)])
        for i in range(pop_size):
            for j in range(sites):
                for k in range(sites):
                    if k != j:
                        for l in range(sites):
                            cov2w1[j][k][l] += \
                                sr.outer_general(Q2[j][k][i], P[l][i]) / pop_size
        print("computed cov2w1")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_cov2w1')), cov2w1)
        # # # # Covariance of second order phenotype matrices with first
        # order phenotypes.
        cov2w1a = np.array(
            [[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)]
                for k in range(sites)] for l in range(sites)])
        cov2w1b = np.array(
            [[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(dm)] for k in range(sites)] for l in range(sites)])
        # will need this when doing third order
        # cov2w1c = array([[[[[0.0 for z in range(dm)] for i in range(dm)]
        # for j in range(dm)] for k in range(sites)] for l in range(sites)])

        for i in range(pop_size):
            for j in range(sites):
                for k in range(sites):
                    if k != j:
                        cov2w1a[j][k] += sr.outer_general(
                            Q2[j][k][i], P[0][i]) / pop_size
                        cov2w1b[j][k] += \
                            sr.outer_general(
                                Q2[j][k][i], P1i1[1][0][i]) / pop_size
                        # will need this when doing third order
                        # cov2w1c[j][k] += sr.outer_general(
                        # Q2[j][k][i], P1D[2][i]) / n
        print("computed cov2w1a,cov2w1b")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_cov2w1a')), cov2w1a)
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_cov2w1b')), cov2w1b)

        # regressions of second order phenotype matrices on first order phenotypes.
        r2on1a = np.array(
            [[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(dm)] for k in range(sites)] for l in range(sites)])
        r2on1b = np.array(
            [[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(dm)] for k in range(sites)] for l in range(sites)])
        # will need this when doing third order
        # r2on1c = array([[[[[0.0 for z in range(dm)] for i in range(dm)] for j
        # in range(dm)] for k in range(sites)] for l in range(sites)])

        for i in range(sites):
            for j in range(sites):
                if j != i:
                    for k in range(dm):
                        for l in range(dm):
                            for m in range(dm):
                                if var[0][m] > 0.0000000001:
                                    r2on1a[i][j][k][l][m] = \
                                        cov2w1a[i][j][k][l][m] / var[0][m]
                                else:
                                    r2on1a[i][j][k][l][m] = 0
                                if varP1i1[1][0][m] > 0.0000000001:
                                    r2on1b[i][j][k][l][m] = \
                                        cov2w1b[i][j][k][l][m] / varP1i1[1][0][m]
                                else:
                                    r2on1b[i][j][k][l][m] = 0
                                # Need for third degree polynomial
                                # if varP1D[2][m] > 0.0000000001:
                                #     r2on1c[i][j][k][l][m] = \
                                # cov2w1c[i][j][k][l][m] / varP1D[2][m]
                                # else:
                                #     r2on1c[i][j][k][l][m] = 0
        print("computed r2on1a, r2on1b")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_r2on1a')), r2on1a)
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_r2on1b')), r2on1b)

        # # # Second order polynomials
        P2 = np.array(
            [[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(pop_size)] for k in range(sites)]
                for l in range(sites)])
        # P1Da = np.array(
        #   [[[0.0 for z in range(dm)] for i in range(pop_size)]
        # for j in range(sites)])
        for i in range(sites):
            for j in range(sites):
                if j != i:
                    for k in range(pop_size):
                        P2[i][j][k] = Q2[i][j][k] - sr.inner_general(
                            r2on1a[i][j], Pa[0][k]) - sr.inner_general(
                            r2on1b[i][j], Pa1i1[1][0][k])
                            # - sr.inner_general(r2on1c[i][j], P1Da[2][k]) # noqa
        # r12on1 = r2on1a[0][1]
        # r12on2i1 = r2on1b[0][1]
        print("computed P2")
        PP12 = P2[0][1]
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_P2')), P2)
        # # # Second order terms with zeros except for the value that is
        # # # present (i.e. orthogonalized within each matrix)
        P2a = np.array(
            [[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(pop_size)] for k in range(sites)] for l in
                range(sites)])
        for h in range(0, pop_size):  # individual
            for i in range(sites):
                for j in range(sites):
                    if j != i:
                        P2a[i][j][h] = sr.inner_general(
                            sr.outer_general(
                                phi2[i][j][h], phi2[i][j][h]), P2[i][j][h])
        print("computed P2a")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_P2a')), P2a)
        # PPa12 = P2a[0][1]
        # # # # Covariances between second order phenotypes
        cov2w2 = np.array(
            [[[[[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(dm)] for k in range(dm)] for l in
                range(sites)] for m in range(sites)] for n in range(sites)]
                for p in range(sites)])
        for i in range(pop_size):
            for j in range(sites):
                for k in range(sites):
                    if k != j:
                        for l in range(sites):
                            for m in range(sites):
                                if m != l:
                                    cov2w2[j][k][l][m] += sr.outer_general(
                                        P2[j][k][i], P2[l][m][i]) / pop_size
        print("computed cov2w2")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_cov2w2')), cov2w2)
        # # Variances of second order phenotypes
        var2 = np.array(
            [[[[0.0 for z in range(dm)]
                for i in range(dm)] for j in range(sites)] for k in range(sites)])
        for i in range(sites):
            for j in range(sites):
                if j != i:
                    for k in range(dm):
                        for l in range(dm):
                            var2[i][j][k][l] = cov2w2[i][j][i][j][k][l][k][l]
        var12 = var2[0][1]
        print("computed var2")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_var2')), var2)
        # #Variances of second order phenotypes (this is another way of computing variances, by using the polynomial itself)
        # # for i in range(sites):
        # #     for j in range(sites):
        # #         if j >> i:
        # #             for k in range(pop_size):
        # #                 for l in range(dm):
        # #                     for m in range(dm):
        # #                         var2[i][j][l][m] += \
        # #   (P2[i][j][k][l][m]**2)/n - (P2m[i][j][l][m]**2)/n

        # # # # regressions of second order phenotypes on one another
        reg2on2 = np.array(
            [[[[[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(dm)] for k in range(dm)] for l in
                range(sites)] for m in range(sites)] for n in range(sites)]
                for p in range(sites)])
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
                                                    numerator = \
                                                        cov2w2[i][j][k][l][m][n][o][p]
                                                    denominator = var2[k][l][o][p]
                                                    reg2on2[i][j][k][l][m][n][o][p] = \
                                                        numerator / denominator
                                                else:
                                                    reg2on2[i][j][k][l][m][n][o][p] = 0
        print("computed reg2on2")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_reg2on2')), reg2on2)
        # # # # Second order phenotypes independent of one another
        P2i2 = np.array(
            [[[[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(pop_size)] for k in range(sites)] for l in
                range(sites)] for m in range(sites)] for n in range(sites)])
        P2i2a = np.array(
            [[[[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(pop_size)] for k in range(sites)] for l in
                range(sites)] for m in range(sites)] for n in range(sites)])
        for i in range(pop_size):
            for j in range(sites):
                for k in range(sites):
                    if k != j:
                        for l in range(sites):
                            for m in range(sites):
                                if m != l:
                                    P2i2[j][k][l][m][i] = \
                                        P2[j][k][i] - sr.inner_general(
                                            reg2on2[j][k][l][m], P2a[l][m][i])
                                    P2i2a[j][k][l][m][i] = sr.inner_general(
                                        sr.outer_general(
                                            phi2[j][k][i],
                                            phi2[j][k][i]), P2i2[j][k][l][m][i])
        print("computed P2i2, P2i2a")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_P2i2')), P2i2)
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_P2i2a')), P2i2a)
        # # # # cov of 2'nd order phi with another independent of the third
        cov2w2i2 = np.array(
            [[[[[[[[[[0.0 for z in range(dm)] for i in range(dm)] for j in range(
                dm)]
                for k in range(dm)] for l
                in range(sites)] for m in range(sites)] for n in range(sites)]
                for p in range(sites)] for q
                in range(sites)] for r in range(sites)])
        var2i2 = np.array(
            [[[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(sites)] for k in range(sites)] for l
                in range(sites)] for m in range(sites)])
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
                                                cov2w2i2[j][k][l][m][n][o] += \
                                                    sr.outer_general(
                                                        P2[j][k][i],
                                                        P2i2[l][m][n][o][i]) / pop_size
                                    for p in range(dm):
                                        for q in range(dm):
                                            var2i2[j][k][l][m][p][q] += (
                                                P2i2[j][k][l][m][i][p][q] ** 2) / pop_size
        print("computed cov2w2i2, var2i2")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_cov2w2i2')), cov2w2i2)
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_var2i2')), var2i2)
        # # # # regressions corresponding to the above
        reg2on2i2 = np.array(
            [[[[[[[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(dm)] for k in range(dm)] for
                l in range(sites)] for m in range(sites)]
                for n in range(sites)] for p in range(sites)] for
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
                                                                numerator = \
                                                                    cov2w2i2[i][j][k][l][k1][l1][m][n][o][p]
                                                                denominator = \
                                                                    var2i2[k][l][k1][l1][o][p]
                                                                reg2on2i2[i][j][k][l][k1][l1][m][n][o][p] = \
                                                                    numerator / denominator
                                                            else:
                                                                reg2on2i2[i][j][k][l][k1][l1][m][n][o][p] = 0
        print("computed reg2on2i2")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_reg2on2i2')), reg2on2i2)
        # # # # 2'nd order phi independent of all others
        P2D = np.array(
            [[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(pop_size)] for k in range(sites)] for l in range(sites)])
        P2Da = np.array(
            [[[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(pop_size)] for k in range(sites)] for l in
                range(sites)])
        for i in range(1):
            for j in range(sites):
                for k in range(sites):
                    if k != j:
                        for l in range(sites):
                            for m in range(sites):
                                if m != l and (
                                        l, m) != (j, k) and (l, m) != (k, j):
                                    for n in range(sites):
                                        for o in range(sites):
                                            if o != n and (n, o) != (l, m) and (n, o) != (m, l) and (n, o) != (j, k) and (n, o) != (k, j):
                                                P2D[j][k][i] = \
                                                    P2[j][k][i] - sr.inner_general(
                                                        reg2on2i2[j][k][l][m][n][o],
                                                        P2i2a[l][m][n][o][i]) - sr.inner_general(
                                                        reg2on2[j][k][n][o],
                                                        P2a[n][o][i])
                                                P2Da[j][k][i] = sr.inner_general(
                                                    sr.outer_general(
                                                        phi2[j][k][i],
                                                        phi2[j][k][i]),
                                                    P2D[j][k][i])
        print("computed P2D, P2Da")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_P2D')), P2D)
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_P2Da')), P2Da)
        # # # #variance in P2D
        var2D = np.array(
            [[[[0.0 for z in range(dm)] for i in range(dm)]
                for j in range(sites)] for k in range(sites)])
        for i in range(pop_size):
            for j in range(sites):
                for k in range(sites):
                    if k != j:
                        for l in range(dm):
                            for m in range(dm):
                                var2D[j][k][l][m] += (P2D[j][k][i][l][m] ** 2) / pop_size
        print("computed var2D")
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_var2D')), var2D)
    # ------------Projecting the trait values into the space of orthogonal ----
    # -------------polynomials ------------------------------------------------
    # initializing arrays
    cov1FP = np.array([[0.0 for z in range(dm)] for i in range(sites)])
    covFP = np.array([[0.0 for z in range(dm)] for i in range(sites)])
    covFPP = np.array([[0.0 for z in range(dm)] for i in range(dm)])
    # rFon1 = array([0.0 for z in range(dm)])
    # rFon2 = array([0.0 for z in range(dm)])
    rFon2i1 = np.array([0.0 for z in range(dm)])
    rFon12 = np.array([[0.0 for z in range(dm)] for i in range(dm)])
    covFw1 = np.array([[0.0 for z in range(dm)] for i in range(sites)])
    covFw1i1 = np.array(
        [[[0.0 for z in range(dm)] for i in range(sites)]
            for j in range(sites)])
    covFw1D = np.array([[0.0 for z in range(dm)] for i in range(sites)])
    rFon1 = np.array([[0.0 for z in range(dm)] for i in range(sites)])
    rFon1i1 = np.array(
        [[[0.0 for z in range(dm)] for i in range(sites)]
            for j in range(sites)])
    rFon1D = np.array([[0.0 for z in range(dm)] for i in range(sites)])
    covFw2 = np.array(
        [[[[0.0 for z in range(dm)] for i in range(dm)]
            for j in range(sites)] for k in range(sites)])
    covFw2i2 = np.array(
        [[[[[[0.0 for z in range(dm)]
            for i in range(dm)] for j in range(sites)]
            for k in range(sites)] for
            l in range(sites)] for m in range(sites)])
    covFw2D = np.array(
        [[[[0.0 for z in range(dm)] for i in range(dm)]
            for j in range(sites)] for k in range(sites)])
    rFon2 = np.array(
        [[[[0.0 for z in range(dm)] for i in range(dm)]
            for j in range(sites)] for k in range(sites)])
    rFon2i2 = np.array(
        [[[[[[0.0 for z in range(dm)] for i in range(dm)]
            for j in range(sites)] for k in range(sites)] for
            l in range(sites)] for m in range(sites)])
    rFon2D = np.array(
        [[[[0.0 for z in range(dm)] for i in range(dm)]
            for j in range(sites)] for k in range(sites)])

    # Need for third degree polynomial
    # covFw3 = np.array(
    #     [[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)])
    # rFon3 = np.array(
    #     [[[0.0 for z in range(dm)] for i in range(dm)] for j in range(dm)])
    # Calculating the mean trait value
    Fm = 0
    for i in range(pop_size):  # individuals
        Fm += F[i] / pop_size
    naming_phenotype = os.path.basename(f2.name)
    np.save(os.path.join(out_dir, naming_phenotype + str('_Fm')), Fm)
    # Covariances of the trait with each element of the 1'st order vectors.
    # We can use the 'dot' operator here to get the inner product of a
    # first and a second rank tensor (a vector and a matrix).
    if poly_order == 'first':
        covFP[0] = np.dot(F, P[0]) / pop_size  # for site 1
        cov1FP[1] = np.dot(F, P[1]) / pop_size
        covFP[1] = np.dot(F, P2i1) / pop_size  # for site 2 independent of 1
        naming = os.path.basename(f2.name)
        np.save(os.path.join(out_dir, naming_phenotype + str('_covFP[0]')), covFP[0])
        naming = os.path.basename(f2.name)
        np.save(os.path.join(out_dir, naming_phenotype + str('_cov1FP[1]')), cov1FP[1])
        naming = os.path.basename(f2.name)
        np.save(os.path.join(out_dir, naming_phenotype + str('_covFP[1]')), covFP[1])

        for i in range(pop_size):
            for j in range(sites):
                for k in range(dm):
                    covFw1[j][k] += F[i] * P[j][i][k] / pop_size
                    covFw1D[j][k] += F[i] * P1D[j][i][k] / pop_size
                for l in range(sites):
                    if l != j:
                        for m in range(dm):
                            covFw1i1[j][l][m] += F[i] * P1i1[j][l][i][m] / pop_size
        naming = os.path.basename(f2.name)
        np.save(os.path.join(out_dir, naming_phenotype + str('_covFw1i1')), covFw1i1)

    if poly_order == 'second':
        #this part is from first order
        covFP[0] = np.dot(F, P[0]) / pop_size  # for site 1
        cov1FP[1] = np.dot(F, P[1]) / pop_size
        covFP[1] = np.dot(F, P2i1) / pop_size  # for site 2 independent of 1
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_covFP[0]')), covFP[0])
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_cov1FP[1]')), cov1FP[1])
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_covFP[1]')), covFP[1])

        for i in range(pop_size):
            for j in range(sites):
                for k in range(dm):
                    covFw1[j][k] += F[i] * P[j][i][k] / pop_size
                    covFw1D[j][k] += F[i] * P1D[j][i][k] / pop_size
                for l in range(sites):
                    if l != j:
                        for m in range(dm):
                            covFw1i1[j][l][m] += F[i] * P1i1[j][l][i][m] / pop_size
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_covFw1i1')), covFw1i1)
        #end of that part

        for i in range(pop_size):
            for j in range(sites):
                for k in range(sites):
                    if k != j:
                        for l in range(dm):
                            for m in range(dm):
                                covFw2[j][k][l][m] += F[i] * P2[j][k][i][l][m] / pop_size
                                covFw2D[j][k][l][m] += \
                                    F[i] * P2D[j][k][i][l][m] / pop_size
                for n in range(sites):
                    for o in range(sites):
                        if n != o:
                            for p in range(dm):
                                for q in range(dm):
                                    covFw2i2[j][k][n][o][p][q] += \
                                        F[i] * P2i2[j][k][n][o][i][p][q] / pop_size
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_covFw2')), covFw2)
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_covFw2D')), covFw2D)
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_covFw2i2')), covFw2i2)

    # Need for third degree polynomial
    # for i in range(pop_size):
    #     for j in range(dm):
    #         for k in range(dm):
    #             for l in range(dm):
    #                 covFw3[j][k][l] += F[i] * P3[i][j][k][l] / pop_size

    # Covariance of the trait with each element of the second order phenotype.
    # note: We can nOT use the 'dot' operator here because
    # we do not have a matrix times a vector.
        for i in range(0, pop_size):  # indiv
            covFPP += (F[i] * PP12[i] / pop_size)
        naming = os.path.basename(f.name)
        np.save(os.path.join(out_dir, naming + str('_covFPP')), covFPP)
    # Regressions of the trait on each element of the first order
    # phenotype vectors.
    if poly_order == 'first':
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

        np.save(os.path.join(out_dir, naming_phenotype + str('_rFon1')), rFon1)
        print("computed rFon1")
        # rFon1D is needed as we're working with 3 sites
        np.save(os.path.join(out_dir, naming_phenotype + str('_rFon1D')), rFon1D)
        print("computed rFon1D")

        for i in range(sites):
            for j in range(sites):
                if j != i:
                    for k in range(dm):
                        if varP1i1[i][j][k] > 0.00000000001:
                            rFon1i1[i][j][k] = covFw1i1[i][j][k] / varP1i1[i][j][k]
                        else:
                            rFon1i1[i][j][k] = 0
        # Contribution of site 1 for each individual.
        # This is the regression of the trait on site 1 times (inner product)
        # the individual's site 1 vector that has been orthogonalized within
        # the vector.
        for i in range(pop_size):
            # test_val = sr.inner_general(rFon1[0],Pa[0][i])
            Fon1[i] = sr.inner_general(rFon1[0], Pa[0][i])

        # Contribution of site 2 independent of 1 for each individual.
        for i in range(pop_size):
            Fon2i1[i] = sr.inner_general(rFon1i1[1][0], Pa1i1[1][0][i])

    # Regressions of the trait on each element of the second order
    # phenotype matrices.
    if poly_order == 'second':
        #this part is from first order
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

        np.save(os.path.join(out_dir, naming + str('_rFon1')), rFon1)
        print("computed rFon1")
        # rFon1D is needed as we're working with 3 sites
        np.save(os.path.join(out_dir, naming + str('_rFon1D')), rFon1D)
        print("computed rFon1D")

        for i in range(sites):
            for j in range(sites):
                if j != i:
                    for k in range(dm):
                        if varP1i1[i][j][k] > 0.00000000001:
                            rFon1i1[i][j][k] = covFw1i1[i][j][k] / varP1i1[i][j][k]
                        else:
                            rFon1i1[i][j][k] = 0
        # Contribution of site 1 for each individual.
        # This is the regression of the trait on site 1 times (inner product)
        # the individual's site 1 vector that has been orthogonalized within
        # the vector.
        for i in range(pop_size):
            # test_val = sr.inner_general(rFon1[0],Pa[0][i])
            Fon1[i] = sr.inner_general(rFon1[0], Pa[0][i])

        # Contribution of site 2 independent of 1 for each individual.
        for i in range(pop_size):
            Fon2i1[i] = sr.inner_general(rFon1i1[1][0], Pa1i1[1][0][i])
        #end of that section from first order

        for i in range(sites):
            for j in range(sites):
                if j != i:
                    for k in range(dm):
                        for l in range(dm):
                            if var2[i][j][k][l] > 0.00000000001:
                                rFon2[i][j][k][l] = \
                                    covFw2[i][j][k][l] / var2[i][j][k][l]
                            else:
                                rFon2[i][j][k][l] = 0
                            if var2D[i][j][k][l] > 0.00000000001:
                                rFon2D[i][j][k][l] = \
                                    covFw2D[i][j][k][l] / var2D[i][j][k][l]
                            else:
                                rFon2D[i][j][k][l] = 0
        # we need rFon2D when doing up to 3rd order and just rFon2
        # when doing up to 2nd order
        # rFon2D is needed as we're working with 3 sites
        np.save(os.path.join(out_dir, naming + str('_rFon2')), rFon2)
        print("computed rFon2")
        # rFon2D is needed as we're working with 3 sites
        np.save(os.path.join(out_dir, naming + str('_rFon2D')), rFon2D)
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
                                            numerator = covFw2i2[i][j][k][l][m][n]
                                            denominator = var2i2[i][j][k][l][m][n]
                                            rFon2i2[i][j][k][l][m][n] = \
                                                numerator / denominator
                                        else:
                                            rFon2i2[i][j][k][l][m][n] = 0
    # Need when doing 3rd order
    # for i in range(dm):
    #     for j in range(dm):
    #         for k in range(dm):
    #             if var3[i][j][k] > 0.0000000001:
    #                 rFon3[i][j][k] = covFw3[i][j][k] / var3[i][j][k]
    #             else:
    #                 rFon3[i][j][k] = 0

    # print("rFon2i1"+str(rFon2i1))
    #
        # # Regressions of the trait on each element
        # of the second order phenotype matrix.
        for i in range(0, dm):  # nucleotide 1
            for j in range(0, dm):  # nucleotide 2
                if var12[i][j] > 0.0000000001:
                    rFon12[i][j] = covFPP[i][j] / var12[i][j]
                else:
                    rFon12[i] = 0
        # # Contribution of the second order phenotype for each individual.
        for i in range(pop_size):
            Fon12[i] = sr.inner_general(rFon2[0][1], P2a[0][1][i])

        for i in range(pop_size):
            if np.fabs(Fon12[i]) < 0.0000000000001:
                Fon12[i] = 0

        # ----------Calculating the expected trait value for each individual
        # ----------given it's phenotype and the regressions calculated
        # -----------above (to check  whether or not everything works).

        # for i in range(pop_size):  # indiv
    	#        Fest[i] = Fm + Fon1[i] + Fon2i1[i] + Fon12[i]
    	#        if abs(Fest[i]) < 0.0000000000001:  # avoiding roundoff error
    	#              Fest[i] = 0   	           # modify or remove for large datasets
        for i in range(pop_size):  # indiv
        	Fest[i] = Fm + Fon1[i] + Fon2i1[i] + Fon12[i]
        	if abs(Fest[i]) < 0.0000000000001:  # avoiding roundoff error
        		Fest[i] = 0   	           # modify or remove for large datasets
        np.save(os.path.join(out_dir, naming + str('_Fest')), Fest)
    # contribution of third order phenotype for each individual......
    # for i in range(pop_size):
    #     Fon3[i] = sr.inner_general(rFon3[0], P3a)
    # Ignoring very small values that would be due to roundoff error.
    # Change or delete this for a large data set.
    for i in range(pop_size):
        if np.fabs(Fon1[i]) < 0.0000000000001:
            Fon1[i] = 0
        if np.fabs(Fon2i1[i]) < 0.0000000000001:
            Fon2i1[i] = 0

    # -----------------------Listing the main results------------------
    print('Regression of trait on site 1')
    print(rFon1)
    print(
        'Regression on 1st order polynomial - orthogonalized within - rFon1D')
    print(rFon1D)
    print('Regression of trait on site 2 independent of 1')
    print(rFon2i1)

    if poly_order == 'second':
        print('Regression of trait on site 2')
        print(rFon2)
        print(
            'Regression on 2nd order polynomial - orthogonalized within - rFon2D')
        print(rFon2D)
        print('Regression on (site 1)x(site 2), independent of first order')
        print(rFon12)

    print('Trait values estimated from regressions')
    print(Fest)
    print("--- %s seconds ---" % (time.time() - start_time))


if __name__ == '__main__':
    orthogonal_polynomial()
