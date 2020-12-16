import numpy as np
import os

from click.testing import CliRunner

from ortho_seq_code.orthogonal_polynomial import orthogonal_polynomial
from ortho_seq_code.cli import cli


def output_file_templates(basename):

    return [
        f"{basename}_mean.npy",
        f"{basename}_P.npy",
        f"{basename}_var.npy",
        f"{basename}_cov.npy",
        f"{basename}_reg11.npy",
        f"{basename}_Pa.npy",
        f"{basename}_P1i1.npy",
        f"{basename}_varP1i1.npy",
        f"{basename}_varP1i1.npy",
        f"{basename}_cov11i1.npy",
        f"{basename}_reg11i1.npy",
        f"{basename}_Pa1i1.npy",
        f"{basename}_P1D.npy",
        f"{basename}_varP1D.npy", #end of building first order space

    ]

def output_file_2ndorder_templates(basename):

    return [

        f"{basename}_phi2.npy",
        f"{basename}_phi2m.npy",
        f"{basename}_Q2.npy",
        f"{basename}_cov2w1.npy",
        f"{basename}_cov2w1a.npy",
        f"{basename}_cov2w1b.npy",
        f"{basename}_r2on1a.npy",
        f"{basename}_r2on1b.npy",
        f"{basename}_P2.npy",
        f"{basename}_P2a.npy",
        f"{basename}_cov2w2.npy",
        f"{basename}_var2.npy",
        f"{basename}_reg2on2.npy",
        f"{basename}_P2i2.npy",
        f"{basename}_P2i2a.npy",
        f"{basename}_cov2w2i2.npy",
        f"{basename}_var2i2.npy",
        f"{basename}_reg2on2i2.npy",
        f"{basename}_P2D.npy",
        f"{basename}_P2Da.npy", #end of building second order space

    ]

def pheno_output_file_templates(basename_pheno):

    return [

        f"{basename_pheno}_Fm.npy", #mean phenotype value
        f"{basename_pheno}_covFP[0].npy",
        f"{basename_pheno}_cov1FP[1].npy",
        f"{basename_pheno}_covFP[1].npy",
        f"{basename_pheno}_covFw1i1.npy",
        f"{basename_pheno}_rFon1.npy",
        f"{basename_pheno}_rFon1D.npy", #end of projections of phenotypes onto first order

    ]

def pheno_output_file_2ndorder_templates(basename_pheno_2ndorder):

    return [

        f"{basename_pheno_2ndorder}_covFw2.npy", #start of projections of phenotypes onto second order
        f"{basename_pheno_2ndorder}_covFw2D.npy",
        f"{basename_pheno_2ndorder}_covFw2i2.npy",
        f"{basename_pheno_2ndorder}_covFPP.npy",
        f"{basename_pheno_2ndorder}_rFon2.npy",
        f"{basename_pheno_2ndorder}_rFon2D.npy", #end of projections of phenotypes onto first order
        f"{basename_pheno_2ndorder}_Fest.npy", #trait values estimated from regressions/projections
    ]

def test_cli(
        protein_seqs_no_padding,
        protein_pheno_no_padding
):
    molecule = 'protein'
    sites = 6
    dm = 20
    pop_size = 6
    poly_order = 'first'
    out_dir = '/tmp'

    runner = CliRunner()

    result = runner.invoke(
        cli, [
            protein_seqs_no_padding,
            "--pheno_file", protein_pheno_no_padding,
            "--molecule",  molecule,
            "--sites", sites,
            "--dm", dm,
            "--pop_size", pop_size,
            "--poly_order", poly_order,
            "--out_dir", out_dir,
        ]
    )

    assert result.exit_code == 0


def test_nucleotide_first_order(
        nucleotide_first_order_data_dir,
        nucleotide_first_order_expected_output_dir,
        nucleotide_params_first_order
):

    orthogonal_polynomial(
        *nucleotide_params_first_order
    )

    basename = os.path.basename(nucleotide_params_first_order.seqs_filename)
    basename_pheno = os.path.basename(nucleotide_params_first_order.pheno_filename)

    for t in output_file_templates(basename): #work related to building the polynomial space
        expected = np.load(os.path.join(nucleotide_first_order_expected_output_dir, t))
        actual = np.load(os.path.join(nucleotide_params_first_order.out_dir, t))
        assert expected.all() == actual.all()

    for s in pheno_output_file_templates(basename_pheno): #work related to projecting phenotypes onto the polynomial space we built
        expected = np.load(os.path.join(nucleotide_first_order_expected_output_dir, s))
        actual = np.load(os.path.join(nucleotide_params_first_order.out_dir, s))
        assert expected.all() == actual.all()

def test_nucleotide_second_order(
        nucleotide_second_order_data_dir,
        nucleotide_second_order_expected_output_dir,
        nucleotide_params_second_order
):

    orthogonal_polynomial(
        *nucleotide_params_second_order
    )

    basename = os.path.basename(nucleotide_params_second_order.seqs_filename)
    basename_pheno_2ndorder = os.path.basename(nucleotide_params_second_order.pheno_filename)

    for t in output_file_2ndorder_templates(basename):
        expected = np.load(os.path.join(nucleotide_second_order_expected_output_dir, t))
        actual = np.load(os.path.join(nucleotide_params_second_order.out_dir, t))
        assert expected.all() == actual.all()

    for s in pheno_output_file_2ndorder_templates(basename_pheno_2ndorder):
        expected = np.load(os.path.join(nucleotide_second_order_expected_output_dir, s))
        actual = np.load(os.path.join(nucleotide_params_second_order.out_dir, s))
        assert expected.all() == actual.all()


def test_protein_first_order(
        protein_data_dir,
        protein_nopad_expected_output_dir,
        protein_params_first_order
):

    orthogonal_polynomial(*protein_params_first_order)

    basename = os.path.basename(protein_params_first_order.seqs_filename)
    basename_pheno = os.path.basename(protein_params_first_order.pheno_filename)


    for t in output_file_templates(basename):
        expected = np.load(os.path.join(protein_nopad_expected_output_dir, t))
        actual = np.load(os.path.join(protein_params_first_order.out_dir, t))
        assert expected.all() == actual.all()

    for s in pheno_output_file_templates(basename_pheno):
        expected = np.load(os.path.join(protein_nopad_expected_output_dir, s))
        actual = np.load(os.path.join(protein_params_first_order.out_dir, s))
        assert expected.all() == actual.all()
