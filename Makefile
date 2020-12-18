SEQS=ortho_seq_code/tests/data/protein/protein_seqs_nopad.txt
PHENO=ortho_seq_code/tests/data/protein/protein_pheno_nopad.txt

run-cli:
	ortho_seq orthogonal-polynomial \
		${SEQS} \
		--pheno_file ${PHENO}  \
		--molecule protein \
		--sites 6 \
		--dm 20 \
		--pop_size 6  \
		--poly_order first \
		--out_dir results/


unit-tests:
	pytest ortho_seq_code --pdb -x

test:
	py.test

coverage:
	coverage run --source sencha --omit="*/test*" --module py.test
	coverage report --show-missing
