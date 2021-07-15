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
		--out_dir results \
		--custom_alphbt None/


unit-tests:
	pytest -v ortho_seq_code

coverage:
	coverage run --source ortho_seq_code --omit="*/test*" --module py.test
	coverage report --show-missing
