SEQS=ortho_seq_code/tests/data/protein/protein_seqs_nopad.txt
PHENO=ortho_seq_code/tests/data/protein/protein_pheno_nopad.txt

run-cli:
	ortho_seq orthogonal-polynomial \
		${SEQS} \
		--pheno_file ${PHENO}  \
		--molecule protein \
		--poly_order first \
		--out_dir results/ \
		--alphbt_input None \
		--min_pct 72


unit-tests:
	pytest -v ortho_seq_code

coverage:
	coverage run --source ortho_seq_code --omit="*/test*" --module py.test
	coverage report --show-missing
