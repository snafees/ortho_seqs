converting sequence data to tensor-valued orthogonal polynomials
# to install 
python setup.py install

To use the code, run the below command::

```
ortho_seq orthogonal-polynomial ./ortho_seq_code/testdata/pho4_r4_s1_site_3-4.csv --phenotype ./ortho_seq_code/testdata/pho4_r4_s1_ddg.txt  --pop-size 1000 --out-dir ../results_ortho_seq_testing
```