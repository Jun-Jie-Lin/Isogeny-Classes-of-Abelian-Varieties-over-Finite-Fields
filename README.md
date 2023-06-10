# Isogeny-Classes-of-Abelian-Varieties-over-Finite-Fields

This repository contains the code used to check the results in the MSc of the same name by me.
The *.sage* files contain data obtained from the LMFDB, namely the L-polynomials of simple abelian varieties for a given prime power q and a given dimension g.

The checks for dimension 3, 4 and 5 were done up to q = 8, 4 and 3 respectively due to hardware constraints.

For dimensions 3, 4 and 5, use the file *Weil polynomials of degree 6, 8* and *10* respectively.

I worked with the *.ipynb*-files, but added a copy of each as *.py*-file.
The first part of these files are about generating all q-Weil polynomials without real roots of the given degree using Theorems 9.2.3, 10.2.3 or 11.2.3, and comparing the result to the built-in function in Sage. The second part determines which of the generated q-Weil polynomials are characteristic polynomials of simple abelian varieties of the given dimension and the corresponding p-rank, as in Theorems 9.6 and 10.4, 11.4. These results are compared to the data from the LMFDB.
