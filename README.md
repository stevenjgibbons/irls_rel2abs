# irls_rel2abs  
Fortran program to convert relative distances (in time or space) to absolute distances using the iteratively reweighted least squares and the linear algebra formulation of VanDecar and Crosson, BSSA, v80, pp. 150-169, 1990 https://doi.org/10.1785/BSSA0800010150.  

Compile program by typing
```
make irls_rel2abs
```
(you will need to have valid paths to LAPACK and BLAS libraries).  

The program *irls_rel2abs* takes in two integer arguments: NVALS, the number of absolute values, and IFIXVL, the number of that value that should be fixed to zero.  
It then reads from standard input lines of the form
```
I   J   rel_val_J_minus_I    weight
```

A complete example is given in the file *run_irls_rel2abs.sh*  
```
#!/bin/sh
echo "Unequal weights"
cat << EOF > ./irls_rel2abs.in
  1    2     1.01    1.0
  2    3     1.11    0.8
  1    3     2.01    0.99
  1    4     3.51    1.00
  2    4     1.97    0.8
  3    4     0.89    0.5
  1    5    10.00    0.5
  2    5     9.00    0.5
  3    5     8.00    0.5
  4    5     7.00    0.5
EOF
./irls_rel2abs 5 1 < irls_rel2abs.in
```


