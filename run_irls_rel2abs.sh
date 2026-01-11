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
echo "Equal weights"
cat << EOF > irls_rel2abs.in
  1    2     1.01    1.0
  2    3     1.11    1.0
  1    3     2.01    1.0
  1    4     3.01    1.0
  2    4     1.97    1.0
  3    4     0.89    1.0
  1    5    10.00    1.0
  2    5     9.00    1.0
  3    5     8.00    1.0
  4    5     7.00    1.0
EOF
./irls_rel2abs 5 1 < irls_rel2abs.in
