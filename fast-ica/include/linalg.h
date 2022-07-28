/*
Copyright 2022 Mattia Orlandi

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#ifndef FAST_ICA_LINALG_H
#define FAST_ICA_LINALG_H

#include "matrix.h"

void to_hessenberg(const Matrix *h, const Matrix *m);
void qr_decomposition(const Matrix *q, const Matrix *r, const Matrix *m);
//void lin_solve(const Matrix *x, const Matrix *a, const Matrix *b);
void solve_eig(const Matrix *eig_vals, const Matrix *eig_vecs, const Matrix *m);
void covariance(const Matrix *cov_mat, const Matrix *mean_vec, const Matrix *m);

#endif //FAST_ICA_LINALG_H
