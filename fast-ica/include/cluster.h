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

#ifndef FAST_ICA_CLUSTER_H
#define FAST_ICA_CLUSTER_H

#include "../include/matrix.h"

/*
 * Basic self-operations (the trailing underscore means that the operation is performed in-place)
 */
fp norm_par(const Matrix *m);
fp mean_par(const Matrix *m);
void row_mean_par(const Matrix *r, const Matrix *m);
void col_mean_par(const Matrix *c, const Matrix *m);
fp std_par(const Matrix *m, fp mean_);
void transpose_par(const Matrix *t, const Matrix *m);

/*
 * Basic operations between two matrices (the trailing underscore means that the operation is performed in-place)
 */
void add_scalar_par(const Matrix *s, const Matrix *m, fp scalar);
void add_scalar_par_(const Matrix *m, fp scalar);
void scale_par(const Matrix *s, const Matrix *m, fp scalar);
void scale_par_(const Matrix *m, fp scalar);
void add_mat_par(const Matrix *s, const Matrix *m1, const Matrix *m2);
void add_mat_par_(const Matrix *m1, const Matrix *m2);
void add_row_par(const Matrix *s, const Matrix *m, const Matrix *r);
void add_row_par_(const Matrix *m, const Matrix *r);
void add_col_par(const Matrix *s, const Matrix *m, const Matrix *c);
void add_col_par_(const Matrix *m, const Matrix *c);
void sub_mat_par(const Matrix *s, const Matrix *m1, const Matrix *m2);
void sub_mat_par_(const Matrix *m1, const Matrix *m2);
void sub_row_par(const Matrix *s, const Matrix *m, const Matrix *r);
void sub_row_par_(const Matrix *m, const Matrix *r);
void sub_col_par(const Matrix *s, const Matrix *m, const Matrix *c);
void sub_col_par_(const Matrix *m, const Matrix *c);
void hadamard_par(const Matrix *s, const Matrix *m1, const Matrix *m2);
void hadamard_par_(const Matrix *m1, const Matrix *m2);
void hadamard_row_par(const Matrix *s, const Matrix *m, const Matrix *r);
void hadamard_row_par_(const Matrix *m, const Matrix *r);
void hadamard_col_par(const Matrix *s, const Matrix *m, const Matrix *c);
void hadamard_col_par_(const Matrix *m, const Matrix *c);
void div_mat_par(const Matrix *s, const Matrix *m1, const Matrix *m2);
void div_mat_par_(const Matrix *m1, const Matrix *m2);
void div_row_par(const Matrix *s, const Matrix *m, const Matrix *r);
void div_row_par_(const Matrix *m, const Matrix *r);
void div_col_par(const Matrix *s, const Matrix *m, const Matrix *c);
void div_col_par_(const Matrix *m, const Matrix *c);
void mat_mul_par(const Matrix *s, const Matrix *m1, const Matrix *m2);
void mat_mul_t1_par(const Matrix *s, const Matrix *m1, const Matrix *m2);
void mat_mul_t2_par(const Matrix *s, const Matrix *m1, const Matrix *m2);
fp dot_par(const Matrix *v1, const Matrix *v2);
void outer_par(const Matrix *s, const Matrix *v1, const Matrix *v2);

/*
 * Matrix manipulation
 */
void tri_up_par(const Matrix *m);

/*
 * Utils
 */
void copy_mat_par(const Matrix *s, const Matrix *m);

#endif //FAST_ICA_CLUSTER_H
