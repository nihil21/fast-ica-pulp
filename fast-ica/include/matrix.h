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

#ifndef FAST_ICA_MATRIX_H
#define FAST_ICA_MATRIX_H

#include <stdbool.h>
#include "pmsis.h"
#include "fp.h"

/*
 * Data type declaration and macros
 */
typedef struct {
    fp *data;
    uint height;
    uint width;
    uint offset;
    bool on_cluster;
} Matrix;

#define MAT_DATA(m) (((Matrix *) (m))->data)
#define MAT_IDX(m, i, j) (((i) * ((m)->offset)) + (j))
#define MAT_CELL(m, i, j) ((MAT_DATA(m))[MAT_IDX(m, i, j)])

/*
 * Initialize and free
 */
Matrix *new_mat(uint height, uint width, bool on_cluster);
void free_mat(Matrix *m);

/*
 * Transfer matrices between Fabric Controller and Cluster
 */
void mat_fc2cl(const Matrix *m_fc, const Matrix *m_cl);
void mat_fc2cl_async(const Matrix *m_fc, const Matrix *m_cl, pi_cl_dma_cmd_t *dma_cmd_data);
void mat_cl2fc(const Matrix *m_cl, const Matrix *m_fc);
void mat_cl2fc_async(const Matrix *m_cl, const Matrix *m_fc, pi_cl_dma_cmd_t *dma_cmd_data);

/*
 * Predefined matrices
 */
void eye(const Matrix *e);
void mat_randint(const Matrix *m, int min, int max);
void mat_rand(const Matrix *m, fp min, fp max);
void mat_randn(const Matrix *m);
void linspace(const Matrix *ls, fp start, fp stop, int n_samples);

/*
 * Basic self-operations (the trailing underscore means that the operation is performed in-place)
 */
fp norm(const Matrix *m);
fp mean(const Matrix *m);
void row_mean(const Matrix *r, const Matrix *m);
void col_mean(const Matrix *c, const Matrix *m);
fp std(const Matrix *m);
void row_std(const Matrix *r, const Matrix *m);
void col_std(const Matrix *c, const Matrix *m);
void transpose(const Matrix *t, const Matrix *m);

/*
 * Basic operations between two matrices (the trailing underscore means that the operation is performed in-place)
 */
void add_scalar(const Matrix *s, const Matrix *m, fp scalar);
void add_scalar_(const Matrix *m, fp scalar);
void scale(const Matrix *s, const Matrix *m, fp scalar);
void scale_(const Matrix *m, fp scalar);
void add_mat(const Matrix *s, const Matrix *m1, const Matrix *m2);
void add_mat_(const Matrix *m1, const Matrix *m2);
void add_row(const Matrix *s, const Matrix *m, const Matrix *r);
void add_row_(const Matrix *m, const Matrix *r);
void add_col(const Matrix *s, const Matrix *m, const Matrix *c);
void add_col_(const Matrix *m, const Matrix *c);
void sub_mat(const Matrix *s, const Matrix *m1, const Matrix *m2);
void sub_mat_(const Matrix *m1, const Matrix *m2);
void sub_row(const Matrix *s, const Matrix *m, const Matrix *r);
void sub_row_(const Matrix *m, const Matrix *r);
void sub_col(const Matrix *s, const Matrix *m, const Matrix *c);
void sub_col_(const Matrix *m, const Matrix *c);
void hadamard(const Matrix *s, const Matrix *m1, const Matrix *m2);
void hadamard_(const Matrix *m1, const Matrix *m2);
void hadamard_row(const Matrix *s, const Matrix *m, const Matrix *r);
void hadamard_row_(const Matrix *m, const Matrix *r);
void hadamard_col(const Matrix *s, const Matrix *m, const Matrix *c);
void hadamard_col_(const Matrix *m, const Matrix *c);
void div_mat(const Matrix *s, const Matrix *m1, const Matrix *m2);
void div_mat_(const Matrix *m1, const Matrix *m2);
void div_row(const Matrix *s, const Matrix *m, const Matrix *r);
void div_row_(const Matrix *m, const Matrix *r);
void div_col(const Matrix *s, const Matrix *m, const Matrix *c);
void div_col_(const Matrix *m, const Matrix *c);
void mat_mul(const Matrix *s, const Matrix *m1, const Matrix *m2);
void mat_mul_t1(const Matrix *s, const Matrix *m1, const Matrix *m2);
void mat_mul_t2(const Matrix *s, const Matrix *m1, const Matrix *m2);
fp dot(const Matrix *v1, const Matrix *v2);
void outer(const Matrix *s, const Matrix *v1, const Matrix *v2);

/*
 * Boolean operations
 */
bool is_vector(const Matrix *m);
bool is_row_vector(const Matrix *m);
bool is_col_vector(const Matrix *m);
bool are_equal(const Matrix *m1, const Matrix *m2, fp tol);
bool is_square(const Matrix *m);

/*
 * Matrix manipulation
 */
Matrix diagonal(const Matrix *m);
void tri_up(const Matrix *m);
Matrix slice(const Matrix *m, uint row_start, uint row_stop, uint col_start, uint col_stop);
Matrix extract_row(const Matrix *m, uint k);
Matrix extract_col(const Matrix *m, uint k);
void compact(Matrix *m, uint row_start, uint col_start, uint orig_len);

/*
 * Utils
 */
void print_mat(const Matrix *m);
// void write_mat(const char *path, const Matrix *m);
void copy_mat(const Matrix *s, const Matrix *m);

#endif //FAST_ICA_MATRIX_H
