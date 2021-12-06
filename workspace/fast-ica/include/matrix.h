//
// Created by nihil on 20/09/21.
//

#ifndef FAST_ICA_MATRIX_H
#define FAST_ICA_MATRIX_H

#include <stdbool.h>
#include "pmsis.h"
#include "fp.h"

/*
 * Data type declaration and macros
 */
typedef struct Matrix {
    fp* data;
    const uint height;
    const uint width;
    bool on_cluster;
} Matrix;
#define MAT_DATA(m) (((Matrix *) (m))->data)
#define MAT_IDX(m, i, j) (((i) * ((m)->width)) + (j))
#define MAT_CELL(m, i, j) ((MAT_DATA(m))[MAT_IDX(m, i, j)])

/*
 * Initialize and free
 */
Matrix *new_mat(uint height, uint width, bool on_cluster);
Matrix *new_vec(uint length, bool on_cluster);  // column vector
void free_mat(Matrix *m);

/*
 * Transfer matrices between Fabric Controller and Cluster
 */
void mat_fc2cl(Matrix *m_fc, Matrix *m_cl);
void mat_fc2cl_async(Matrix *m_fc, Matrix *m_cl, pi_cl_dma_cmd_t *dma_cmd_data);
void mat_cl2fc(Matrix *m_cl, Matrix *m_fc);
void mat_cl2fc_async(Matrix *m_cl, Matrix *m_fc, pi_cl_dma_cmd_t *dma_cmd_data);

/*
 * Predefined matrices
 */
Matrix *eye(uint n, bool on_cluster);
Matrix *mat_randint(uint height, uint width, int min, int max, bool on_cluster);
Matrix *mat_rand(uint height, uint width, fp min, fp max, bool on_cluster);
Matrix *mat_randn(uint height, uint width, bool on_cluster);
Matrix *linspace(fp start, fp stop, int n_samples, bool on_cluster);
Matrix *from_array(const fp *data, uint height, uint width, bool on_cluster);

/*
 * Basic self-operations (the trailing underscore means that the operation is performed in-place)
 */
fp norm(Matrix *m);
fp mean(Matrix *m);
Matrix *row_mean(Matrix *m);
Matrix *col_mean(Matrix *m);
fp std(Matrix *m);
Matrix *row_std(Matrix *m);
Matrix *col_std(Matrix *m);
Matrix *transpose(Matrix *m);
Matrix *scale(Matrix *m, fp scalar);
void scale_(Matrix *m, fp scalar);

/*
 * Basic operations between two matrices (the trailing underscore means that the operation is performed in-place)
 */
Matrix *add_mat(Matrix *m1, Matrix *m2);
void add_mat_(Matrix *m1, Matrix *m2);
Matrix *add_row(Matrix *m, Matrix *r);
void add_row_(Matrix *m, Matrix *r);
Matrix *add_col(Matrix *m, Matrix *c);
void add_col_(Matrix *m, Matrix *c);
Matrix *sub_mat(Matrix *m1, Matrix *m2);
void sub_mat_(Matrix *m1, Matrix *m2);
Matrix *sub_row(Matrix *m, Matrix *r);
void sub_row_(Matrix *m, Matrix *r);
Matrix *sub_col(Matrix *m, Matrix *c);
void sub_col_(Matrix *m, Matrix *c);
Matrix *add_scalar(Matrix *m, fp scalar);
void add_scalar_(Matrix *m, fp scalar);
Matrix *hadamard(Matrix *m1, Matrix *m2);
void hadamard_(Matrix *m1, Matrix *m2);
Matrix *hadamard_row(Matrix *m, Matrix *r);
void hadamard_row_(Matrix *m, Matrix *r);
Matrix *hadamard_col(Matrix *m, Matrix *c);
void hadamard_col_(Matrix *m, Matrix *c);
Matrix *div_mat(Matrix *m1, Matrix *m2);
void div_mat_(Matrix *m1, Matrix *m2);
Matrix *div_row(Matrix *m, Matrix *r);
void div_row_(Matrix *m, Matrix *r);
Matrix *div_col(Matrix *m, Matrix *c);
void div_col_(Matrix *m, Matrix *c);
Matrix *mat_mul(Matrix *m1, Matrix *m2);
Matrix *mat_mul_t1(Matrix* m1, Matrix* m2);
Matrix *mat_mul_t2(Matrix* m1, Matrix* m2);
fp dot(Matrix *v1, Matrix *v2);
Matrix *outer(Matrix *v1, Matrix *v2);

/*
 * Boolean operations
 */
bool is_vector(Matrix *m);
bool is_row_vector(Matrix *m);
bool is_col_vector(Matrix *m);
bool are_equal(Matrix *m1, Matrix *m2, fp tol);
bool is_square(Matrix *m);

/*
 * Matrix manipulation
 */
Matrix *diagonal(Matrix *m);
void tri_up(Matrix *m);
Matrix *read_slice(Matrix *m, uint row_start, uint row_stop, uint col_start, uint col_stop);
void write_slice(Matrix *m1, Matrix *m2, uint row_start, uint col_start);
Matrix *extract_row(Matrix *m, uint k);
Matrix *extract_col(Matrix *m, uint k);
void paste_row(Matrix *m, Matrix *r, uint k);
void paste_col(Matrix *m, Matrix *c, uint k);

/*
 * Utils
 */
void print_mat(Matrix *m);
// void write_mat(const char *path, Matrix *m);
Matrix *copy_mat(Matrix *m);

#endif //FAST_ICA_MATRIX_H
