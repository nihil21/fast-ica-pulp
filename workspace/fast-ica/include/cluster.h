//
// Created by nihil on 16/11/21.
//

#ifndef FAST_ICA_CLUSTER_H
#define FAST_ICA_CLUSTER_H

#include "../include/matrix.h"

typedef struct ClOpMMM {
    Matrix *op1;
    Matrix *op2;
    Matrix *res;
} ClOpMMM;

typedef struct ClOpMSM {
    Matrix *op1;
    fp *op2;
    Matrix *res;
} ClOpMSM;

typedef struct ClOpMSS {
    Matrix *op1;
    fp *op2;
    fp *res;
} ClOpMSS;

/*
 * Basic self-operations (the trailing underscore means that the operation is performed in-place)
 */
void norm_par(void *arg);
void mean_par(void *arg);
void row_mean_par(void *arg);
void col_mean_par(void *arg);
void std_par(void *arg);
void row_std_par(void *arg);
void col_std_par(void *arg);
void transpose_par(void *arg);
void scale_par(void *arg);
void scale_par_(void *arg);

/*
 * Basic operations between two matrices (the trailing underscore means that the operation is performed in-place)
 */
void add_mat_par(void *arg);
void add_mat_par_(void *arg);
void add_row_par(void *arg);
void add_row_par_(void *arg);
void add_col_par(void *arg);
void add_col_par_(void *arg);
void sub_mat_par(void *arg);
void sub_mat_par_(void *arg);
void sub_row_par(void *arg);
void sub_row_par_(void *arg);
void sub_col_par(void *arg);
void sub_col_par_(void *arg);
void add_scalar_par(void *arg);
void add_scalar_par_(void *arg);
void hadamard_par(void *arg);
void hadamard_par_(void *arg);
void hadamard_row_par(void *arg);
void hadamard_row_par_(void *arg);
void hadamard_col_par(void *arg);
void hadamard_col_par_(void *arg);
void div_mat_par(void *arg);
void div_mat_par_(void *arg);
void div_row_par(void *arg);
void div_row_par_(void *arg);
void div_col_par(void *arg);
void div_col_par_(void *arg);
void mat_mul_par(void *arg);
void mat_mul_t1_par(void *arg);
void mat_mul_t2_par(void *arg);

#endif //FAST_ICA_CLUSTER_H