//
// Created by nihil on 16/11/21.
//

#include "../include/cluster.h"

PI_L1 static fp buffer[CORES];

/*
 * Get the norm of a given matrix
 */
void norm_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMSS *) arg)->op1;
    fp *s = ((ClOpMSS *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    buffer[core_id] = 0;
    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            buffer[core_id] += POW(MAT_CELL(m, i, j), 2);
        }
    }
    pi_cl_team_barrier();
    
    if (core_id == 0) {
        for (int k = 0; k < CORES; k++)
            *s += buffer[k];
    }
}

/*
 * Get the mean of a given matrix
 */
void mean_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMSS *) arg)->op1;
    fp *s = ((ClOpMSS *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    buffer[core_id] = 0;
    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            buffer[core_id] += MAT_CELL(m, i, j);
        }
    }
    pi_cl_team_barrier();
    
    if (core_id == 0) {
        for (int k = 0; k < CORES; k++)
            *s += buffer[k];
    }
}

/*
 * Get the mean of a given matrix along rows
 */
void row_mean_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *r = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    // Memory access not contiguous
    for (int j = j_start; j < j_end; j++) {
        fp acc = 0;
        for (int i = 0; i < m->height; i++) {
            acc += MAT_CELL(m, i, j);
        }
        MAT_CELL(r, 0, j) = acc / (fp) m->height;
    }
}

/*
 * Get the mean of a given matrix along columns
 */
void col_mean_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *c = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        buffer[core_id] = 0;
        for (int j = j_start; j < j_end; j++) {
            buffer[core_id] += MAT_CELL(m, i, j);
        }
        pi_cl_team_barrier();
        
        if (core_id == 0) {
            fp acc = 0;
            for (int k = 0; k < CORES; k++)
                acc += buffer[k];
            MAT_CELL(c, i, 0) = acc / (fp) m->width;
        }
    }
}

/*
 * Get the standard deviation of a given matrix
 */
void std_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMSS *) arg)->op1;
    fp *mean_ = ((ClOpMSS *) arg)->op2;
    fp *s = ((ClOpMSS *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    buffer[core_id] = 0;
    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            buffer[core_id] += POW(MAT_CELL(m, i, j) - *mean_, 2);
        }
    }
    pi_cl_team_barrier();
    
    if (core_id == 0) {
        for (int k = 0; k < CORES; k++)
            *s += buffer[k];
    }
}

/*
 * Get the standard deviation of a given matrix along rows
 */
void row_std_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *mean_ = ((ClOpMMM *) arg)->op2;
    Matrix *r = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    // Memory access not contiguous
    for (int j = j_start; j < j_end; j++) {
        fp acc = 0;
        for (int i = 0; i < m->height; i++) {
            acc += POW(MAT_CELL(m, i, j) - MAT_CELL(mean_, 0, j), 2);
        }
        MAT_CELL(r, 0, j) = SQRT(acc / (fp) m->height);
    }
}

/*
 * Get the standard deviation of a given matrix along columns
 */
void col_std_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *mean_ = ((ClOpMMM *) arg)->op2;
    Matrix *c = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        buffer[core_id] = 0;
        for (int j = j_start; j < j_end; j++) {
            buffer[core_id] += POW(MAT_CELL(m, i, j) - MAT_CELL(mean_, i, 0), 2);
        }
        pi_cl_team_barrier();
        
        if (core_id == 0) {
            fp acc = 0;
            for (int k = 0; k < CORES; k++)
                acc += buffer[k];
            MAT_CELL(c, i, 0) = SQRT(acc / (fp) m->width);
        }
    }
}

/*
 * Transpose a matrix
 */
void transpose_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *t = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(t, j, i) = MAT_CELL(m, i, j);
        }
    }
}

/*
 * Scales a matrix by a scalar constant
 */
void scale_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMSM *) arg)->op1;
    fp *scalar = ((ClOpMSM *) arg)->op2;
    Matrix *s = ((ClOpMSM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) * (*scalar);
        }
    }
}

/*
 * Scales a matrix by a scalar constant (in-place operation, it modifies the matrix)
 */
void scale_par_(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMSM *) arg)->op1;
    fp *scalar = ((ClOpMSM *) arg)->op2;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) *= (*scalar);
        }
    }
}

/*
 * Perform element-wise addition between two matrices with same shape
 */
void add_mat_par(void *arg) {
    // Unpack arguments
    Matrix *m1 = ((ClOpMMM *) arg)->op1;
    Matrix *m2 = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m1->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) + MAT_CELL(m2, i, j);
        }
    }
}

/*
 * Perform element-wise addition between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void add_mat_par_(void *arg) {
    // Unpack arguments
    Matrix *m1 = ((ClOpMMM *) arg)->op1;
    Matrix *m2 = ((ClOpMMM *) arg)->op2;

    int core_id = pi_core_id();
    uint j_chunk = (m1->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m1, i, j) += MAT_CELL(m2, i, j);
        }
    }
}

/*
 * Perform row-wise addition between a matrix and a row vector with same width
 */
void add_row_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *r = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) + MAT_CELL(r, 0, i);
        }
    }
}

/*
 * Perform row-wise addition between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void add_row_par_(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *r = ((ClOpMMM *) arg)->op2;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) += MAT_CELL(r, 0, i);
        }
    }
}

/*
 * Perform column-wise addition between a matrix and a column vector with same height
 */
void add_col_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *c = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) + MAT_CELL(c, i, 0);
        }
    }
}

/*
 * Perform column-wise addition between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void add_col_par_(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *c = ((ClOpMMM *) arg)->op2;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) += MAT_CELL(c, i, 0);
        }
    }
}

/*
 * Perform element-wise subtraction between two matrices with same shape
 */
void sub_mat_par(void *arg) {
    // Unpack arguments
    Matrix *m1 = ((ClOpMMM *) arg)->op1;
    Matrix *m2 = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m1->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) - MAT_CELL(m2, i, j);
        }
    }
}

/*
 * Perform element-wise subtraction between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void sub_mat_par_(void *arg) {
    // Unpack arguments
    Matrix *m1 = ((ClOpMMM *) arg)->op1;
    Matrix *m2 = ((ClOpMMM *) arg)->op2;

    int core_id = pi_core_id();
    uint j_chunk = (m1->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m1, i, j) -= MAT_CELL(m2, i, j);
        }
    }
}

/*
 * Perform row-wise subtraction between a matrix and a row vector with same width
 */
void sub_row_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *r = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) - MAT_CELL(r, 0, i);
        }
    }
}

/*
 * Perform row-wise subtraction between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void sub_row_par_(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *r = ((ClOpMMM *) arg)->op2;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) -= MAT_CELL(r, 0, i);
        }
    }
}

/*
 * Perform column-wise subtraction between a matrix and a column vector with same height
 */
void sub_col_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *c = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) - MAT_CELL(c, i, 0);
        }
    }
}

/*
 * Perform column-wise subtraction between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void sub_col_par_(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *c = ((ClOpMMM *) arg)->op2;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) -= MAT_CELL(c, i, 0);
        }
    }
}

/*
 * Add scalar to matrix
 */
void add_scalar_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMSM *) arg)->op1;
    fp *scalar = ((ClOpMSM *) arg)->op2;
    Matrix *s = ((ClOpMSM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) + (*scalar);
        }
    }
}

/*
 * Add scalar to matrix (in-place operation, it modifies the matrix)
 */
void add_scalar_par_(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMSM *) arg)->op1;
    fp *scalar = ((ClOpMSM *) arg)->op2;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) += (*scalar);
        }
    }
}

/*
 * Perform element-wise product (i.e. Hadamard) between two matrices with same shape
 */
void hadamard_par(void *arg) {
    // Unpack arguments
    Matrix *m1 = ((ClOpMMM *) arg)->op1;
    Matrix *m2 = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m1->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) * MAT_CELL(m2, i, j);
        }
    }
}

/*
 * Perform element-wise product (i.e. Hadamard) between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void hadamard_par_(void *arg) {
    // Unpack arguments
    Matrix *m1 = ((ClOpMMM *) arg)->op1;
    Matrix *m2 = ((ClOpMMM *) arg)->op2;

    int core_id = pi_core_id();
    uint j_chunk = (m1->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m1, i, j) *= MAT_CELL(m2, i, j);
        }
    }
}

/*
 * Perform row-wise product (i.e. Hadamard) between a matrix and a row vector with same width
 */
void hadamard_row_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *r = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) * MAT_CELL(r, 0, i);
        }
    }
}

/*
 * Perform row-wise product (i.e. Hadamard) between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void hadamard_row_par_(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *r = ((ClOpMMM *) arg)->op2;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) *= MAT_CELL(r, 0, i);
        }
    }
}

/*
 * Perform column-wise product (i.e. Hadamard) between a matrix and a column vector with same height
 */
void hadamard_col_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *c = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) * MAT_CELL(c, i, 0);
        }
    }
}

/*
 * Perform column-wise product (i.e. Hadamard) between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void hadamard_col_par_(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *c = ((ClOpMMM *) arg)->op2;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) *= MAT_CELL(c, i, 0);
        }
    }
}

/*
 * Perform element-wise division between two matrices with same shape
 */
void div_mat_par(void *arg) {
    // Unpack arguments
    Matrix *m1 = ((ClOpMMM *) arg)->op1;
    Matrix *m2 = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m1->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) / MAT_CELL(m2, i, j);
        }
    }
}

/*
 * Perform element-wise division between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void div_mat_par_(void *arg) {
    // Unpack arguments
    Matrix *m1 = ((ClOpMMM *) arg)->op1;
    Matrix *m2 = ((ClOpMMM *) arg)->op2;

    int core_id = pi_core_id();
    uint j_chunk = (m1->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m1, i, j) /= MAT_CELL(m2, i, j);
        }
    }
}

/*
 * Perform row-wise division between a matrix and a row vector with same width
 */
void div_row_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *r = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) / MAT_CELL(r, 0, i);
        }
    }
}

/*
 * Perform row-wise division between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void div_row_par_(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *r = ((ClOpMMM *) arg)->op2;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) /= MAT_CELL(r, 0, i);
        }
    }
}

/*
 * Perform column-wise division between a matrix and a column vector with same height
 */
void div_col_par(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *c = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) / MAT_CELL(c, i, 0);
        }
    }
}

/*
 * Perform column-wise division between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void div_col_par_(void *arg) {
    // Unpack arguments
    Matrix *m = ((ClOpMMM *) arg)->op1;
    Matrix *c = ((ClOpMMM *) arg)->op2;

    int core_id = pi_core_id();
    uint j_chunk = (m->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) /= MAT_CELL(c, i, 0);
        }
    }
}

/*
 * Perform matrix multiplication between an AxB matrix and a BxC matrix
 */
void mat_mul_par(void *arg) {
    // Unpack arguments
    Matrix *m1 = ((ClOpMMM *) arg)->op1;
    Matrix *m2 = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m2->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m2->width ? j_start + j_chunk : m2->width;

    // i-k-j loop to optimize memory access
    for (uint i = 0; i < m1->height; i++) {
        for (uint k = 0; k < m1->width; k++) {
            for (uint j = j_start; j < j_end; j++) {
                MAT_CELL(s, i, j) += MAT_CELL(m1, i, k) * MAT_CELL(m2, k, j);
            }
        }
    }
}

/*
 * Perform matrix multiplication between a BxA matrix and a BxC matrix, without transposing the first matrix
 */
void mat_mul_t1_par(void *arg) {
    // Unpack arguments
    Matrix *m1 = ((ClOpMMM *) arg)->op1;
    Matrix *m2 = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m2->width + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m2->width ? j_start + j_chunk : m2->width;

    // k-i-j loop to optimize memory access
    for (uint k = 0; k < m1->height; k++) {
        for (uint i = 0; i < m1->width; i++) {
            for (uint j = j_start; j < j_end; j++) {
                MAT_CELL(s, i, j) += MAT_CELL(m1, k, i) * MAT_CELL(m2, k, j);
            }
        }
    }
}

/*
 * Perform matrix multiplication between an AxB matrix and a CxB matrix, without transposing the second matrix
 */
void mat_mul_t2_par(void *arg) {
    // Unpack arguments
    Matrix *m1 = ((ClOpMMM *) arg)->op1;
    Matrix *m2 = ((ClOpMMM *) arg)->op2;
    Matrix *s = ((ClOpMMM *) arg)->res;

    int core_id = pi_core_id();
    uint j_chunk = (m2->height + CORES - 1) / CORES;
    uint j_start = core_id * j_chunk;
    uint j_end = j_start + j_chunk < m2->height ? j_start + j_chunk : m2->height;

    // i-j-k loop to optimize memory access
    for (uint i = 0; i < m1->height; i++) {
        for (uint j = j_start; j < j_end; j++) {
            for (uint k = 0; k < m1->width; k++) {
                MAT_CELL(s, i, j) += MAT_CELL(m1, i, k) * MAT_CELL(m2, j, k);
            }
        }
    }
}