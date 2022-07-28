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

#include "../include/cluster.h"
#include <stdio.h>

PI_L1 static fp buffer[CORES];  // buffer for reduction
PI_L1 static fp g_acc;  // global accumulator

/*
 * Get the norm of a given matrix
 */
fp norm_par(const Matrix *m) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    fp acc = 0;
    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            acc += POW(MAT_CELL(m, i, j), 2);
        }
    }
    buffer[core_id] = acc;

    pi_cl_team_barrier();

    g_acc = 0;
    if (core_id == 0) {
        for (int k = 0; k < CORES; k++)
            g_acc += buffer[k];
    }

    pi_cl_team_barrier();
    
    return g_acc;
}

/*
 * Get the mean of a given matrix
 */
fp mean_par(const Matrix *m) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    fp acc = 0;
    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            acc += MAT_CELL(m, i, j);
        }
    }
    buffer[core_id] = acc;

    pi_cl_team_barrier();
    
    g_acc = 0;
    if (core_id == 0) {
        for (int k = 0; k < CORES; k++)
            g_acc += buffer[k];
    }

    pi_cl_team_barrier();

    return g_acc;
}

/*
 * Get the mean of a given matrix along rows
 */
void row_mean_par(const Matrix *r, const Matrix *m) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    // Memory access not contiguous
    for (int j = j_start; j < j_end; j++) {
        fp acc = 0;
        for (int i = 0; i < m->height; i++) {
            acc += MAT_CELL(m, i, j);
        }
        MAT_CELL(r, 0, j) = acc / (fp) m->height;
    }

    pi_cl_team_barrier();
}

/*
 * Get the mean of a given matrix along columns
 */
void col_mean_par(const Matrix *c, const Matrix *m) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        fp acc = 0;
        for (int j = j_start; j < j_end; j++) {
            acc += MAT_CELL(m, i, j);
        }
        buffer[core_id] = acc;

        pi_cl_team_barrier();
        
        if (core_id == 0) {
            acc = 0;
            for (int k = 0; k < CORES; k++)
                acc += buffer[k];
            MAT_CELL(c, i, 0) = acc / (fp) m->width;
        }
    }

    pi_cl_team_barrier();
}

/*
 * Get the standard deviation of a given matrix
 */
fp std_par(const Matrix *m, fp mean_) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    fp acc = 0;
    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            acc += POW(MAT_CELL(m, i, j) - mean_, 2);
        }
    }
    buffer[core_id] = acc;

    pi_cl_team_barrier();
    
    g_acc = 0;
    if (core_id == 0) {
        for (int k = 0; k < CORES; k++)
            g_acc += buffer[k];
    }

    pi_cl_team_barrier();

    return g_acc;
}

/*
 * Transpose a matrix
 */
void transpose_par(const Matrix *t, const Matrix *m) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(t, j, i) = MAT_CELL(m, i, j);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Add scalar to matrix
 */
void add_scalar_par(const Matrix *s, const Matrix *m, const fp scalar) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) + scalar;
        }
    }

    pi_cl_team_barrier();
}

/*
 * Add scalar to matrix (in-place operation, it modifies the matrix)
 */
void add_scalar_par_(const Matrix *m, const fp scalar) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) += scalar;
        }
    }

    pi_cl_team_barrier();
}

/*
 * Scales a matrix by a scalar constant
 */
void scale_par(const Matrix *s, const Matrix *m, const fp scalar) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) * scalar;
        }
    }

    pi_cl_team_barrier();
}

/*
 * Scales a matrix by a scalar constant (in-place operation, it modifies the matrix)
 */
void scale_par_(const Matrix *m, const fp scalar) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) *= scalar;
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform element-wise addition between two matrices with same shape
 */
void add_mat_par(const Matrix *s, const Matrix *m1, const Matrix *m2) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m1->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) + MAT_CELL(m2, i, j);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform element-wise addition between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void add_mat_par_(const Matrix *m1, const Matrix *m2) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m1->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m1, i, j) += MAT_CELL(m2, i, j);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform row-wise addition between a matrix and a row vector with same width
 */
void add_row_par(const Matrix *s, const Matrix *m, const Matrix *r) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) + MAT_CELL(r, 0, i);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform row-wise addition between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void add_row_par_(const Matrix *m, const Matrix *r) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) += MAT_CELL(r, 0, i);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform column-wise addition between a matrix and a column vector with same height
 */
void add_col_par(const Matrix *s, const Matrix *m, const Matrix *c) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) + MAT_CELL(c, i, 0);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform column-wise addition between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void add_col_par_(const Matrix *m, const Matrix *c) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) += MAT_CELL(c, i, 0);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform element-wise subtraction between two matrices with same shape
 */
void sub_mat_par(const Matrix *s, const Matrix *m1, const Matrix *m2) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m1->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) - MAT_CELL(m2, i, j);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform element-wise subtraction between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void sub_mat_par_(const Matrix *m1, const Matrix *m2) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m1->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m1, i, j) -= MAT_CELL(m2, i, j);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform row-wise subtraction between a matrix and a row vector with same width
 */
void sub_row_par(const Matrix *s, const Matrix *m, const Matrix *r) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) - MAT_CELL(r, 0, i);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform row-wise subtraction between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void sub_row_par_(const Matrix *m, const Matrix *r) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) -= MAT_CELL(r, 0, i);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform column-wise subtraction between a matrix and a column vector with same height
 */
void sub_col_par(const Matrix *s, const Matrix *m, const Matrix *c) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) - MAT_CELL(c, i, 0);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform column-wise subtraction between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void sub_col_par_(const Matrix *m, const Matrix *c) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) -= MAT_CELL(c, i, 0);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform element-wise product (i.e. Hadamard) between two matrices with same shape
 */
void hadamard_par(const Matrix *s, const Matrix *m1, const Matrix *m2) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m1->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) * MAT_CELL(m2, i, j);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform element-wise product (i.e. Hadamard) between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void hadamard_par_(const Matrix *m1, const Matrix *m2) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m1->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m1, i, j) *= MAT_CELL(m2, i, j);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform row-wise product (i.e. Hadamard) between a matrix and a row vector with same width
 */
void hadamard_row_par(const Matrix *s, const Matrix *m, const Matrix *r) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) * MAT_CELL(r, 0, i);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform row-wise product (i.e. Hadamard) between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void hadamard_row_par_(const Matrix *m, const Matrix *r) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) *= MAT_CELL(r, 0, i);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform column-wise product (i.e. Hadamard) between a matrix and a column vector with same height
 */
void hadamard_col_par(const Matrix *s, const Matrix *m, const Matrix *c) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) * MAT_CELL(c, i, 0);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform column-wise product (i.e. Hadamard) between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void hadamard_col_par_(const Matrix *m, const Matrix *c) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) *= MAT_CELL(c, i, 0);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform element-wise division between two matrices with same shape
 */
void div_mat_par(const Matrix *s, const Matrix *m1, const Matrix *m2) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m1->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) / MAT_CELL(m2, i, j);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform element-wise division between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void div_mat_par_(const Matrix *m1, const Matrix *m2) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m1->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m1->width ? j_start + j_chunk : m1->width;

    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m1, i, j) /= MAT_CELL(m2, i, j);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform row-wise division between a matrix and a row vector with same width
 */
void div_row_par(const Matrix *s, const Matrix *m, const Matrix *r) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) / MAT_CELL(r, 0, i);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform row-wise division between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void div_row_par_(const Matrix *m, const Matrix *r) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) /= MAT_CELL(r, 0, i);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform column-wise division between a matrix and a column vector with same height
 */
void div_col_par(const Matrix *s, const Matrix *m, const Matrix *c) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j) / MAT_CELL(c, i, 0);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform column-wise division between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void div_col_par_(const Matrix *m, const Matrix *c) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;

    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(m, i, j) /= MAT_CELL(c, i, 0);
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform matrix multiplication between an AxB matrix and a BxC matrix
 */
void mat_mul_par(const Matrix *s, const Matrix *m1, const Matrix *m2) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m2->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m2->width ? j_start + j_chunk : m2->width;

    // i-k-j loop to optimize memory access
    for (int i = 0; i < m1->height; i++) {
        for (int k = 0; k < m1->width; k++) {
            for (int j = j_start; j < j_end; j++) {
                MAT_CELL(s, i, j) += MAT_CELL(m1, i, k) * MAT_CELL(m2, k, j);
            }
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform matrix multiplication between a BxA matrix and a BxC matrix, without transposing the first matrix
 */
void mat_mul_t1_par(const Matrix *s, const Matrix *m1, const Matrix *m2) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m2->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m2->width ? j_start + j_chunk : m2->width;

    // k-i-j loop to optimize memory access
    for (int k = 0; k < m1->height; k++) {
        for (int i = 0; i < m1->width; i++) {
            for (int j = j_start; j < j_end; j++) {
                MAT_CELL(s, i, j) += MAT_CELL(m1, k, i) * MAT_CELL(m2, k, j);
            }
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform matrix multiplication between an AxB matrix and a CxB matrix, without transposing the second matrix
 */
void mat_mul_t2_par(const Matrix *s, const Matrix *m1, const Matrix *m2) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m2->height + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m2->height ? j_start + j_chunk : m2->height;
    
    // i-j-k loop to optimize memory access
    for (int i = 0; i < m1->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            for (int k = 0; k < m1->width; k++) {
                MAT_CELL(s, i, j) += MAT_CELL(m1, i, k) * MAT_CELL(m2, j, k);
            }
        }
    }

    pi_cl_team_barrier();
}

/*
 * Perform dot product between two vectors
 */
fp dot_par(const Matrix *v1, const Matrix *v2) {
    const uint len = v1->height == 1 ? v1->width : v1->height;

    const uint core_id = pi_core_id();
    const uint i_chunk = (len + CORES - 1) / CORES;
    const uint i_start = core_id * i_chunk;
    const uint i_end = i_start + i_chunk < len ? i_start + i_chunk : len;

    fp acc = 0;
    for (int i = i_start; i < i_end; i++) {
        acc += MAT_DATA(v1)[i] * MAT_DATA(v2)[i];
    }
    buffer[core_id] = acc;
    
    pi_cl_team_barrier();

    g_acc = 0;
    if (core_id == 0) {
        for (int k = 0; k < CORES; k++)
            g_acc += buffer[k];
    }

    pi_cl_team_barrier();

    return g_acc;
}

/*
 * Perform outer product between two vectors
 */
void outer_par(const Matrix *s, const Matrix *v1, const Matrix *v2) {
    const uint len1 = v1->height == 1 ? v1->width : v1->height;
    const uint len2 = v2->height == 1 ? v2->width : v2->height;

    const uint core_id = pi_core_id();
    const uint j_chunk = (len2 + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < len2 ? j_start + j_chunk : len2;

    for (int i = 0; i < len1; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_DATA(v1)[i] * MAT_DATA(v2)[j];
        }
    }

    pi_cl_team_barrier();
}

/*
 * Set cells below the main diagonal to zero (in-place operation, it modifies the matrix)
 */
void tri_up_par(const Matrix *m) {
    const uint core_id = pi_core_id();
    const uint i_chunk = (m->height + CORES - 1) / CORES;
    const uint i_start = core_id * i_chunk;
    const uint i_end = i_start + i_chunk < m->height ? i_start + i_chunk : m->height;

    for (int i = i_start; i < i_end; i++) {
        for (int j = 0; j < i; j++) {
            MAT_CELL(m, i, j) = 0;
        }
    }

    pi_cl_team_barrier();
}

/*
 * Copy a matrix into a new matrix
 */
void copy_mat_par(const Matrix *s, const Matrix *m) {
    const uint core_id = pi_core_id();
    const uint j_chunk = (m->width + CORES - 1) / CORES;
    const uint j_start = core_id * j_chunk;
    const uint j_end = j_start + j_chunk < m->width ? j_start + j_chunk : m->width;
    
    for (int i = 0; i < m->height; i++) {
        for (int j = j_start; j < j_end; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j);
        }
    }

    pi_cl_team_barrier();
}
