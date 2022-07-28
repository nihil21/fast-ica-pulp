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

#include <string.h>
#include "../include/matrix.h"
#include "../include/random.h"
#include "../include/cluster.h"
#include "../include/utils.h"

extern struct pi_device cluster_dev;

/*
 * Initialize a zero matrix given its dimensions
 */
Matrix *new_mat(const uint height, const uint width, const bool on_cluster) {
#if defined(DEBUG)
    assert(height > 0 && width > 0, "Matrix height and width should be greater than 0.");
#endif
    Matrix *m;
    if (on_cluster) {
        // Allocate memory for the struct (Cluster's L1 memory)
        m = pi_l1_malloc(&cluster_dev, sizeof(Matrix));
#if defined(DEBUG)
        assert(m != NULL, "Could not allocate Matrix struct (Cluster's L1 memory).");
#endif

        // Allocate memory for the array field (Cluster's L1 memory)
        m->data = pi_l1_malloc(&cluster_dev, height * width * sizeof(fp));
#if defined(DEBUG)
        assert(m->data != NULL, "Could not allocate Matrix data (Cluster's L1 memory).");
#endif
        memset(m->data, 0, height * width * sizeof(fp));
    } else {
        // Allocate memory for the struct (Fabric Controller's L2 memory)
        m = pi_l2_malloc(sizeof(Matrix));
#if defined(DEBUG)
        assert(m != NULL, "Could not allocate Matrix struct (Fabric Controller's L2 memory).");
#endif

        // Allocate memory for the array field (Fabric Controller's L2 memory)
        m->data = pi_l2_malloc(height * width * sizeof(fp));
#if defined(DEBUG)
        assert(m->data != NULL, "Could not allocate Matrix data (Fabric Controller's L2 memory).");
#endif
        memset(m->data, 0, height * width * sizeof(fp));  // initialize to zero
    }

    // Set other fields
    m->height = height;
    m->width = width;
    m->offset = width;
    m->on_cluster = on_cluster;

    return m;
}

/*
 * Free the memory allocated for a given matrix
 */
void free_mat(Matrix *m) {
    if (m != NULL) {
        if (m->data != NULL) {
            if (m->on_cluster)
                pi_l1_free(&cluster_dev, m->data, m->height * m->width * sizeof(fp));
            else
                pi_l2_free(m->data, m->height * m->width * sizeof(fp));
            m->data = NULL;
        }
        if (m->on_cluster)
            pi_l1_free(&cluster_dev, m, sizeof(Matrix));
        else
            pi_l2_free(m, sizeof(Matrix));
        m = NULL;
    }
}

/*
 * Move matrices from Fabric Controller to Cluster
 */
void mat_fc2cl(const Matrix *m_fc, const Matrix *m_cl) {
    // Transfer data field
    pi_cl_dma_cmd_t copy;
    pi_cl_dma_cmd((uint32_t) MAT_DATA(m_fc), (uint32_t) MAT_DATA(m_cl), m_fc->height * m_fc->width * sizeof(fp), PI_CL_DMA_DIR_EXT2LOC, &copy);
    pi_cl_dma_cmd_wait(&copy);
}

/*
 * Move matrices from Fabric Controller to Cluster (async)
 */
void mat_fc2cl_async(const Matrix *m_fc, const Matrix *m_cl, pi_cl_dma_cmd_t *dma_cmd_data) {
    // Transfer data field
    pi_cl_dma_cmd((uint32_t) MAT_DATA(m_fc), (uint32_t) MAT_DATA(m_cl), m_fc->height * m_fc->width * sizeof(fp), PI_CL_DMA_DIR_EXT2LOC, dma_cmd_data);
}

/*
 * Move matrices from Cluster to Fabric Controller
 */
void mat_cl2fc(const Matrix *m_cl, const Matrix *m_fc) {
    // Transfer data field
    pi_cl_dma_cmd_t copy;
    pi_cl_dma_cmd((uint32_t) MAT_DATA(m_fc), (uint32_t) MAT_DATA(m_cl), m_cl->height * m_cl->width * sizeof(fp), PI_CL_DMA_DIR_LOC2EXT, &copy);
    pi_cl_dma_cmd_wait(&copy);
}

/*
 * Move matrices from Cluster to Fabric Controller (async)
 */
void mat_cl2fc_async(const Matrix *m_cl, const Matrix *m_fc, pi_cl_dma_cmd_t *dma_cmd_data) {
    // Transfer data field
    pi_cl_dma_cmd((uint32_t) MAT_DATA(m_fc), (uint32_t) MAT_DATA(m_cl), m_cl->height * m_cl->width * sizeof(fp), PI_CL_DMA_DIR_LOC2EXT, dma_cmd_data);
}

/*
 * Fill the identity matrix
 */
void eye(const Matrix *e) {
#if defined(DEBUG)
    assert(is_square(e), "The matrix e must be square.");
#endif
    memset(e->data, 0, e->height * e->width * sizeof(fp));

    for (int k = 0; k < e->height; k++)
        MAT_CELL(e, k, k) = 1;
}

/*
 * Fill a matrix with random uniform integers in a given range
 */
void mat_randint(const Matrix *m, const int min, const int max) {
    for (int i = 0; i < m->height; i++) {
        for (int j = 0; j < m->width; j++) {
            MAT_CELL(m, i, j) = (fp) uniform_randint_range(min, max);
        }
    }
}

/*
 * Fill a matrix with random uniform numbers in a given range
 */
void mat_rand(const Matrix *m, const fp min, const fp max) {
    for (int i = 0; i < m->height; i++) {
        for (int j = 0; j < m->width; j++) {
            MAT_CELL(m, i, j) = uniform_rand_range(min, max);
        }
    }
}

/*
 * Fill a matrix with random normal numbers
 */
void mat_randn(const Matrix *m) {
    for (int i = 0; i < m->height; i++) {
        for (int j = 0; j < m->width; j++) {
            MAT_CELL(m, i, j) = standard_rand();
        }
    }
}

/*
 * Build a linear space with the given range and number of samples
 */
void linspace(const Matrix *ls, const fp start, const fp stop, const int n_samples) {
#if defined(DEBUG)
    assert(is_vector(ls), "The ls matrix should be a vector.");
    assert(start < stop, "The stop argument should be greater than the start argument.");
#endif

    fp step = (stop - start) / (fp) (n_samples - 1);
    if (ls->height == 1) {  // row vector
        for (int i = 0; i < n_samples; i++)
            MAT_CELL(ls, 0, i) = start + step * (fp) i;
    } else {  // column vector
        for (int i = 0; i < n_samples; i++)
            MAT_CELL(ls, i, 0) = start + step * (fp) i;
    }
}

/*
 * Compute the L2 norm of a given matrix
 */
fp norm(const Matrix *m) {
    fp acc = 0;
    if (m->on_cluster) {
        acc = norm_par(m);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                acc += MAT_CELL(m, i, j) * MAT_CELL(m, i, j);
            }
        }
    }
    
    return SQRT(acc);
}

/*
 * Get the mean of a given matrix
 */
fp mean(const Matrix *m) {
    fp acc = 0;
    if (m->on_cluster) {
        acc = mean_par(m);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                acc += MAT_CELL(m, i, j);
            }
        }
    }
    
    return acc / (fp) (m->height * m->width);
}

/*
 * Get the mean of a given matrix along rows
 */
void row_mean(const Matrix *r, const Matrix *m) {
#if defined(DEBUG)
    assert(r->on_cluster == m->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(is_row_vector(r) && r->width == m->width, "The r matrix should be a row vector with the same width as m.");
#endif

    if (r->on_cluster) {
        row_mean_par(r, m);  // compute it in parallel
    } else {
        // Memory access not contiguous
        for (int j = 0; j < m->width; j++) {
            fp acc = 0;
            for (int i = 0; i < m->height; i++) {
                acc += MAT_CELL(m, i, j);
            }
            MAT_CELL(r, 0, j) = acc / (fp) m->height;
        }
    }
}

/*
 * Get the mean of a given matrix along columns
 */
void col_mean(const Matrix *c, const Matrix *m) {
#if defined(DEBUG)
    assert(c->on_cluster == m->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(is_col_vector(c) && c->height == m->height, "The c matrix should be a column vector with the same height as m.");
#endif

    if (c->on_cluster) {
        col_mean_par(c, m);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            fp acc = 0;
            for (int j = 0; j < m->width; j++) {
                acc += MAT_CELL(m, i, j);
            }
            MAT_CELL(c, i, 0) = acc / (fp) m->width;
        }
    }
}

/*
 * Get the standard deviation of a given matrix
 */
fp std(const Matrix *m) {
    // Compute mean
    fp mean_ = mean(m);

    fp acc = 0;
    if (m->on_cluster) {
        acc = std_par(m, mean_);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                acc += POW(MAT_CELL(m, i, j) - mean_, 2);
            }
        }
    }
    
    return SQRT(acc / (fp) (m->height * m->width));
}

/*
 * Get the standard deviation of a given matrix along rows
 */
void row_std(const Matrix *r, const Matrix *m) {
#if defined(DEBUG)
    assert(!r->on_cluster && !m->on_cluster, "This operation is supported only for matrices in L2 memory.");
    assert(is_row_vector(r) && r->width == m->width, "The r matrix should be a row vector with the same width as m.");
#endif

    // Allocate matrix
    Matrix *mean_ = new_mat(1, m->width, false);
    // Compute mean
    row_mean(mean_, m);

    // Memory access not contiguous
    for (int j = 0; j < m->width; j++) {
        fp acc = 0;
        for (int i = 0; i < m->height; i++) {
            acc += POW(MAT_CELL(m, i, j) - MAT_CELL(mean_, 0, j), 2);
        }
        MAT_CELL(r, 0, j) = SQRT(acc / (fp) m->height);
    }

    // Free memory
    free_mat(mean_);
}

/*
 * Get the standard deviation of a given matrix along columns
 */
void col_std(const Matrix *c, const Matrix *m) {
#if defined(DEBUG)
    assert(!c->on_cluster && !m->on_cluster, "This operation is supported only for matrices in L2 memory.");
    assert(is_col_vector(c) && c->height == m->height, "The c matrix should be a column vector with the same height as m.");
#endif

    // Allocate matrix
    Matrix *mean_ = new_mat(m->height, 1, false);
    // Compute mean
    col_mean(mean_, m);

    for (int i = 0; i < m->height; i++) {
        fp acc = 0;
        for (int j = 0; j < m->width; j++) {
            acc += POW(MAT_CELL(m, i, j) - MAT_CELL(mean_, i, 0), 2);
        }
        MAT_CELL(c, i, 0) = SQRT(acc / (fp) m->width);
    }

    // Free memory
    free_mat(mean_);
}

/*
 * Transpose a matrix
 */
void transpose(const Matrix *t, const Matrix* m) {
#if defined(DEBUG)
    assert(t->on_cluster == m->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(t->height == m->width && t->width == m->height, "The t matrix should have height and width swapped w.r.t. m.");
#endif

    if (t->on_cluster) {
        transpose_par(t, m);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(t, j, i) = MAT_CELL(m, i, j);
            }
        }
    }
}

/*
 * Add scalar to matrix
 */
void add_scalar(const Matrix *s, const Matrix *m, const fp scalar) {
#if defined(DEBUG)
    assert(s->on_cluster == m->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(s->height == m->height && s->width == m->width, "The two matrices should have the same shape.");
#endif

    if (s->on_cluster) {
        add_scalar_par(s, m, scalar);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) + scalar;
            }
        }
    }
}

/*
 * Add scalar to matrix (in-place operation, it modifies the matrix)
 */
void add_scalar_(const Matrix *m, const fp scalar) {
    if (m->on_cluster) {
        add_scalar_par_(m, scalar);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(m, i, j) += scalar;
            }
        }
    }
}

/*
 * Scales a matrix by a scalar constant
 */
void scale(const Matrix *s, const Matrix *m, const fp scalar) {
#if defined(DEBUG)
    assert(s->on_cluster == m->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(s->height == m->height && s->width == m->width, "The two matrices should have the same shape.");
#endif
    if(s->on_cluster) {
        scale_par(s, m, scalar);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) * scalar;
            }
        }
    }
}

/*
 * Scales a matrix by a scalar constant (in-place operation, it modifies the matrix)
 */
void scale_(const Matrix *m, const fp scalar) {
    if(m->on_cluster) {
        scale_par_(m, scalar);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(m, i, j) *= scalar;
            }
        }
    }
}

/*
 * Perform element-wise addition between two matrices with same shape
 */
void add_mat(const Matrix *s, const Matrix *m1, const Matrix *m2) {
#if defined(DEBUG)
    assert(s->on_cluster == m1->on_cluster && m1->on_cluster == m2->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(s->height == m1->height && s->width == m1->width && m1->height == m2->height && m1->width == m2->width, "The three matrices should have the same shape.");
#endif

    if (s->on_cluster) {
        add_mat_par(s, m1, m2);  // compute it in parallel
    } else {
        for (int i = 0; i < m1->height; i++) {
            for (int j = 0; j < m1->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) + MAT_CELL(m2, i, j);
            }
        }
    }
}

/*
 * Perform element-wise addition between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void add_mat_(const Matrix *m1, const Matrix *m2) {
#if defined(DEBUG)
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(m1->height == m2->height && m1->width == m2->width, "The two matrices should have the same shape.");
#endif

    if (m1->on_cluster) {
        add_mat_par_(m1, m2);  // compute it in parallel
    } else {
        for (int i = 0; i < m1->height; i++) {
            for (int j = 0; j < m1->width; j++) {
                MAT_CELL(m1, i, j) += MAT_CELL(m2, i, j);
            }
        }
    }
}

/*
 * Perform row-wise addition between a matrix and a row vector with same width
 */
void add_row(const Matrix *s, const Matrix *m, const Matrix *r) {
#if defined(DEBUG)
    assert(s->on_cluster == m->on_cluster && m->on_cluster == r->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(s->height == m->height && s->width == m->width, "The s and m matrices should have the same shape.");
    assert(is_row_vector(r) && r->width == m->width, "The r matrix should be a row vector with the same width as m.");
#endif
    
    if (s->on_cluster) {
        add_row_par(s, m, r);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) + MAT_CELL(r, 0, i);
            }
        }
    }
}

/*
 * Perform row-wise addition between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void add_row_(const Matrix *m, const Matrix *r) {
#if defined(DEBUG)
    assert(m->on_cluster == r->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(is_row_vector(r) && r->width == m->width, "The r matrix should be a row vector with the same width as m.");
#endif

    if (m->on_cluster) {
        add_row_par_(m, r);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(m, i, j) += MAT_CELL(r, 0, i);
            }
        }
    }
}

/*
 * Perform column-wise addition between a matrix and a column vector with same height
 */
void add_col(const Matrix *s, const Matrix *m, const Matrix *c) {
#if defined(DEBUG)
    assert(s->on_cluster == m->on_cluster && m->on_cluster == c->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(s->height == m->height && s->width == m->width, "The s and m matrices should have the same shape.");
    assert(is_col_vector(c) && c->height == m->height, "The c matrix should be a column vector with the same height as m.");
#endif

    if (m->on_cluster) {
        add_col_par(s, m, s);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) + MAT_CELL(c, i, 0);
            }
        }
    }
}

/*
 * Perform column-wise addition between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void add_col_(const Matrix *m, const Matrix *c) {
#if defined(DEBUG)
    assert(m->on_cluster == c->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(is_col_vector(c) && c->height == m->height, "The c matrix should be a column vector with the same height as m.");
#endif

    if (m->on_cluster) {
        add_col_par_(m, c);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(m, i, j) += MAT_CELL(c, i, 0);
            }
        }
    }
}

/*
 * Perform element-wise subtraction between two matrices with same shape
 */
void sub_mat(const Matrix *s, const Matrix *m1, const Matrix *m2) {
#if defined(DEBUG)
    assert(s->on_cluster == m1->on_cluster && m1->on_cluster == m2->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(s->height == m1->height && s->width == m1->width && m1->height == m2->height && m1->width == m2->width, "The three matrices should have the same shape.");
#endif

    if (s->on_cluster) {
        sub_mat_par(s, m1, m2);  // compute it in parallel
    } else {
        for (int i = 0; i < m1->height; i++) {
            for (int j = 0; j < m1->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) - MAT_CELL(m2, i, j);
            }
        }
    }
}

/*
 * Perform element-wise subtraction between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void sub_mat_(const Matrix *m1, const Matrix *m2) {
#if defined(DEBUG)
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(m1->height == m2->height && m1->width == m2->width, "The two matrices should have the same shape.");
#endif

    if (m1->on_cluster) {
        sub_mat_par_(m1, m2);  // compute it in parallel
    } else {
        for (int i = 0; i < m1->height; i++) {
            for (int j = 0; j < m1->width; j++) {
                MAT_CELL(m1, i, j) -= MAT_CELL(m2, i, j);
            }
        }
    }
}

/*
 * Perform row-wise subtraction between a matrix and a row vector with same width
 */
void sub_row(const Matrix *s, const Matrix *m, const Matrix *r) {
#if defined(DEBUG)
    assert(s->on_cluster == m->on_cluster && m->on_cluster == r->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(s->height == m->height && s->width == m->width, "The s and m matrices should have the same shape.");
    assert(is_row_vector(r) && r->width == m->width, "The r matrix should be a row vector with the same width as m.");
#endif
    
    if (s->on_cluster) {
        sub_row_par(s, m, r);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) - MAT_CELL(r, 0, i);
            }
        }
    }
}

/*
 * Perform row-wise subtraction between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void sub_row_(const Matrix *m, const Matrix *r) {
#if defined(DEBUG)
    assert(m->on_cluster == r->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(is_row_vector(r) && r->width == m->width, "The r matrix should be a row vector with the same width as m.");
#endif

    if (m->on_cluster) {
        sub_row_par_(m, r);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(m, i, j) -= MAT_CELL(r, 0, i);
            }
        }
    }
}

/*
 * Perform column-wise subtraction between a matrix and a column vector with same height
 */
void sub_col(const Matrix *s, const Matrix *m, const Matrix *c) {
#if defined(DEBUG)
    assert(s->on_cluster == m->on_cluster && m->on_cluster == c->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(s->height == m->height && s->width == m->width, "The s and m matrices should have the same shape.");
    assert(is_col_vector(c) && c->height == m->height, "The c matrix should be a column vector with the same height as m.");
#endif

    if (s->on_cluster) {
        sub_col_par(s, m, c);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) - MAT_CELL(c, i, 0);
            }
        }
    }
}

/*
 * Perform column-wise subtraction between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void sub_col_(const Matrix *m, const Matrix *c) {
#if defined(DEBUG)
    assert(m->on_cluster == c->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(is_col_vector(c) && c->height == m->height, "The c matrix should be a column vector with the same height as m.");
#endif

    if (m->on_cluster) {
        sub_col_par_(m, c);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(m, i, j) -= MAT_CELL(c, i, 0);
            }
        }
    }
}

/*
 * Perform element-wise product (i.e. Hadamard) between two matrices with same shape
 */
void hadamard(const Matrix *s, const Matrix *m1, const Matrix *m2) {
#if defined(DEBUG)
    assert(s->on_cluster == m1->on_cluster && m1->on_cluster == m2->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(s->height == m1->height && s->width == m1->width && m1->height == m2->height && m1->width == m2->width, "The three matrices should have the same shape.");
#endif

    if (s->on_cluster) {
        hadamard_par(s, m1, m2);  // compute it in parallel
    } else {
        for (int i = 0; i < m1->height; i++) {
            for (int j = 0; j < m1->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) * MAT_CELL(m2, i, j);
            }
        }
    }
}

/*
 * Perform element-wise product (i.e. Hadamard) between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void hadamard_(const Matrix *m1, const Matrix *m2) {
#if defined(DEBUG)
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(m1->height == m2->height && m1->width == m2->width, "The two matrices should have the same shape.");
#endif

    if (m1->on_cluster) {
        hadamard_par_(m1, m2);  // compute it in parallel
    } else {
        for (int i = 0; i < m1->height; i++) {
            for (int j = 0; j < m1->width; j++) {
                MAT_CELL(m1, i, j) *= MAT_CELL(m2, i, j);
            }
        }
    }
}

/*
 * Perform row-wise product (i.e. Hadamard) between a matrix and a row vector with same width
 */
void hadamard_row(const Matrix *s, const Matrix *m, const Matrix *r) {
#if defined(DEBUG)
    assert(s->on_cluster == m->on_cluster && m->on_cluster == r->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(s->height == m->height && s->width == m->width, "The s and m matrices should have the same shape.");
    assert(is_row_vector(r) && r->width == m->width, "The r matrix should be a row vector with the same width as m.");
#endif
    
    if (s->on_cluster) {
        hadamard_row_par(s, m, r);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) * MAT_CELL(r, 0, i);
            }
        }
    }
}

/*
 * Perform row-wise product (i.e. Hadamard) between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void hadamard_row_(const Matrix *m, const Matrix *r) {
#if defined(DEBUG)
    assert(m->on_cluster == r->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(is_row_vector(r) && r->width == m->width, "The r matrix should be a row vector with the same width as m.");
#endif

    if (m->on_cluster) {
        hadamard_row_par_(m, r);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(m, i, j) *= MAT_CELL(r, 0, i);
            }
        }
    }
}

/*
 * Perform column-wise product (i.e. Hadamard) between a matrix and a column vector with same height
 */
void hadamard_col(const Matrix *s, const Matrix *m, const Matrix *c) {
#if defined(DEBUG)
    assert(s->on_cluster == m->on_cluster && m->on_cluster == c->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(s->height == m->height && s->width == m->width, "The s and m matrices should have the same shape.");
    assert(is_col_vector(c) && c->height == m->height, "The c matrix should be a column vector with the same height as m.");
#endif

    if (s->on_cluster) {
        hadamard_col_par(s, m, c);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) * MAT_CELL(c, i, 0);
            }
        }
    }
}

/*
 * Perform column-wise product (i.e. Hadamard) between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void hadamard_col_(const Matrix *m, const Matrix *c) {
#if defined(DEBUG)
    assert(m->on_cluster == c->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(is_col_vector(c) && c->height == m->height, "The c matrix should be a column vector with the same height as m.");
#endif

    if (m->on_cluster) {
        hadamard_col_par_(m, c);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(m, i, j) *= MAT_CELL(c, i, 0);
            }
        }
    }
}

/*
 * Perform element-wise division between two matrices with same shape
 */
void div_mat(const Matrix *s, const Matrix *m1, const Matrix *m2) {
#if defined(DEBUG)
    assert(s->on_cluster == m1->on_cluster && m1->on_cluster == m2->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(s->height == m1->height && s->width == m1->width && m1->height == m2->height && m1->width == m2->width, "The three matrices should have the same shape.");
#endif

    if (s->on_cluster) {
        div_mat_par(s, m1, m2);  // compute it in parallel
    } else {
        for (int i = 0; i < m1->height; i++) {
            for (int j = 0; j < m1->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) / MAT_CELL(m2, i, j);
            }
        }
    }
}

/*
 * Perform element-wise division between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void div_mat_(const Matrix *m1, const Matrix *m2) {
#if defined(DEBUG)
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(m1->height == m2->height && m1->width == m2->width, "The two matrices should have the same shape.");
#endif

    if (m1->on_cluster) {
        div_mat_par_(m1, m2);  // compute it in parallel
    } else {
        for (int i = 0; i < m1->height; i++) {
            for (int j = 0; j < m1->width; j++) {
                MAT_CELL(m1, i, j) /= MAT_CELL(m2, i, j);
            }
        }
    }
}

/*
 * Perform row-wise division between a matrix and a row vector with same width
 */
void div_row(const Matrix *s, const Matrix *m, const Matrix *r) {
#if defined(DEBUG)
    assert(s->on_cluster == m->on_cluster && m->on_cluster == r->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(s->height == m->height && s->width == m->width, "The s and m matrices should have the same shape.");
    assert(is_row_vector(r) && r->width == m->width, "The r matrix should be a row vector with the same width as m.");
#endif

    if (m->on_cluster) {
        div_row_par(s, m, r);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) / MAT_CELL(r, 0, i);
            }
        }
    }
}

/*
 * Perform row-wise division between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void div_row_(const Matrix *m, const Matrix *r) {
#if defined(DEBUG)
    assert(m->on_cluster == r->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(is_row_vector(r) && r->width == m->width, "The r matrix should be a row vector with the same width as m.");
#endif

    if (m->on_cluster) {
        div_row_par_(m, r);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(m, i, j) /= MAT_CELL(r, 0, i);
            }
        }
    }
}

/*
 * Perform column-wise division between a matrix and a column vector with same height
 */
void div_col(const Matrix *s, const Matrix *m, const Matrix *c) {
#if defined(DEBUG)
    assert(s->on_cluster == m->on_cluster && m->on_cluster == c->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(s->height == m->height && s->width == m->width, "The s and m matrices should have the same shape.");
    assert(is_col_vector(c) && c->height == m->height, "The c matrix should be a column vector with the same height as m.");
#endif

    if (m->on_cluster) {
        div_col_par(s, m, c);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) / MAT_CELL(c, i, 0);
            }
        }
    }
}

/*
 * Perform column-wise division between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void div_col_(const Matrix *m, const Matrix *c) {
#if defined(DEBUG)
    assert(m->on_cluster == c->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(is_col_vector(c) && c->height == m->height, "The c matrix should be a column vector with the same height as m.");
#endif

    if (m->on_cluster) {
        div_col_par_(m, c);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(m, i, j) /= MAT_CELL(c, i, 0);
            }
        }
    }
}

/*
 * Perform matrix multiplication between an AxB matrix and a BxC matrix
 */
void mat_mul(const Matrix *s, const Matrix *m1, const Matrix *m2) {
#if defined(DEBUG)
    assert(s->on_cluster == m1->on_cluster && m1->on_cluster == m2->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(m1->width == m2->height, "The width of the first matrix and the height of the second matrix should be equal.");
    assert(s->height == m1->height && s->width == m2->width, "The s matrix should have a compatible shape.");
#endif

    if (s->on_cluster) {
        if (pi_core_id() == 0) {
            // Initialize s data to zero
            memset(s->data, 0, s->height * s->width * sizeof(fp));
        }
        pi_cl_team_barrier();
        mat_mul_par(s, m1, m2);  // compute it in parallel
    } else {
        // Initialize s data to zero
        memset(s->data, 0, s->height * s->width * sizeof(fp));
        // i-k-j loop to optimize memory access
        for (int i = 0; i < m1->height; i++) {
            for (int k = 0; k < m1->width; k++) {
                for (int j = 0; j < m2->width; j++) {
                    MAT_CELL(s, i, j) += MAT_CELL(m1, i, k) * MAT_CELL(m2, k, j);
                }
            }
        }
    }
}

/*
 * Perform matrix multiplication between a BxA matrix and a BxC matrix, without transposing the first matrix
 */
void mat_mul_t1(const Matrix *s, const Matrix *m1, const Matrix *m2) {
#if defined(DEBUG)
    assert(s->on_cluster == m1->on_cluster && m1->on_cluster == m2->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(m1->height == m2->height, "The height of the two matrices should be equal.");
    assert(s->height == m1->width && s->width == m2->width, "The s matrix should have a compatible shape.");
#endif

    if (s->on_cluster) {
        if (pi_core_id() == 0) {
            // Initialize s data to zero
            memset(s->data, 0, s->height * s->width * sizeof(fp));
        }
        pi_cl_team_barrier();
        mat_mul_t1_par(s, m1, m2);  // compute it in parallel
    } else {
        // Initialize s data to zero
        memset(s->data, 0, s->height * s->width * sizeof(fp));
        // k-i-j loop to optimize memory access
        for (int k = 0; k < m1->height; k++) {
            for (int i = 0; i < m1->width; i++) {
                for (int j = 0; j < m2->width; j++) {
                    MAT_CELL(s, i, j) += MAT_CELL(m1, k, i) * MAT_CELL(m2, k, j);
                }
            }
        }
    }
}

/*
 * Perform matrix multiplication between an AxB matrix and a CxB matrix, without transposing the second matrix
 */
void mat_mul_t2(const Matrix *s, const Matrix *m1, const Matrix *m2) {
#if defined(DEBUG)
    assert(s->on_cluster == m1->on_cluster && m1->on_cluster == m2->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(m1->width == m2->width, "The width of the two matrices should be equal.");
    assert(s->height == m1->height && s->width == m2->height, "The s matrix should have a compatible shape.");
#endif

    if (s->on_cluster) {
        if (pi_core_id() == 0) {
            // Initialize s data to zero
            memset(s->data, 0, s->height * s->width * sizeof(fp));
        }
        pi_cl_team_barrier();
        mat_mul_t2_par(s, m1, m2);  // compute it in parallel
    } else {
        // Initialize s data to zero
        memset(s->data, 0, s->height * s->width * sizeof(fp));
        // i-j-k loop to optimize memory access
        for (int i = 0; i < m1->height; i++) {
            for (int j = 0; j < m2->height; j++) {
                for (int k = 0; k < m1->width; k++) {
                    MAT_CELL(s, i, j) += MAT_CELL(m1, i, k) * MAT_CELL(m2, j, k);
                }
            }
        }
    }
}

/*
 * Perform dot product between two vectors
 */
fp dot(const Matrix *v1, const Matrix *v2) {
#if defined(DEBUG)
    assert(v1->on_cluster == v2->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(is_vector(v1) && is_vector(v2), "The two matrices should be vectors (i.e. either their height or width should be equal to 1).");
    assert(MAX(v1->height, v1->width) == MAX(v2->height, v2->width), "The two vectors should have the same length.");
#endif

    fp acc = 0;
    if (v1->on_cluster) {
        acc = dot_par(v1, v2);  // compute it in parallel
    } else {
        const uint len = v1->height == 1 ? v1->width : v1->height;
        for (int i = 0; i < len; i++) {
            acc += MAT_DATA(v1)[i] * MAT_DATA(v2)[i];
        }
    }

    return acc;
}

/*
 * Perform outer product between two vectors
 */
void outer(const Matrix *s, const Matrix *v1, const Matrix *v2) {
#if defined(DEBUG)
    assert(s->on_cluster == v1->on_cluster && v1->on_cluster == v2->on_cluster, "The three matrices should be allocated on the same memory.");
    assert(is_vector(v1) && is_vector(v2), "The v1 and v2 matrices should be vectors (i.e. either their height or width should be equal to 1).");
    assert(s->height == MAX(v1->height, v1->width) && s->width == MAX(v2->height, v2->width), "The s matrix should have a compatible shape.");
#endif

    if (v1->on_cluster) {
        if (pi_core_id() == 0) {
            // Initialize s data to zero
            memset(s->data, 0, s->height * s->width * sizeof(fp));
        }
        pi_cl_team_barrier();
        outer_par(s, v1, v2);  // compute it in parallel
    } else {
        // Initialize s data to zero
        memset(s->data, 0, s->height * s->width * sizeof(fp));
        const uint len1 = MAX(v1->height, v1->width);
        const uint len2 = MAX(v2->height, v2->width);
        for (int i = 0; i < len1; i++) {
            for (int j = 0; j < len2; j++) {
                MAT_CELL(s, i, j) = MAT_DATA(v1)[i] * MAT_DATA(v2)[j];
            }
        }
    }
}

/*
 * Check if two matrices are equal
 */
bool are_equal(const Matrix *m1, const Matrix *m2, const fp tol) {
    if (m1->height != m2->height || m1->width != m2->width)
        return false;

    for (int i = 0; i < m1->height; i++) {
        for (int j = 0; j < m1->width; j++) {
            if (MAT_CELL(m1, i, j) - MAT_CELL(m2, i, j) > tol)
                return false;
        }
    }
    return true;
}

/*
 * Check if a matrix is square
 */
bool is_square(const Matrix *m) {
    return m->height == m->width;
}

/*
 * Check if a matrix is a vector
 */
bool is_vector(const Matrix *m) {
    return m->height == 1 || m->width == 1;
}

/*
 * Check if a matrix is a row vector
 */
bool is_row_vector(const Matrix *m) {
    return m->height == 1;
}

/*
 * Check if a matrix is a column vector
 */
bool is_col_vector(const Matrix *m) {
    return m->width == 1;
}

/*
 * Get the main diagonal of a given matrix and return a new column vector (allocated in the same memory)
 */
Matrix diagonal(const Matrix *m) {
    // Create matrix d
    Matrix d = {
        .data = m->data,
        .height = m->height,
        .width = 1,  // column vector
        .offset = m->width + 1,  // same offset as m + 1
        .on_cluster = m->on_cluster
    };

    return d;
}

/*
 * Set cells below the main diagonal to zero (in-place operation, it modifies the matrix)
 */
void tri_up(const Matrix *m) {
    if (m->on_cluster) {
        tri_up_par(m);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < i; j++) {
                MAT_CELL(m, i, j) = 0;
            }
        }
    }
}

/*
 * Slice given matrix and return a new matrix (allocated in the same memory)
 */
Matrix slice(const Matrix *m, const uint row_start, const uint row_stop, const uint col_start, const uint col_stop) {
    uint row_span = row_stop - row_start + 1;
    uint col_span = col_stop - col_start + 1;
    
#if defined(DEBUG)
    assert(0 <= row_start && row_start < m->height, "The row_start argument is out of bounds.");
    assert(0 <= row_stop && row_stop < m->height, "The row_stop argument is out of bounds.");
    assert(0 <= col_start && col_start < m->width, "The col_start argument is out of bounds.");
    assert(0 <= col_stop && col_stop < m->width, "The col_stop argument is out of bounds.");
    assert(row_start <= row_stop, "The row_stop argument must be greater than or equal to row_start.");
    assert(col_start <= col_stop, "The col_stop argument must be greater than or equal to col_start.");
    assert(m->width == m->offset, "The width of the matrix must be equal to its offset.");
#endif

    // Get pointer to (row_start, col_start) cell
    fp *new_pt = MAT_DATA(m) + MAT_IDX(m, row_start, col_start);
    // Create matrix s
    Matrix s = {
        .data = new_pt,
        .height = row_span,
        .width = col_span,
        .offset = m->width,  // same offset as m
        .on_cluster = m->on_cluster
    };

    return s;
}

/*
 * Extract the k-th row, creating a new vector (allocated in the same memory)
 */
Matrix extract_row(const Matrix *m, const uint k) {
#if defined(DEBUG)
    assert(0 <= k && k < m->height, "Index is out of bounds for rows.");
#endif
    
    return slice(m, k, k, 0, m->width - 1);
}

/*
 * Extract the k-th column, creating a new vector (allocated in the same memory)
 */
Matrix extract_col(const Matrix *m, const uint k) {
#if defined(DEBUG)
    assert(0 <= k && k < m->width, "Index is out of bounds for columns.");
#endif

    return slice(m, 0, m->height - 1, k, k);
}

/*
 * Compact a matrix obtained by slicing, namely a matrix whose offset is different from its width (it modifies the original matrix)
 */
void compact(Matrix *m, const uint row_start, const uint col_start, uint const orig_len) {
#if defined(DEBUG)
    assert(m->width != m->offset, "The width of the matrix must be different from its offset.");
#endif

    // Get original pointer
    fp *old_pt = MAT_DATA(m) - MAT_IDX(m, row_start, col_start) + 1;

    if (pi_core_id() == 0) {
        fp *cur_pt = old_pt;
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                *cur_pt = MAT_CELL(m, i, j);
                cur_pt++;
            }
        }

        // Set remaining memory to zero
        memset(cur_pt, 0, (orig_len - m->height * m->width) * sizeof(fp));
    }

    pi_cl_team_barrier();

    // Update struct
    m->data = old_pt;
    m->offset = m->width;
}

/*
 * Prints a matrix to standard output
 */
void print_mat(const Matrix *m) {
    for (int i = 0; i < m->height; i++) {
        for (int j = 0; j < m->width; j++) {
            printf("%.5f ", MAT_CELL(m, i, j));
        }
        printf("\n");
    }
}

/*
 * Write a matrix to a binary file
 */
/*
void write_mat(const char *path, const Matrix *m) {
    // Open file
    FILE *file = fopen(path, "wb");
    assert(file != NULL, "Cannot open or create file.");

    // Write number of dimensions, followed by height and width
    int n_dim = 2;
    assert(fwrite(&n_dim, sizeof(int), 1, file) == 1, "Could not write to file.");
    assert(fwrite(&m->height, sizeof(int), 1, file) == 1, "Could not write to file.");
    assert(fwrite(&m->width, sizeof(int), 1, file) == 1, "Could not write to file.");

    // Write flattened matrix data
    int size = m->height * m->width;
    assert(fwrite(m->data, sizeof(fp), size, file) == size, "Could not write to file.");

    // Close file
    assert(fclose(file) == 0, "Could not close file.");
}
*/

/*
 * Copy a matrix into a new matrix
 */
void copy_mat(const Matrix *s, const Matrix *m) {
#if defined(DEBUG)
    assert(s->on_cluster == m->on_cluster, "The two matrices should be allocated on the same memory.");
    assert(s->height == m->height && s->width == m->width, "The two matrices should have the same shape.");
#endif

    if (s->on_cluster) {
        copy_mat_par(s, m);  // compute it in parallel
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j);
            }
        }
    }
}
