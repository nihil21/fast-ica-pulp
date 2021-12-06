//
// Created by nihil on 20/09/21.
//

#include <string.h>
#include "../include/matrix.h"
#include "../include/utils.h"
#include "../include/random.h"
#include "../include/cluster.h"

extern struct pi_device cluster_dev;

/*
 * Initialize a zero matrix given its dimensions
 */
Matrix *new_mat(uint height, uint width, bool on_cluster) {
    assert(height > 0 && width > 0, "Matrix height and width should be greater than 0.");

    Matrix *m;
    if (on_cluster) {
        // Allocate memory for the struct (Cluster's L1 memory)
        m = pi_l1_malloc(&cluster_dev, sizeof(Matrix));
        assert(m != NULL, "Could not allocate Matrix struct (Cluster's L1 memory).");

        // Allocate memory for the array field (Cluster's L1 memory)
        m->data = pi_l1_malloc(&cluster_dev, height * width * sizeof(fp));
        assert(m->data != NULL, "Could not allocate Matrix data (Cluster's L1 memory).");
        memset(m->data, 0, height * width * sizeof(fp));
    } else {
        // Allocate memory for the struct (Fabric Controller's L2 memory)
        m = pi_l2_malloc(sizeof(Matrix));
        assert(m != NULL, "Could not allocate Matrix struct (Fabric Controller's L2 memory).");

        // Allocate memory for the array field (Fabric Controller's L2 memory)
        m->data = pi_l2_malloc(height * width * sizeof(fp));
        assert(m->data != NULL, "Could not allocate Matrix data (Fabric Controller's L2 memory).");
        memset(m->data, 0, height * width * sizeof(fp));  // initialize to zero
    }

    // Set other fields
    *(uint *) &m->height = height;  // cast away constness
    *(uint *) &m->width = width;  // cast away constness
    m->on_cluster = on_cluster;

    return m;
}

/*
 * Initialize a zero column vector given its length
 */
Matrix *new_vec(uint length, bool on_cluster) {
    return new_mat(length, 1, on_cluster);
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
void mat_fc2cl(Matrix *m_fc, Matrix *m_cl) {
    // Transfer data field
    pi_cl_dma_cmd_t dma_cmd_data;
    pi_cl_dma_cmd((uint32_t) MAT_DATA(m_fc), (uint32_t) MAT_DATA(m_cl), m_fc->height * m_fc->width * sizeof(fp), PI_CL_DMA_DIR_EXT2LOC, &dma_cmd_data);
    pi_cl_dma_cmd_wait(&dma_cmd_data);
}

/*
 * Move matrices from Fabric Controller to Cluster (async)
 */
void mat_fc2cl_async(Matrix *m_fc, Matrix *m_cl, pi_cl_dma_cmd_t *dma_cmd_data) {
    // Transfer data field
    pi_cl_dma_cmd((uint32_t) MAT_DATA(m_fc), (uint32_t) MAT_DATA(m_cl), m_fc->height * m_fc->width * sizeof(fp), PI_CL_DMA_DIR_EXT2LOC, dma_cmd_data);
}

/*
 * Move matrices from Cluster to Fabric Controller
 */
void mat_cl2fc(Matrix *m_cl, Matrix *m_fc) {
    // Transfer data field
    pi_cl_dma_cmd_t dma_cmd_data;
    pi_cl_dma_cmd((uint32_t) MAT_DATA(m_fc), (uint32_t) MAT_DATA(m_cl), m_cl->height * m_cl->width * sizeof(fp), PI_CL_DMA_DIR_LOC2EXT, &dma_cmd_data);
    pi_cl_dma_cmd_wait(&dma_cmd_data);
}

/*
 * Move matrices from Cluster to Fabric Controller (async)
 */
void mat_cl2fc_async(Matrix *m_cl, Matrix *m_fc, pi_cl_dma_cmd_t *dma_cmd_data) {
    // Transfer data field
    pi_cl_dma_cmd((uint32_t) MAT_DATA(m_fc), (uint32_t) MAT_DATA(m_cl), m_cl->height * m_cl->width * sizeof(fp), PI_CL_DMA_DIR_LOC2EXT, dma_cmd_data);
}

/*
 * Allocate the identity matrix
 */
Matrix *eye(uint n, bool on_cluster) {
    Matrix *i = new_mat(n, n, on_cluster);

    for (int k = 0; k < n; k++) {
        MAT_CELL(i, k, k) = 1;
    }

    return i;
}

/*
 * Allocate a matrix with the given dimensions and fill it with random uniform integers in a given range
 */
Matrix *mat_randint(uint height, uint width, int min, int max, bool on_cluster) {
    Matrix *m = new_mat(height, width, on_cluster);

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            MAT_CELL(m, i, j) = (fp) uniform_randint_range(min, max);
        }
    }

    return m;
}

/*
 * Allocate a matrix with the given dimensions and fill it with random uniform numbers in a given range
 */
Matrix *mat_rand(uint height, uint width, fp min, fp max, bool on_cluster) {
    Matrix *m = new_mat(height, width, on_cluster);

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            MAT_CELL(m, i, j) = uniform_rand_range(min, max);
        }
    }

    return m;
}

/*
 * Allocate a matrix with the given dimensions and fill it with random normal numbers
 */
Matrix *mat_randn(uint height, uint width, bool on_cluster) {
    Matrix *m = new_mat(height, width, on_cluster);

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            MAT_CELL(m, i, j) = standard_rand();
        }
    }

    return m;
}

/*
 * Build a linear space with the given range and number of samples
 */
Matrix *linspace(fp start, fp stop, int n_samples, bool on_cluster) {
    assert(start < stop, "The stop argument should be greater than the start argument.");

    Matrix *ls = new_vec(n_samples, on_cluster);

    fp step = (stop - start) / (fp) (n_samples - 1);
    for (int i = 0; i < n_samples; i++) {
        MAT_CELL(ls, i, 0) = start + step * (fp) i;
    }

    return ls;
}

/*
 * Create a matrix with given height and width from an array of data
 */
Matrix *from_array(const fp *data, uint height, uint width, bool on_cluster) {
    Matrix *s = new_mat(height, width, on_cluster);

    int d_idx = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            MAT_CELL(s, i, j) = data[d_idx++];
        }
    }

    return s;
}

/*
 * Compute the L2 norm of a given matrix
 */
fp norm(Matrix *m) {
    fp acc = 0;
    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMSS arg = {m, NULL, &acc};
        // Compute norm in parallel
        pi_cl_team_fork(CORES, norm_par, (void *) &arg);
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
fp mean(Matrix *m) {
    fp acc = 0;
    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMSS arg = {m, NULL, &acc};
        // Compute mean in parallel
        pi_cl_team_fork(CORES, mean_par, (void *) &arg);
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
Matrix *row_mean(Matrix *m) {
    Matrix *r = new_mat(1, m->width, m->on_cluster);  // row vector

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, NULL, r};
        // Compute mean in parallel
        pi_cl_team_fork(CORES, row_mean_par, (void *) &arg);
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

    return r;
}

/*
 * Get the mean of a given matrix along columns
 */
Matrix *col_mean(Matrix *m) {
    Matrix *c = new_vec(m->height, m->on_cluster);  // column vector

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, NULL, c};
        // Compute mean in parallel
        pi_cl_team_fork(CORES, col_mean_par, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            fp acc = 0;
            for (int j = 0; j < m->width; j++) {
                acc += MAT_CELL(m, i, j);
            }
            MAT_CELL(c, i, 0) = acc / (fp) m->width;
        }
    }

    return c;
}

/*
 * Get the standard deviation of a given matrix
 */
fp std(Matrix *m) {
    fp mean_ = mean(m);
    fp acc = 0;
    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMSS arg = {m, &mean_, &acc};
        // Compute standard deviation in parallel
        pi_cl_team_fork(CORES, std_par, (void *) &arg);
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
Matrix *row_std(Matrix *m) {
    Matrix *r = new_mat(1, m->width, m->on_cluster);  // row vector
    Matrix *mean_ = row_mean(m);

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, mean_, r};
        // Compute standard deviation in parallel
        pi_cl_team_fork(CORES, row_std_par, (void *) &arg);
    } else {
        // Memory access not contiguous
        for (int j = 0; j < m->width; j++) {
            fp acc = 0;
            for (int i = 0; i < m->height; i++) {
                acc += POW(MAT_CELL(m, i, j) - MAT_CELL(mean_, 0, j), 2);
            }
            MAT_CELL(r, 0, j) = SQRT(acc / (fp) m->height);
        }
    }

    free_mat(mean_);

    return r;
}

/*
 * Get the standard deviation of a given matrix along columns
 */
Matrix *col_std(Matrix *m) {
    Matrix *c = new_vec(m->height, m->on_cluster);  // column vector
    Matrix *mean_ = col_mean(m);

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, mean_, c};
        // Compute standard deviation in parallel
        pi_cl_team_fork(CORES, col_std_par, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            fp acc = 0;
            for (int j = 0; j < m->width; j++) {
                acc += POW(MAT_CELL(m, i, j) - MAT_CELL(mean_, i, 0), 2);
            }
            MAT_CELL(c, i, 0) = SQRT(acc / (fp) m->width);
        }
    }

    free_mat(mean_);

    return c;
}

/*
 * Transpose a matrix
 */
Matrix *transpose(Matrix* m) {
    Matrix *t = new_mat(m->width, m->height, m->on_cluster);

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, NULL, t};
        // Compute transpose in parallel
        pi_cl_team_fork(CORES, transpose_par, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(t, j, i) = MAT_CELL(m, i, j);
            }
        }
    }

    return t;
}

/*
 * Scales a matrix by a scalar constant
 */
Matrix *scale(Matrix *m, fp scalar) {
    Matrix *s = new_mat(m->height, m->width, m->on_cluster);

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMSM arg = {m, &scalar, s};
        // Compute scaled matrix in parallel
        pi_cl_team_fork(CORES, scale_par, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) * scalar;
            }
        }
    }

    return s;
}

/*
 * Scales a matrix by a scalar constant (in-place operation, it modifies the matrix)
 */
void scale_(Matrix *m, fp scalar) {
    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMSM arg = {m, &scalar, NULL};
        // Compute scaled matrix in parallel
        pi_cl_team_fork(CORES, scale_par_, (void *) &arg);
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
Matrix *add_mat(Matrix *m1, Matrix *m2) {
    assert(m1->height == m2->height && m1->width == m2->width, "The two matrices should have the same shape.");
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m1->height, m1->width, m1->on_cluster);

    if (m1->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m1, m2, s};
        // Perform addition in parallel
        pi_cl_team_fork(CORES, add_mat_par, (void *) &arg);
    } else {
        for (int i = 0; i < m1->height; i++) {
            for (int j = 0; j < m1->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) + MAT_CELL(m2, i, j);
            }
        }
    }

    return s;
}

/*
 * Perform element-wise addition between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void add_mat_(Matrix *m1, Matrix *m2) {
    assert(m1->height == m2->height && m1->width == m2->width, "The two matrices should have the same shape.");
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");

    if (m1->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m1, m2, NULL};
        // Perform addition in parallel
        pi_cl_team_fork(CORES, add_mat_par_, (void *) &arg);
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
Matrix *add_row(Matrix *m, Matrix *r) {
    assert(is_row_vector(r), "The second matrix should be a row vector (i.e. with height equal to 1).");
    assert(r->width == m->width, "The width of the first matrix and the length of the row vector (i.e. its width) should be the same.");
    assert(m->on_cluster == r->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m->height, m->width, m->on_cluster);

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, r, s};
        // Perform addition in parallel
        pi_cl_team_fork(CORES, add_row_par, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) + MAT_CELL(r, 0, i);
            }
        }
    }

    return s;
}

/*
 * Perform row-wise addition between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void add_row_(Matrix *m, Matrix *r) {
    assert(is_row_vector(r), "The second matrix should be a row vector (i.e. with height equal to 1).");
    assert(r->width == m->width, "The width of the first matrix and the length of the row vector (i.e. its width) should be the same.");
    assert(m->on_cluster == r->on_cluster, "The two matrices should be allocated on the same memory.");

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, r, NULL};
        // Perform addition in parallel
        pi_cl_team_fork(CORES, add_row_par_, (void *) &arg);
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
Matrix *add_col(Matrix *m, Matrix *c) {
    assert(is_col_vector(c), "The second matrix should be a column vector (i.e. with width equal to 1).");
    assert(c->height == m->height, "The height of the first matrix and the length of the column vector (i.e. its height) should be the same.");
    assert(m->on_cluster == c->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m->height, m->width, m->on_cluster);

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, c, s};
        // Perform addition in parallel
        pi_cl_team_fork(CORES, add_col_par, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) + MAT_CELL(c, i, 0);
            }
        }
    }

    return s;
}

/*
 * Perform column-wise addition between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void add_col_(Matrix *m, Matrix *c) {
    assert(is_col_vector(c), "The second matrix should be a column vector (i.e. with width equal to 1).");
    assert(c->height == m->height, "The height of the first matrix and the length of the column vector (i.e. its height) should be the same.");
    assert(m->on_cluster == c->on_cluster, "The two matrices should be allocated on the same memory.");

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, c, NULL};
        // Perform addition in parallel
        pi_cl_team_fork(CORES, add_col_par_, (void *) &arg);
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
Matrix *sub_mat(Matrix *m1, Matrix *m2) {
    assert(m1->height == m2->height && m1->width == m2->width, "The two matrices should have the same shape.");
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m1->height, m1->width, m1->on_cluster);

    if (m1->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m1, m2, s};
        // Perform subtraction in parallel
        pi_cl_team_fork(CORES, sub_mat_par, (void *) &arg);
    } else {
        for (int i = 0; i < m1->height; i++) {
            for (int j = 0; j < m1->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) - MAT_CELL(m2, i, j);
            }
        }
    }

    return s;
}

/*
 * Perform element-wise subtraction between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void sub_mat_(Matrix *m1, Matrix *m2) {
    assert(m1->height == m2->height && m1->width == m2->width, "The two matrices should have the same shape.");
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");

    if (m1->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m1, m2, NULL};
        // Perform subtraction in parallel
        pi_cl_team_fork(CORES, sub_mat_par_, (void *) &arg);
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
Matrix *sub_row(Matrix *m, Matrix *r) {
    assert(is_row_vector(r), "The second matrix should be a row vector (i.e. with height equal to 1).");
    assert(r->width == m->width, "The width of the first matrix and the length of the row vector (i.e. its width) should be the same.");
    assert(m->on_cluster == r->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m->height, m->width, m->on_cluster);

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, r, s};
        // Perform subtraction in parallel
        pi_cl_team_fork(CORES, sub_row_par, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) - MAT_CELL(r, 0, i);
            }
        }
    }

    return s;
}

/*
 * Perform row-wise subtraction between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void sub_row_(Matrix *m, Matrix *r) {
    assert(is_row_vector(r), "The second matrix should be a row vector (i.e. with height equal to 1).");
    assert(r->width == m->width, "The width of the first matrix and the length of the row vector (i.e. its width) should be the same.");
    assert(m->on_cluster == r->on_cluster, "The two matrices should be allocated on the same memory.");

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, r, NULL};
        // Perform subtraction in parallel
        pi_cl_team_fork(CORES, sub_row_par_, (void *) &arg);
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
Matrix *sub_col(Matrix *m, Matrix *c) {
    assert(is_col_vector(c), "The second matrix should be a column vector (i.e. with width equal to 1).");
    assert(c->height == m->height, "The height of the first matrix and the length of the column vector (i.e. its height) should be the same.");
    assert(m->on_cluster == c->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m->height, m->width, m->on_cluster);

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, c, s};
        // Perform subtraction in parallel
        pi_cl_team_fork(CORES, sub_col_par, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) - MAT_CELL(c, i, 0);
            }
        }
    }

    return s;
}

/*
 * Perform column-wise subtraction between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void sub_col_(Matrix *m, Matrix *c) {
    assert(is_col_vector(c), "The second matrix should be a column vector (i.e. with width equal to 1).");
    assert(c->height == m->height, "The height of the first matrix and the length of the column vector (i.e. its height) should be the same.");
    assert(m->on_cluster == c->on_cluster, "The two matrices should be allocated on the same memory.");

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, c, NULL};
        // Perform subtraction in parallel
        pi_cl_team_fork(CORES, sub_col_par_, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(m, i, j) -= MAT_CELL(c, i, 0);
            }
        }
    }
}

/*
 * Add scalar to matrix
 */
Matrix *add_scalar(Matrix *m, fp scalar) {
    Matrix *s = new_mat(m->height, m->width, m->on_cluster);

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMSM arg = {m, &scalar, s};
        // Perform addition in parallel
        pi_cl_team_fork(CORES, add_scalar_par, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) + scalar;
            }
        }
    }

    return s;
}

/*
 * Add scalar to matrix (in-place operation, it modifies the matrix)
 */
void add_scalar_(Matrix *m, fp scalar) {
    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMSM arg = {m, &scalar, NULL};
        // Perform addition in parallel
        pi_cl_team_fork(CORES, add_scalar_par_, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(m, i, j) += scalar;
            }
        }
    }
}

/*
 * Perform element-wise product (i.e. Hadamard) between two matrices with same shape
 */
Matrix *hadamard(Matrix *m1, Matrix *m2) {
    assert(m1->height == m2->height && m1->width == m2->width, "The two matrices should have the same shape.");
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m1->height, m1->width, m1->on_cluster);

    if (m1->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m1, m2, s};
        // Perform multiplication in parallel
        pi_cl_team_fork(CORES, hadamard_par, (void *) &arg);
    } else {
        for (int i = 0; i < m1->height; i++) {
            for (int j = 0; j < m1->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) * MAT_CELL(m2, i, j);
            }
        }
    }

    return s;
}

/*
 * Perform element-wise product (i.e. Hadamard) between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void hadamard_(Matrix *m1, Matrix *m2) {
    assert(m1->height == m2->height && m1->width == m2->width, "The two matrices should have the same shape.");
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");

    if (m1->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m1, m2, NULL};
        // Perform multiplication in parallel
        pi_cl_team_fork(CORES, hadamard_par_, (void *) &arg);
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
Matrix *hadamard_row(Matrix *m, Matrix *r) {
    assert(is_row_vector(r), "The second matrix should be a row vector (i.e. with height equal to 1).");
    assert(r->width == m->width, "The width of the first matrix and the length of the row vector (i.e. its width) should be the same.");
    assert(m->on_cluster == r->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m->height, m->width, m->on_cluster);

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, r, s};
        // Perform multiplication in parallel
        pi_cl_team_fork(CORES, hadamard_row_par, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) * MAT_CELL(r, 0, i);
            }
        }
    }

    return s;
}

/*
 * Perform row-wise product (i.e. Hadamard) between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void hadamard_row_(Matrix *m, Matrix *r) {
    assert(is_row_vector(r), "The second matrix should be a row vector (i.e. with height equal to 1).");
    assert(r->width == m->width, "The width of the first matrix and the length of the row vector (i.e. its width) should be the same.");
    assert(m->on_cluster == r->on_cluster, "The two matrices should be allocated on the same memory.");

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, r, NULL};
        // Perform multiplication in parallel
        pi_cl_team_fork(CORES, hadamard_row_par_, (void *) &arg);
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
Matrix *hadamard_col(Matrix *m, Matrix *c) {
    assert(is_col_vector(c), "The second matrix should be a column vector (i.e. with width equal to 1).");
    assert(c->height == m->height, "The height of the first matrix and the length of the column vector (i.e. its height) should be the same.");
    assert(m->on_cluster == c->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m->height, m->width, m->on_cluster);

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, c, s};
        // Perform multiplication in parallel
        pi_cl_team_fork(CORES, hadamard_col_par, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) * MAT_CELL(c, i, 0);
            }
        }
    }

    return s;
}

/*
 * Perform column-wise product (i.e. Hadamard) between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void hadamard_col_(Matrix *m, Matrix *c) {
    assert(is_col_vector(c), "The second matrix should be a column vector (i.e. with width equal to 1).");
    assert(c->height == m->height, "The height of the first matrix and the length of the column vector (i.e. its height) should be the same.");
    assert(m->on_cluster == c->on_cluster, "The two matrices should be allocated on the same memory.");

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, c, NULL};
        // Perform multiplication in parallel
        pi_cl_team_fork(CORES, hadamard_col_par_, (void *) &arg);
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
Matrix *div_mat(Matrix *m1, Matrix *m2) {
    assert(m1->height == m2->height && m1->width == m2->width, "The two matrices should have the same shape.");
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m1->height, m1->width, m1->on_cluster);
    
    if (m1->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m1, m2, s};
        // Perform division in parallel
        pi_cl_team_fork(CORES, div_mat_par, (void *) &arg);
    } else {
        for (int i = 0; i < m1->height; i++) {
            for (int j = 0; j < m1->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m1, i, j) / MAT_CELL(m2, i, j);
            }
        }
    }

    return s;
}

/*
 * Perform element-wise division between two matrices with same shape (in-place operation, it modifies the first matrix)
 */
void div_mat_(Matrix *m1, Matrix *m2) {
    assert(m1->height == m2->height && m1->width == m2->width, "The two matrices should have the same shape.");
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");

    if (m1->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m1, m2, NULL};
        // Perform division in parallel
        pi_cl_team_fork(CORES, div_mat_par_, (void *) &arg);
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
Matrix *div_row(Matrix *m, Matrix *r) {
    assert(is_row_vector(r), "The second matrix should be a row vector (i.e. with height equal to 1).");
    assert(r->width == m->width, "The width of the first matrix and the length of the row vector (i.e. its width) should be the same.");
    assert(m->on_cluster == r->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m->height, m->width, m->on_cluster);

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, r, s};
        // Perform division in parallel
        pi_cl_team_fork(CORES, div_row_par, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) / MAT_CELL(r, 0, i);
            }
        }
    }

    return s;
}

/*
 * Perform row-wise division between a matrix and a row vector with same width (in-place operation, it modifies the first matrix)
 */
void div_row_(Matrix *m, Matrix *r) {
    assert(is_row_vector(r), "The second matrix should be a row vector (i.e. with height equal to 1).");
    assert(r->width == m->width, "The width of the first matrix and the length of the row vector (i.e. its width) should be the same.");
    assert(m->on_cluster == r->on_cluster, "The two matrices should be allocated on the same memory.");

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, r, NULL};
        // Perform division in parallel
        pi_cl_team_fork(CORES, div_row_par_, (void *) &arg);
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
Matrix *div_col(Matrix *m, Matrix *c) {
    assert(is_col_vector(c), "The second matrix should be a column vector (i.e. with width equal to 1).");
    assert(c->height == m->height, "The height of the first matrix and the length of the column vector (i.e. its height) should be the same.");
    assert(m->on_cluster == c->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m->height, m->width, m->on_cluster);

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, c, s};
        // Perform division in parallel
        pi_cl_team_fork(CORES, div_col_par, (void *) &arg);
    } else {
        for (int i = 0; i < m->height; i++) {
            for (int j = 0; j < m->width; j++) {
                MAT_CELL(s, i, j) = MAT_CELL(m, i, j) / MAT_CELL(c, i, 0);
            }
        }
    }

    return s;
}

/*
 * Perform column-wise division between a matrix and a column vector with same height (in-place operation, it modifies the first matrix)
 */
void div_col_(Matrix *m, Matrix *c) {
    assert(is_col_vector(c), "The second matrix should be a column vector (i.e. with width equal to 1).");
    assert(c->height == m->height, "The height of the first matrix and the length of the column vector (i.e. its height) should be the same.");
    assert(m->on_cluster == c->on_cluster, "The two matrices should be allocated on the same memory.");

    if (m->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m, c, NULL};
        // Perform division in parallel
        pi_cl_team_fork(CORES, div_col_par_, (void *) &arg);
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
Matrix *mat_mul(Matrix *m1, Matrix *m2) {
    assert(m1->width == m2->height, "The width of the first matrix and the height of the second matrix should be equal.");
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m1->height, m2->width, m1->on_cluster);

    if (m1->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m1, m2, s};
        // Perform GEMM in parallel
        pi_cl_team_fork(CORES, mat_mul_par, (void *) &arg);
    } else {
        // i-k-j loop to optimize memory access
        for (int i = 0; i < m1->height; i++) {
            for (int k = 0; k < m1->width; k++) {
                for (int j = 0; j < m2->width; j++) {
                    MAT_CELL(s, i, j) += MAT_CELL(m1, i, k) * MAT_CELL(m2, k, j);
                }
            }
        }
    }

    return s;
}

/*
 * Perform matrix multiplication between a BxA matrix and a BxC matrix, without transposing the first matrix
 */
Matrix *mat_mul_t1(Matrix* m1, Matrix* m2) {
    assert(m1->height == m2->height, "The height of the two matrices should be equal.");
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m1->width, m2->width, m1->on_cluster);

    if (m1->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m1, m2, s};
        // Perform GEMM in parallel
        pi_cl_team_fork(CORES, mat_mul_t1_par, (void *) &arg);
    } else {
        // k-i-j loop to optimize memory access
        for (int k = 0; k < m1->height; k++) {
            for (int i = 0; i < m1->width; i++) {
                for (int j = 0; j < m2->width; j++) {
                    MAT_CELL(s, i, j) += MAT_CELL(m1, k, i) * MAT_CELL(m2, k, j);
                }
            }
        }
    }

    return s;
}

/*
 * Perform matrix multiplication between an AxB matrix and a CxB matrix, without transposing the second matrix
 */
Matrix *mat_mul_t2(Matrix* m1, Matrix* m2) {
    assert(m1->width == m2->width, "The width of the two matrices should be equal.");
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s = new_mat(m1->height, m2->height, m1->on_cluster);

    if (m1->on_cluster) {
        // Pack operands and result in a struct
        ClOpMMM arg = {m1, m2, s};
        // Perform GEMM in parallel
        pi_cl_team_fork(CORES, mat_mul_t2_par, (void *) &arg);
    } else {
        // i-j-k loop to optimize memory access
        for (int i = 0; i < m1->height; i++) {
            for (int j = 0; j < m2->height; j++) {
                for (int k = 0; k < m1->width; k++) {
                    MAT_CELL(s, i, j) += MAT_CELL(m1, i, k) * MAT_CELL(m2, j, k);
                }
            }
        }
    }

    return s;
}

/*
 * Perform dot product between two vectors
 */
fp dot(Matrix *v1, Matrix *v2) {
    assert(is_vector(v1) && is_vector(v2), "The two matrices should be vectors (i.e. either their height or width should be equal to 1).");
    assert(v1->on_cluster == v2->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s;

    if (is_row_vector(v1)) {  // V1: (1, N)
        if (is_row_vector(v2)) {  // V2: (1, N)
            assert(v1->width == v2->width, "The two vectors should have the same length.");
            s = mat_mul_t2(v1, v2);  // (1, N) @ (1, N).T = (1, N) @ (N, 1) = (1, 1)
        } else {  // V2: (N, 1)
            assert(v1->width == v2->height, "The two vectors should have the same length.");
            s = mat_mul(v1, v2);  // (1, N) @ (N, 1) = (1, 1)
        }
    } else {  // V1: (N, 1)
        if (is_row_vector(v2)) {  // V2: (1, N)
            assert(v1->height == v2->width, "The two vectors should have the same length.");
            s = mat_mul(v2, v1);  // (1, N) @ (N, 1) = (1, 1)
        } else {  // V2: (N, 1)
            assert(v1->height == v2->height, "The two vectors should have the same length.");
            s = mat_mul_t1(v1, v2);  // (N, 1).T @ (N, 1) = (1, N) @ (N, 1) = (1, 1)
        }
    }

    fp res = MAT_CELL(s, 0, 0);
    free_mat(s);
    return res;
}

/*
 * Perform outer product between two vectors
 */
Matrix *outer(Matrix *v1, Matrix *v2) {
    assert(is_vector(v1) && is_vector(v2), "The two matrices should be vectors (i.e. either their height or width should be equal to 1).");
    assert(v1->on_cluster == v2->on_cluster, "The two matrices should be allocated on the same memory.");
    Matrix *s;

    if (is_row_vector(v1)) {  // V1: (1, N)
        if (is_row_vector(v2)) {  // V2: (1, M)
            s = mat_mul_t1(v1, v2);  // (1, N).T @ (1, M) = (N, 1) @ (1, M) = (N, M)
        } else {  // V2: (M, 1)
            Matrix *v2_t = transpose(v2);  // (M, 1).T = (1, M)
            s = mat_mul_t1(v1, v2_t);  // (1, N).T @ (1, M) = (N, 1) @ (1, M) = (N, M)
            free_mat(v2_t);
        }
    } else {  // V1: (N, 1)
        if (is_row_vector(v2)) {  // V2: (1, M)
            s = mat_mul(v1, v2);  // (N, 1) @ (1, M) = (N, M)
        } else {  // V2: (M, 1)
            s = mat_mul_t2(v1, v2);  // (N, 1) @ (M, 1).T = (N, 1) @ (1, M) = (N, M)
        }
    }

    return s;
}

/*
 * Check if two matrices are equal
 */
bool are_equal(Matrix *m1, Matrix *m2, fp tol) {
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
bool is_square(Matrix *m) {
    return m->height == m->width;
}

/*
 * Check if a matrix is a vector
 */
bool is_vector(Matrix *m) {
    return m->height == 1 || m->width == 1;
}

/*
 * Check if a matrix is a row vector
 */
bool is_row_vector(Matrix *m) {
    return m->height == 1;
}

/*
 * Check if a matrix is a column vector
 */
bool is_col_vector(Matrix *m) {
    return m->width == 1;
}

/*
 * Copy the values on the diagonal of a matrix into a new column vector
 */
Matrix *diagonal(Matrix *m) {
    // Get minimum between height and width
    int dim = (m->height > m->width ? m->width : m->height);
    Matrix *d = new_vec(dim, m->on_cluster);

    for (int i = 0; i < dim; i++) {
        MAT_CELL(d, i, 0) = MAT_CELL(m, i, i);
    }

    return d;
}

/*
 * Set cells below the main diagonal to zero (in-place operation, it modifies the matrix)
 */
void tri_up(Matrix *m) {
    for (int i = 0; i < m->height; i++) {
        for (int j = 0; j < i; j++) {
            MAT_CELL(m, i, j) = 0;
        }
    }
}

/*
 * Slice given matrix and return a new matrix
 */
Matrix *read_slice(Matrix *m, uint row_start, uint row_stop, uint col_start, uint col_stop) {
    // Input check
    assert(0 <= row_start && row_start < m->height, "The row_start argument is out of bounds.");
    assert(0 <= row_stop && row_stop < m->height, "The row_stop argument is out of bounds.");
    assert(0 <= col_start && col_start < m->width, "The col_start argument is out of bounds.");
    assert(0 <= col_stop && col_stop < m->width, "The col_stop argument is out of bounds.");
    assert(row_start <= row_stop, "The row_stop argument must be greater than or equal to row_start.");
    assert(col_start <= col_stop, "The col_stop argument must be greater than or equal to col_start.");

    Matrix *s = new_mat(row_stop - row_start + 1, col_stop - col_start + 1, m->on_cluster);

    for (int i = row_start; i <= row_stop; i++) {
        for (int j = col_start; j <= col_stop; j++) {
            MAT_CELL(s, i - row_start, j - col_start) = MAT_CELL(m, i, j);
        }
    }

    return s;
}

/*
 * Write into sliced matrix, modifying it
 */
void write_slice(Matrix *m1, Matrix *m2, uint row_start, uint col_start) {
    int row_stop = row_start + m2->height - 1;
    int col_stop = col_start + m2->width - 1;
    // Input check
    assert(0 <= row_start && row_start < m1->height, "The row_start argument is out of bounds.");
    assert(0 <= row_stop && row_stop < m1->height, "The row_stop argument is out of bounds.");
    assert(0 <= col_start && col_start < m1->width, "The col_start argument is out of bounds.");
    assert(0 <= col_stop && col_stop < m1->width, "The col_stop argument is out of bounds.");

    for (int i = row_start; i <= row_stop; i++) {
        for (int j = col_start; j <= col_stop; j++) {
            MAT_CELL(m1, i, j) = MAT_CELL(m2, i - row_start, j - col_start);
        }
    }
}

/*
 * Extract the k-th row, creating a new vector
 */
Matrix *extract_row(Matrix *m, uint k) {
    assert(0 <= k && k < m->height, "Index is out of bounds for rows.");
    Matrix *r = new_mat(1, m->width, m->on_cluster);

    for (int i = 0; i < m->width; i++) {
        MAT_CELL(r, 0, i) = MAT_CELL(m, k, i);
    }

    return r;
}

/*
 * Extract the k-th column, creating a new vector
 */
Matrix *extract_col(Matrix *m, uint k) {
    assert(0 <= k && k < m->width, "Index is out of bounds for columns.");
    Matrix *c = new_vec(m->height, m->on_cluster);

    for (int i = 0; i < m->height; i++) {
        MAT_CELL(c, i, 0) = MAT_CELL(m, i, k);
    }

    return c;
}

/*
 * Copy the values from a vector into the specified matrix row (it modifies the matrix)
 */
void paste_row(Matrix *m, Matrix *r, uint k) {
    assert(0 <= k && k < m->height, "Index is out of bounds for rows.");
    assert(is_row_vector(r), "The second matrix should be a row vector (i.e. with height equal to 1).");
    assert(m->width == r->width, "The width of the first matrix and the length of the row vector (i.e. its width) should be the same.");

    for (int i = 0; i < m->width; i++) {
        MAT_CELL(m, k, i) = MAT_CELL(r, 0, i);
    }
}

/*
 * Copy the values from a vector into the specified matrix column (it modifies the matrix)
 */
void paste_col(Matrix *m, Matrix *c, uint k) {
    assert(0 <= k && k < m->width, "Index is out of bounds for columns.");
    assert(is_col_vector(c), "The second matrix should be a column vector (i.e. with width equal to 1).");
    assert(m->height == c->height, "The height of the first matrix and the length of the column vector (i.e. its height) should be the same.");

    for (int i = 0; i < m->height; i++) {
        MAT_CELL(m, i, k) = MAT_CELL(c, i, 0);
    }
}

/*
 * Prints a matrix to standard output
 */
void print_mat(Matrix *m) {
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
void write_mat(const char *path, Matrix *m) {
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
Matrix *copy_mat(Matrix *m) {
    Matrix *s = new_mat(m->height, m->width, m->on_cluster);

    for (int i = 0; i < m->height; i++) {
        for (int j = 0; j < m->width; j++) {
            MAT_CELL(s, i, j) = MAT_CELL(m, i, j);
        }
    }

    return s;
}