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

#include "../include/linalg.h"
#include "../include/stat_alloc.h"
#include "../include/random.h"
#include "../include/cluster.h"

/*
 * Compute a Householder reflection for a given real vector
 */
void generate_householder(const Matrix *s, const Matrix *m, fp *tau) {
    fp m0 = MAT_CELL(m, 0, 0);
    scale_par(s, m, 1 / (m0 + (fp) SGN(m0) * norm_par(m)));
    MAT_CELL(s, 0, 0) = 1;
    *tau = 2 / dot_par(s, s);
}

/*
 * Reduce the given real matrix to Hessenberg form by means of Householder reflections
 */
void to_hessenberg(const Matrix *h, const Matrix *m) {
    int n = m->height;
    // Copy matrix M into H
    copy_mat_par(h, m);

    // Householder vector
    Matrix v_k = {
        .data = prep_tmp1_data_l1,
        .width = 1,
        .offset = 1,
        .on_cluster = true
    };
    // Tmp matrices
    Matrix tmp1 = {
        .data = prep_tmp2_data_l1,
        .on_cluster = true
    };
    Matrix tmp2 = {
        .data = prep_tmp3_data_l1,
        .on_cluster = true
    };

    for (int k = 0; k < n - 2; k++) {
        // Extract k-th column from H
        Matrix h_k = slice(h, k + 1, n - 1, k, k);  // column vector

        // Generate k-th Householder reflection, which can be applied as follows: PX = X - tau * v @ v.T @ X
        fp tau;
        v_k.height = n - k - 1;
        generate_householder(&v_k, &h_k, &tau);

        // Extract (n - k - 1)x(n - k) submatrix
        h_k = slice(h, k + 1, n - 1, k, n - 1);
        
        // Apply Householder from the left:
        // H_k = H_k - tau * v_k @ v_k.T @ H_k
        tmp1.height = 1;  // row vector
        tmp1.width = h_k.width;
        tmp1.offset = h_k.width;
        tmp2.height = v_k.height;
        tmp2.width = tmp1.width;
        tmp2.offset = tmp1.width;
        mat_mul_t1_par(&tmp1, &v_k, &h_k);
        outer_par(&tmp2, &v_k, &tmp1);
        scale_par_(&tmp2, tau);
        sub_mat_par_(&h_k, &tmp2);

        // Extract (n)x(n - k - 1) submatrix
        h_k = slice(h, 0, n - 1, k + 1, n - 1);

        // Apply Householder from the right:
        // H_k = H_k - tau * (H_k @ v_k) @ v_k.T
        tmp1.height = h_k.height;
        tmp1.width = 1;  // column vector
        tmp1.offset = 1;
        mat_mul_par(&tmp1, &h_k, &v_k);
        tmp2.height = tmp1.height;
        tmp2.width = v_k.height;
        tmp2.offset = v_k.height;
        mat_mul_t2_par(&tmp2, &tmp1, &v_k);
        scale_par_(&tmp2, tau);
        sub_mat_par_(&h_k, &tmp2);
    }
}

/*
 * QR decomposition using Householder reflections
 */
void qr_decomposition(const Matrix *q, const Matrix *r, const Matrix *m) {
    int n = m->height;
    // Copy matrix M into R
    copy_mat_par(r, m);
    // Initialize Q to identity
    eye(q);

    // Householder vector
    Matrix v_k = {
        .data = prep_tmp1_data_l1,
        .width = 1,
        .offset = 1,
        .on_cluster = true
    };
    // H matrix
    Matrix h = {
        .data = prep_tmp2_data_l1,
        .height = n,
        .width = n,
        .offset = n,
        .on_cluster = true
    };
    // Tmp matrices
    Matrix tmp = {
        .data = prep_tmp3_data_l1,
        .on_cluster = true
    };

    for (int k = 0; k < n; k++) {
        // Extract k-th column
        Matrix r_k = slice(r, k, n - 1, k, k);

        // Generate Householder reflection
        fp tau;
        v_k.height = n - k;
        generate_householder(&v_k, &r_k, &tau);

        tmp.height = v_k.height;
        tmp.width = v_k.height;
        tmp.offset = v_k.height;

        // Initialize H to identity
        eye(&h);
        // Compute H_k
        Matrix h_k = slice(&h, k, n - 1, k, n - 1);
        outer_par(&tmp, &v_k, &v_k);
        scale_par_(&tmp, tau);
        sub_mat_par_(&h_k, &tmp);

        tmp.height = n;
        tmp.width = n;
        tmp.offset = n;

        // Q(k + 1) = Q(k) @ H(k)
        mat_mul_par(&tmp, q, &h);
        copy_mat_par(q, &tmp);

        // R(k + 1) = H(k) @ R(k)
        mat_mul_par(&tmp, &h, r);
        copy_mat_par(r, &tmp);
    }
    
    // Make R upper triangular
    tri_up_par(r);
}

/*
 * Solve a linear system Ux = y, where U is upper triangular, by means of back-substitution
 */
void back_substitution(const Matrix *x, const Matrix *u, const Matrix *y) {
    // NB: both x and y are column vectors
    int n_eq = y->height;

    // Iterate over the cells of y backwards
    for (int i = n_eq - 1; i >= 0; i--) {
        fp back_substitute = 0;
        // Iterate over the subsequent cells and accumulate back substitutions
        for (int j = i + 1; j < n_eq; j++) {
            back_substitute += MAT_CELL(x, j, 0) * MAT_CELL(u, i, j);
        }
        // Compute i-th solution
        MAT_CELL(x, i, 0) = (MAT_CELL(y, i, 0) - back_substitute) / MAT_CELL(u, i, i);
    }
}

/*
 * Solve a linear system Ax = b by means of QR decomposition:
 * A = QR => QRx = b and, thanks to the orthogonality of Q, Rx = Q.Tb => Rx = y
 */
void lin_solve(const Matrix *x, const Matrix *a, const Matrix *b) {
    uint n = a->height;
    // Decompose matrix A using QR
    Matrix q = {
        .data = prep_tmp6_data_l1,
        .height = n,
        .width = n,
        .offset = n,
        .on_cluster = true
    };
    Matrix r = {
        .data = prep_tmp7_data_l1,
        .height = n,
        .width = n,
        .offset = n,
        .on_cluster = true
    };
    qr_decomposition(&q, &r, a);

    // Multiply the transpose of Q and b
    Matrix y = {
        .data = prep_tmp1_data_l1,
        .height = n,
        .width = 1,
        .offset = 1,
        .on_cluster = true
    };
    mat_mul_t1_par(&y, &q, b);
    // Solve Rx = y by means of back-substitution
    back_substitution(x, &r, &y);
}

/*
 * Compute the eigenvector associated to an eigenvalue by means of inverse iteration
 */
void inv_iter(const Matrix* eig_vec, fp eig_val, const Matrix *m, fp tol, int max_iter) {
    uint n = m->height;
    // Perturb lambda to prevent the computed matrix from becoming singular
    fp lambda = eig_val + (fp) uniform_rand() * 1e-6f;
    // Compute M' = M - lambda * I
    Matrix lambda_i = {
        .data = prep_tmp3_data_l1,
        .height = n,
        .width = n,
        .offset = n,
        .on_cluster = true
    };
    eye(&lambda_i);
    scale_par_(&lambda_i, lambda);
    Matrix m_prime = {
        .data = prep_tmp4_data_l1,
        .height = n,
        .width = n,
        .offset = n,
        .on_cluster = true
    };
    sub_mat_par(&m_prime, m, &lambda_i);

    // Initialize vector randomly
    mat_randn(eig_vec);
    Matrix prev = {
        .data = prep_tmp5_data_l1,
        .height = n,
        .width = 1,
        .offset = 1,
        .on_cluster = true
    };

    uint i = 0;
    do {
        // Save previous estimate and compute the new one
        copy_mat_par(&prev, eig_vec);
        lin_solve(eig_vec, &m_prime, &prev);

        // If the first entry of the current estimate is negative,
        // swap the sign to improve convergence
        if (MAT_CELL(eig_vec, 0, 0) < 0) {
            scale_par_(eig_vec, -1);
        }

        // Normalize estimate
        fp v_norm = norm_par(eig_vec);
        scale_par_(eig_vec, 1 / v_norm);

        i++;
    } while ((!are_equal(eig_vec, &prev, tol)) && i < max_iter);
}

/*
 * Compute the eigenvalues of a square matrix by means of QR decomposition with shift and Hessenberg reduction
 */
void solve_eig_vals(const Matrix *eig_vals, const Matrix *m, fp tol, int max_iter) {
    uint n = m->height;

    // Tri-diagonalize matrix using Hessenberg
    Matrix t = {
        .data = prep_tmp4_data_l1,
        .height = n,
        .width = n,
        .offset = n,
        .on_cluster = true
    };
    to_hessenberg(&t, m);

    int k = n - 1;
    uint i = 0;
    while (k > 0 && i < max_iter) {
        // Obtain the shift from the lower right corner of the matrix.
        Matrix mu = {
            .data = prep_tmp5_data_l1,
            .height = k + 1,
            .width = k + 1,
            .offset = k + 1,
            .on_cluster = true
        };
        eye(&mu);
        scale_par_(&mu, MAT_CELL(&t, k, k));

        // Shift T matrix and perform QR on shifted matrix
        sub_mat_par_(&t, &mu);
        Matrix q = {
            .data = prep_tmp6_data_l1,
            .height = t.height,
            .width = t.width,
            .offset = t.width,
            .on_cluster = true
        };
        Matrix r = {
            .data = prep_tmp7_data_l1,
            .height = t.height,
            .width = t.width,
            .offset = t.width,
            .on_cluster = true
        };
        qr_decomposition(&q, &r, &t);

        // Multiply R*Q and shift back result
        mat_mul_par(&t, &r, &q);
        add_mat_par_(&t, &mu);

        if (ABS(MAT_CELL(&t, k, k - 1)) < tol) {
            MAT_CELL(eig_vals, k, 0) = MAT_CELL(&t, k, k);
            t = slice(&t, 0, k - 1, 0, k - 1);
            compact(&t, 0, 0, n * n);
            k--;
        }
        i++;
    }
    MAT_CELL(eig_vals, 0, 0) = MAT_CELL(&t, 0, 0);
}

/*
 * Compute the eigenvectors of a square matrix given the eigenvalues
 */
void solve_eig_vecs(const Matrix *eig_vecs, const Matrix *eig_vals, const Matrix *m, fp tol, int max_iter) {
    // Iterate over eigenvalues
    for (int i = 0; i < eig_vals->height; i++) {
        // Extract current eigenvalue
        fp eig_val = MAT_CELL(eig_vals, i, 0);  // eig_vals is a column vector
        // Compute i-th eigenvector (column-wise)
        Matrix eig_vec = extract_col(eig_vecs, i);
        inv_iter(&eig_vec, eig_val, m, tol, max_iter);
    }
}

/*
 * Compute the eigenvalues and eigenvectors of a square matrix by means of QR algorithm
 */
void solve_eig(const Matrix *eig_vals, const Matrix *eig_vecs, const Matrix *m) {
    fp tol = 1e-7;
    solve_eig_vals(eig_vals, m, tol, MAX_ITER);
    solve_eig_vecs(eig_vecs, eig_vals, m, tol, MAX_ITER);
}

/*
 * Given a (n_features, n_samples) matrix, compute the (n_features, n_features) covariance matrix
 */
void covariance(const Matrix *cov_mat, const Matrix *mean_vec, const Matrix *x) {
    const uint n_samples = x->width;

    // 1. Center data
    col_mean_par(mean_vec, x);
    sub_col_par_(x, mean_vec);

    // 2. Perform mat_mul with itself
    mat_mul_t2_par(cov_mat, x, x);
    
    // 3. Normalize
    fp fact = 1 / ((fp) n_samples - 1);
    scale_par_(cov_mat, fact);
}
