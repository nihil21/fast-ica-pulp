//
// Created by nihil on 04/10/21.
//

#include "pmsis.h"
#include "../include/tuple.h"
#include "../include/utils.h"

extern struct pi_device cluster_dev;

/*
 * Creates a tuple from a pair of matrices (Cluster's L1 memory)
 */
Tuple *new_tuple(Matrix *m1, Matrix *m2) {
    Tuple *tuple;
    tuple = pi_l1_malloc(&cluster_dev, sizeof(Tuple));
    assert(tuple != NULL, "Could not allocate tuple.");

    // Set fields
    tuple->m1 = m1;
    tuple->m2 = m2;

    return tuple;
}

/*
 * Free the memory allocated for a given tuple (Cluster's L1 memory)
 */
void free_tuple(Tuple *tuple, bool free_members) {
    if (tuple != NULL) {
        // Optionally free Matrix members
        if (free_members) {
            // Free first matrix
            if (tuple->m1 != NULL) {
                free_mat(tuple->m1);
            }
            // Free second matrix
            if (tuple->m2 != NULL) {
                free_mat(tuple->m2);
            }
        }

        // Free tuple
        pi_l1_free(&cluster_dev, tuple, sizeof(Tuple));
        tuple = NULL;
    }
}