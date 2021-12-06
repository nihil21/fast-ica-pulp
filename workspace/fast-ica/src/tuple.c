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
    assert(m1->on_cluster == m2->on_cluster, "The two matrices should be allocated on the same memory.");
    Tuple *tuple;
    
    if (m1->on_cluster) {
        tuple = pi_l1_malloc(&cluster_dev, sizeof(Tuple));
    } else {
        tuple = pi_l2_malloc(sizeof(Tuple));
    }
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
    bool on_cluster = tuple->m1->on_cluster;

    if (tuple != NULL) {
        // Optionally free matrix members
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
        if (on_cluster) {
            pi_l1_free(&cluster_dev, tuple, sizeof(Tuple));
        } else {
            pi_l2_free(tuple, sizeof(Tuple));
        }
        tuple = NULL;
    }
}