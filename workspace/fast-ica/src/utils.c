//
// Created by nihil on 21/09/21.
//

#include <stdio.h>
#include "pmsis.h"
#include "../include/utils.h"

/*
 * Verify if condition holds, otherwise terminate program
 */
void assert(bool condition, char *message) {
    if (!condition) {
        printf("%s\n", message);
        pmsis_exit(-1);
    }
}