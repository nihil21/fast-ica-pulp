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

#include "../include/sorting.h"

/*
 * Swap two elements
 */
void swap(fp *a, fp *b){
    fp tmp = *a;
    *a = *b;
    *b = tmp;
}

/*
 * Swap two indexes
 */
void swap_i(int *a, int *b){
    fp tmp = *a;
    *a = *b;
    *b = tmp;
}

/*
 * Compare two elements
 */
int compare(fp a, fp b){
    if (a > b)
        return 1;
    else if (a == b)
        return 0;

    return -1;
}

/*
 * Recursive sub-routine of QuickSort
 */
void quick_sort_rec(fp v[], int *sort_idx, int start, int end, int desc_fact){
    int i, j, i_pivot;
    fp pivot;

    if(start < end) {
        i = start;
        j = end;
        i_pivot = end;
        pivot = v[i_pivot];

        do {
            while (i < j && desc_fact * compare(v[i], pivot) <= 0)
                i++;
            while (j > i && desc_fact * compare(v[j], pivot) >= 0)
                j--;
            if (i < j) {
                swap(&(v[i]), &(v[j]));
                swap_i(&(sort_idx[i]), &(sort_idx[j]));
            }
        } while (i < j);

        if (i != i_pivot && desc_fact * compare(v[i], v[i_pivot])) {
            swap(&(v[i]), &(v[i_pivot]));
            swap_i(&(sort_idx[i]), &(sort_idx[i_pivot]));
            i_pivot = i;
        }

        if (start < i_pivot - 1)
            quick_sort_rec(v, sort_idx, start, i_pivot - 1, desc_fact);
        if (i_pivot + 1 < end)
            quick_sort_rec(v, sort_idx, i_pivot + 1, end, desc_fact);
    }
}

/*
 * QuickSort implementation
 */
void quick_sort(fp v[], int sort_idx[], int len, bool desc) {
    int desc_fact = desc ? -1 : 1;
    
    // Create array of indexes
    for (int i = 0; i < len; i++)
        sort_idx[i] = i;

    quick_sort_rec(v, sort_idx, 0, len - 1, desc_fact);
}
