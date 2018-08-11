#include <stdio.h>
#include <stdlib.h>

static int cmp_int_ascending(const void *x, const void *y) 
{
    const int xx = *(const int *) x;
    const int yy = *(const int *) y;
    if (xx < yy) return -1;
    if (xx > yy) return  1;
    return 0;
}
/*!
 * @brief Searches a sorted array, values, for the key.
 * @note This has been modified from ISTI's ISCL and made specific for
 *       arrays sorted in increasing order.
 * @param[in] key     Item to search for in values.
 * @param[in] values  Array sorted in ascending order.  This has dimension [n].
 * @param[in] n       Number of elements in values.
 * @param[out] indx   On successful exit this is the index in vaues that
 *                    matches the key.
 * @param[out] ierr   0 indicates success.
 * @author Ben Baker
 * @copyright ISTI distributed under the Apache 2 license.
 */
void xcloc_sort_bsearch32i(const int key, const int values[],
                           const int n, int *indx, int *ierr)
{
    int *item = NULL;
    size_t np = (size_t) n;
    *ierr = 0;
    *indx = 0;
    if (np < 1)
    {
        *ierr = 1;
        return;
    }
    // out of bounds and edge case
    if (key < values[0])
    {
        *ierr = 1;
        return;
    }
    if (key > values[n-1])
    {
        *ierr = 1;
        return;
    }
    // look for it
    item = (int *) bsearch((const void *) &key, (const void *) values,
                            np, sizeof(int),
                            cmp_int_ascending);
    if (item == NULL)
    {
        *ierr = 1;
        *indx =-1;
    }
    else
    {
        *indx = item - values;
    }
    return;
}
