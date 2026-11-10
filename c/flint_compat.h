/*
 * Compatibility with FLINT < 3.1, which has no _fmpz_vec_dot_general (used by
 * root-unitary's power_sums.c).  The Makefile includes this header before
 * weilpoly.c.  The function below follows the FLINT 3.1 documentation:
 *
 *     res = initial + (-1)^subtract * sum_{i < len} a[i] b[i]
 *
 * with initial = 0 when it is NULL, and b reversed when reverse is set.
 * Define K3FQ_FLINT_COMPAT to use it with newer FLINT as well (for testing).
 */
#ifndef K3FQ_FLINT_COMPAT_H
#define K3FQ_FLINT_COMPAT_H

#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>

#if __FLINT_RELEASE < 30100 || defined(K3FQ_FLINT_COMPAT)
static void k3fq_fmpz_vec_dot_general(fmpz_t res, const fmpz_t initial, int subtract,
                                      const fmpz * a, const fmpz * b, int reverse, slong len)
{
    fmpz_t s;
    slong i;
    fmpz_init(s);
    for (i = 0; i < len; i++)
        fmpz_addmul(s, a + i, reverse ? b + len - 1 - i : b + i);
    if (subtract)
        fmpz_neg(s, s);
    if (initial != NULL)
        fmpz_add(s, s, initial);
    fmpz_swap(res, s);
    fmpz_clear(s);
}
#define _fmpz_vec_dot_general k3fq_fmpz_vec_dot_general
#endif

#endif
