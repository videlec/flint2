/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void
_fmpq_poly_cosh_series(fmpz * g, fmpz_t gden, 
                       const fmpz * h, const fmpz_t hden, slong hlen, slong n)
{
    fmpz * t;
    fmpz_t tden;

    t = _fmpz_vec_init(n);
    fmpz_init(tden);

    /* cosh(x) = (exp(x)+exp(-x))/2 */
    _fmpq_poly_exp_series(g, gden, h, hden, hlen, n);
    _fmpq_poly_inv_series(t, tden, g, gden, n, n);
    _fmpq_poly_add(g, gden, g, gden, n, t, tden, n);
    _fmpq_poly_scalar_div_ui(g, gden, g, gden, n, UWORD(2));
    _fmpq_poly_canonicalise(g, gden, n);

    _fmpz_vec_clear(t, n);
    fmpz_clear(tden);
}

void fmpq_poly_cosh_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
{
    if (n == 0)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (poly->length == 0 || n == 1)
    {
        fmpq_poly_one(res);
        return;
    }

    if (!fmpz_is_zero(poly->coeffs))
    {
        flint_printf("Exception (fmpq_poly_cosh_series). Constant term != 0.\n");
        flint_abort();
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_cosh_series(res->coeffs, res->den,
            poly->coeffs, poly->den, poly->length, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_cosh_series(t->coeffs, t->den,
            poly->coeffs, poly->den, poly->length, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}
