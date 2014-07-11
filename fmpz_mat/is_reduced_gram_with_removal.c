/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "fmpz_mat.h"
#include "fmpq_mat.h"
#include "fmpz_lll.h"

int
fmpz_mat_is_reduced_gram_with_removal(const fmpz_mat_t A, double delta,
                                      double eta, const fmpz_t gs_B, int newd)
{
    slong i, j, k, d = A->r;
    fmpz_mat_t Atmp;
    fmpz_lll_t flfp;
    fmpq_mat_t r, mu;
    fmpq *s;
    fmpq_lll_t fl;
    fmpq_t tmp, gs_Bq;

    if (d == 0 || d == 1)
        return 1;

    fmpq_mat_init(r, d, d);
    fmpq_mat_init(mu, d, d);

    s = _fmpq_vec_init(d);

    fmpq_lll_context_init(fl, delta, eta);

    fmpq_init(tmp);
    fmpq_init(gs_Bq);

    fmpz_set(fmpq_numref(gs_Bq), gs_B);
    fmpz_one(fmpq_denref(gs_Bq));

    fmpz_set(fmpq_mat_entry_num(r, 0, 0), fmpz_mat_entry(A, 0, 0));
    if (newd == 0 && fmpz_cmp(fmpz_mat_entry(A, 0, 0), gs_B) < 0)
    {
        fmpq_mat_clear(r);
        fmpq_mat_clear(mu);
        fmpq_lll_context_clear(fl);
        fmpq_clear(tmp);
        fmpq_clear(gs_Bq);
        _fmpq_vec_clear(s, d);
        return 0;
    }

    for (i = 1; i < newd; i++)
    {
        fmpz_set(fmpq_numref(s), fmpz_mat_entry(A, i, i));
        fmpz_one(fmpq_denref(s));
        for (j = 0; j <= i - 1; j++)
        {
            fmpz_set(fmpq_mat_entry_num(r, i, j), fmpz_mat_entry(A, i, j));
            for (k = 0; k <= j - 1; k++)
            {
                fmpq_submul(fmpq_mat_entry(r, i, j), fmpq_mat_entry(mu, j, k),
                            fmpq_mat_entry(r, i, k));
            }
            fmpq_div(fmpq_mat_entry(mu, i, j), fmpq_mat_entry(r, i, j),
                     fmpq_mat_entry(r, j, j));
            fmpq_abs(tmp, fmpq_mat_entry(mu, i, j));
            if (fmpq_cmp(tmp, fl->eta) > 0) /* check size reduction */
            {
                fmpq_mat_clear(r);
                fmpq_mat_clear(mu);
                fmpq_lll_context_clear(fl);
                fmpq_clear(tmp);
                fmpq_clear(gs_Bq);
                _fmpq_vec_clear(s, d);
                return 0;
            }
            fmpq_set(s + j + 1, s + j);
            fmpq_submul(s + j + 1, fmpq_mat_entry(mu, i, j),
                        fmpq_mat_entry(r, i, j));
        }
        fmpq_set(fmpq_mat_entry(r, i, i), s + i);
        if (i >= newd && fmpq_cmp(fmpq_mat_entry(r, i, i), gs_Bq) < 0)  /* check removals */
        {
            fmpq_mat_clear(r);
            fmpq_mat_clear(mu);
            fmpq_lll_context_clear(fl);
            fmpq_clear(tmp);
            fmpq_clear(gs_Bq);
            _fmpq_vec_clear(s, d);
            return 0;
        }
        fmpq_mul(tmp, fl->delta, fmpq_mat_entry(r, i - 1, i - 1));
        if (fmpq_cmp(tmp, s + i - 1) > 0)   /* check Lovasz condition */
        {
            fmpq_mat_clear(r);
            fmpq_mat_clear(mu);
            fmpq_lll_context_clear(fl);
            fmpq_clear(tmp);
            fmpq_clear(gs_Bq);
            _fmpq_vec_clear(s, d);
            return 0;
        }
    }

    fmpz_mat_init_set(Atmp, A);
    fmpz_lll_context_init(flfp, delta, eta, GRAM, EXACT, 0);
    fmpz_lll_d(Atmp, NULL, flfp);

    fmpz_set(fmpq_mat_entry_num(r, 0, 0), fmpz_mat_entry(Atmp, 0, 0));
    if (newd == 0 && fmpz_cmp(fmpz_mat_entry(Atmp, 0, 0), gs_B) < 0)
    {
        fmpz_mat_clear(Atmp);
        fmpq_mat_clear(r);
        fmpq_mat_clear(mu);
        fmpq_lll_context_clear(fl);
        fmpq_clear(tmp);
        fmpq_clear(gs_Bq);
        _fmpq_vec_clear(s, d);
        return 0;
    }

    for (i = 1; i < d; i++)
    {
        fmpz_set(fmpq_numref(s), fmpz_mat_entry(Atmp, i, i));
        fmpz_one(fmpq_denref(s));
        for (j = 0; j <= i - 1; j++)
        {
            fmpz_set(fmpq_mat_entry_num(r, i, j), fmpz_mat_entry(Atmp, i, j));
            fmpz_one(fmpq_mat_entry_den(r, i, j));
            for (k = 0; k <= j - 1; k++)
            {
                fmpq_submul(fmpq_mat_entry(r, i, j), fmpq_mat_entry(mu, j, k),
                            fmpq_mat_entry(r, i, k));
            }
            fmpq_div(fmpq_mat_entry(mu, i, j), fmpq_mat_entry(r, i, j),
                     fmpq_mat_entry(r, j, j));
            fmpq_set(s + j + 1, s + j);
            fmpq_submul(s + j + 1, fmpq_mat_entry(mu, i, j),
                        fmpq_mat_entry(r, i, j));
        }
        fmpq_set(fmpq_mat_entry(r, i, i), s + i);
        if (i >= newd && fmpq_cmp(fmpq_mat_entry(r, i, i), gs_Bq) < 0)  /* check removals */
        {
            fmpz_mat_clear(Atmp);
            fmpq_mat_clear(r);
            fmpq_mat_clear(mu);
            fmpq_lll_context_clear(fl);
            fmpq_clear(tmp);
            fmpq_clear(gs_Bq);
            _fmpq_vec_clear(s, d);
            return 0;
        }
    }
    fmpz_mat_clear(Atmp);

    fmpq_mat_clear(r);
    fmpq_mat_clear(mu);
    fmpq_lll_context_clear(fl);
    fmpq_clear(tmp);
    fmpq_clear(gs_Bq);
    _fmpq_vec_clear(s, d);
    return 1;
}
