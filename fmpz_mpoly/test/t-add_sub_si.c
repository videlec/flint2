/*
    Copyright (C) 2017 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("add/sub_si....");
    fflush(stdout);

    /* Check (f + a) - a = f */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g, h;
       ordering_t ord;
       slong c;
       slong nvars, len, exp_bound, coeff_bits, exp_bits;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 20) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);
       fmpz_mpoly_init(h, ctx);

       len = n_randint(state, 100);

       exp_bits = n_randint(state, FLINT_BITS -
                     mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars) - 1) + 1;
       exp_bound = n_randbits(state, exp_bits);
       coeff_bits = n_randint(state, 200);

       fmpz_mpoly_randtest(f, state, len, exp_bound, coeff_bits, ctx);
       fmpz_mpoly_randtest(g, state, len, exp_bound, coeff_bits, ctx);
       fmpz_mpoly_randtest(h, state, len, exp_bound, coeff_bits, ctx);

       for (j = 0; j < 10; j++)
       {
          len = f->length;

          c = z_randtest(state);

          fmpz_mpoly_add_si(g, f, c, ctx);

          fmpz_mpoly_sub_si(h, g, c, ctx);

          result = fmpz_mpoly_equal(f, h, ctx);

          if (!result)
          {
             const char * vars[20];
             vars[0] = "x1", vars[1] = "x2", vars[2] = "x3", vars[3] = "x4",
             vars[4] = "x5", vars[5] = "x6", vars[6] = "x7", vars[7] = "x8",
             vars[8] = "x9", vars[9] = "x10", vars[10] = "x11", vars[11] = "x12",
             vars[12] = "x13", vars[13] = "x14", vars[14] = "x15", vars[15] = "x16",
             vars[16] = "x17", vars[17] = "x18", vars[18] = "x19", vars[19] = "x20";

             printf("FAIL\n");

             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                                  len, exp_bits, exp_bound, coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(h, vars, ctx); printf("\n\n");
          
             flint_printf("c = %wd\n", c);

             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
       fmpz_mpoly_clear(h, ctx);  
    }

    /* Check aliasing */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
       fmpz_mpoly_ctx_t ctx;
       fmpz_mpoly_t f, g;
       ordering_t ord;
       slong c;
       slong nvars, len, exp_bound, coeff_bits, exp_bits;

       ord = mpoly_ordering_randtest(state);
       nvars = n_randint(state, 20) + 1;

       fmpz_mpoly_ctx_init(ctx, nvars, ord);

       fmpz_mpoly_init(f, ctx);
       fmpz_mpoly_init(g, ctx);

       len = n_randint(state, 100);

       exp_bits = n_randint(state, FLINT_BITS -
                     mpoly_ordering_isdeg(ord)*FLINT_BIT_COUNT(nvars) - 1) + 1;
       exp_bound = n_randbits(state, exp_bits);
       coeff_bits = n_randint(state, 200);

       fmpz_mpoly_randtest(f, state, len, exp_bound, coeff_bits, ctx);
       fmpz_mpoly_set(g, f, ctx);

       for (j = 0; j < 10; j++)
       {
          len = f->length;

          c = z_randtest(state);

          fmpz_mpoly_add_si(f, f, c, ctx);

          fmpz_mpoly_sub_si(f, f, c, ctx);

          result = fmpz_mpoly_equal(f, g, ctx);

          if (!result)
          {
             const char * vars[20];
             vars[0] = "x1", vars[1] = "x2", vars[2] = "x3", vars[3] = "x4",
             vars[4] = "x5", vars[5] = "x6", vars[6] = "x7", vars[7] = "x8",
             vars[8] = "x9", vars[9] = "x10", vars[10] = "x11", vars[11] = "x12",
             vars[12] = "x13", vars[13] = "x14", vars[14] = "x15", vars[15] = "x16",
             vars[16] = "x17", vars[17] = "x18", vars[18] = "x19", vars[19] = "x20";

             printf("FAIL\n");
             
             printf("Aliasing test\n");
             printf("ord = "); mpoly_ordering_print(ord);
             printf(", len = %ld, exp_bits = %ld, exp_bound = %lx, "
                                      "coeff_bits = %ld, nvars = %ld\n\n",
                                  len, exp_bits, exp_bound, coeff_bits, nvars);

             fmpz_mpoly_print_pretty(f, vars, ctx); printf("\n\n");
             fmpz_mpoly_print_pretty(g, vars, ctx); printf("\n\n");
          
             flint_printf("c = %wd\n", c);

             flint_abort();
          }
       }

       fmpz_mpoly_clear(f, ctx);  
       fmpz_mpoly_clear(g, ctx);  
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

