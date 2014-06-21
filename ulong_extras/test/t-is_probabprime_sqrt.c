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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   long i, result;
   mp_limb_t d;
   mpz_t d_m;
   slong pow;
   ulong bits;
   FLINT_TEST_INIT(state);
   

   flint_printf("is_probabprime_sqrt....");
   fflush(stdout);
   
   for (i = 1UL; i < 1000000000000L; i++) /* Test that primes pass the test */
   {
      if (!n_is_prime(i) && n_is_probabprime_sqrt(i))
      {
         printf("%lu is declared prime\n", i);
      }

      if ((i & 1073741823) == 0)
         printf("i = %ld\n", i);
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
