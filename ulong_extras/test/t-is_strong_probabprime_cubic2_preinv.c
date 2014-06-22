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

    Copyright (C) 2009, 2014 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   long result;
   mp_limb_t n, ninv, t;
   FLINT_TEST_INIT(state);
   
   flint_printf("is_strong_probabprime_cubic2_preinv....");
   fflush(stdout);
   
   for (n = 3UL; n < 1000000000000L; n+=2) /* Test that primes pass the test */
   {
      ninv = n_preinvert_limb(n);
      
      t = n - 1;
      while ((t & 1) == 0)
         t >>= 1;

      if (!n_is_prime(n) && (n % 3) == 1 && (n % 5) != 0 && (n % 7) != 0
         && n_is_strong_probabprime2_preinv(n, ninv, 2, t)
         && n_is_strong_probabprime2_preinv(n, ninv, 3, t)
         && n_is_strong_probabprime2_preinv(n, ninv, 5, t)
         && n_is_strong_probabprime_cubic2_preinv(n, ninv, 7))
         /* && n_is_strong_probabprime_cubic2_preinv(n, ninv, 7) */
         /* && n_is_strong_probabprime2_preinv(n, ninv, 2, t) */
      {
         printf("%lu is declared prime\n", n);
      }

      if ((n & 16777215) == 0 || ((n + 1) & 16777215) == 0)
         printf("n = %ld\n", n);
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
