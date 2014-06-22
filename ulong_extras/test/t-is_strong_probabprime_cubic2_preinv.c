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

int main(int argc, char **argv)
{
   long s1, s2;
   mp_limb_t n, ninv;
   FLINT_TEST_INIT(state);
   
   flint_printf("is_strong_probabprime_cubic2_preinv....");
   fflush(stdout);

   s1 = atol(argv[1]) | 1;
   s2 = atol(argv[2]);
   
   printf("%ld %ld\n", s1, s2);

   for (n = s1; n < s2; n+=2) /* Test that primes pass the test */
   {
      ninv = n_preinvert_limb(n);
      
      if (!n_is_prime(n) && (n % 3) != 0 && (n % 5) != 0 && (n % 7) != 0
         && n_is_strong_probabprime_cubic2_preinv(n, ninv, 2)
         && n_is_strong_probabprime_cubic2_preinv(n, ninv, 3)
         && n_is_strong_probabprime_cubic2_preinv(n, ninv, 5)
         && n_is_strong_probabprime_cubic2_preinv(n, ninv, 7))
      {
         printf("%lu is declared prime\n", n);
      }

      if ((n & 268435455) == 0 || ((n + 1) & 268435455) == 0)
         printf("n = %ld\n", n);
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
