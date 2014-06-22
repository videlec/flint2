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

    Copyright (C) 2014 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int n_is_strong_probabprime_cubic2_preinv(mp_limb_t n, mp_limb_t ninv, mp_limb_t a) 
{
    slong i, j;
    mp_limb_t p0, p1, t;
    mp_limb_t pow, pow2;

    t = n - 1;
    
    i = 0;
    do {
       umul_ppmm(p1, p0, t, UWORD(6148914691236517205));

       if (p0 + p1 + 1 == 0) /* is divisible by 3 */
       {
          i++;
          t = p1 + 1;
       } else 
          break;
    } while (1);

    pow = n_powmod2_preinv(a, t, n, ninv);

    if (pow == 1)
       return 1;

    for (j = 0; j < i; j++)
    {
       pow2 = n_mulmod2_preinv(pow, pow, n, ninv); /* raise to power 3 */
       pow2 = n_mulmod2_preinv(pow2, pow, n, ninv);

       if (pow2 == 1)
          break;
       else
          pow = pow2;
    }

    if (j == i)
       return 0;

    pow2 = n_mulmod2_preinv(pow, pow, n, ninv);
    pow = n_addmod(pow2, pow, n);

    if (pow == n - 1)
       return 1;

    return 0;
}

