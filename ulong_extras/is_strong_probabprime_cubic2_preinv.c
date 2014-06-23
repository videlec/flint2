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
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2014 Dana Jacobsen
    Copyright (C) 2014 William Hart

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int n_is_strong_probabprime_cubic2_preinv(mp_limb_t n, mp_limb_t ninv, mp_limb_t a) 
{
    slong i2, i3, j2, j3;
    mp_limb_t p0, p1, t2, t3;
    mp_limb_t pow, pow2, pow3;
    mp_limb_t p2, p3;

    if (n <= 1)
       return 0;
    
    if (a >= n)
       a = n_mod2_preinv(a, n, ninv);

    t2 = n - 1;

    i2 = i3 = 0;
    
    while ((t2 & 1) == 0)
    {
       i2++;
       t2 >>= 1;
    }
    
    t3 = t2;
    
    p3 = 1;
    
    do {
       umul_ppmm(p1, p0, t3, UWORD(6148914691236517205));

       if (p0 + p1 + 1 == 0) /* is divisible by 3 */
       {
          i3++;
          t3 = p1 + 1;
          p3 *= 3;
       } else 
          break;
    } while (1);

    pow3 = n_powmod2_preinv(a, t3, n, ninv);

    pow2 = n_powmod2_preinv(pow3, p3, n, ninv);
 
    p2 = 1;
    if (pow2 != 1)
    {
       for (j2 = 0; j2 < i2 && pow2 != n - 1; j2++)
       {
          p2 *= 2;
          pow2 = n_mulmod2_preinv(pow2, pow2, n, ninv);
       }

       if (j2 == i2)
       {
          return 0; /* not a strong probable prime to base a */
       }
    }

    /* pow2 is now either +/- 1 mod n */

    pow = n_powmod2_preinv(pow3, p2, n, ninv);

    if (pow == pow2)
       return 1;

    for (j3 = 0; j3 < i3; j3++)
    {
       pow3 = n_mulmod2_preinv(pow, pow, n, ninv); /* raise to power 3 */
       pow3 = n_mulmod2_preinv(pow3, pow, n, ninv);

       if (pow3 == pow2)
          break;
       else
          pow = pow3;
    }

    pow3 = n_mulmod2_preinv(pow, pow, n, ninv);
    if (pow2 == 1)
       pow = n_addmod(pow3, pow, n);
    else
       pow = n_submod(pow3, pow, n);

    if (pow == n - 1)
       return 1;

    return 0;
}

