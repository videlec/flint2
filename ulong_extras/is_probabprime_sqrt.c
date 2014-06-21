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

******************************************************************************/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_sqrtmod_correct(mp_limb_t a, mp_limb_t p) 
{
    slong i, r, m;
    mp_limb_t p1, k, b, g, bpow, gpow, res, t;
    mp_limb_t pinv;

    if (a <= 1)
    {
        return 1; /* 0^2 = 0 and 1^2 = 1 mod p */
    }

    pinv = n_preinvert_limb(p);

    t = p - 1;
    while ((t & 1) == 0)
       t >>= 1;
    
    if (!n_is_strong_probabprime2_preinv(p, pinv, a, t))
       return 0;

    if ((p & UWORD(3)) == 3)
    {
        res = n_powmod2_ui_preinv(a, (p + 1)/4, p, pinv); /* p == 2^B - 1 isn't prime */

        return a == n_mulmod2_preinv(res, res, p, pinv);
    }

    if ((p & UWORD(7)) == 5)
    {
       b = n_powmod2_ui_preinv(a, (p + 3)/8, p, pinv); /* p == 2^B - 3 isn't prime */
       g = n_mulmod2_preinv(b, b, p, pinv);

       if (g == a)
       {
          return 1;
       }

       g = n_powmod2_ui_preinv(2, (p - 1)/4, p, pinv);
       res = n_mulmod2_preinv(g, b, p, pinv);
        
       return a == n_mulmod2_preinv(res, res, p, pinv);
    }

    r = 0;
    p1 = p - 1;

    do {
        p1 >>= UWORD(1); 
        r++;
    } while ((p1 & UWORD(1)) == 0);

    b = n_powmod2_ui_preinv(a, p1, p, pinv);

    for (k = 3; ; k+=2) /* 2 is a quadratic residue mod p = 8k + 1 */
    {
        if (n_jacobi_unsigned(k, p) == -1) break;
    }
    
    g = n_powmod2_ui_preinv(k, p1, p, pinv);
    res = n_powmod2_ui_preinv(a, (p1 + 1) / 2, p, pinv);

    while (b != 1 && r > 1)
    {
        bpow = b;
        m = 0;
        do
        {
            bpow = n_mulmod2_preinv(bpow, bpow, p, pinv);
            m++;
        } while (m < r && bpow != 1);

        if (m == r)
        {
           return 0;
        }

        gpow = g;
        for (i = 1; i < r - m; i++)
        {
            gpow = n_mulmod2_preinv(gpow, gpow, p, pinv);
        }
        res = n_mulmod2_preinv(res, gpow, p, pinv);
        g = n_mulmod2_preinv(gpow, gpow, p, pinv);
        b = n_mulmod2_preinv(b, g, p, pinv);
        r = m;
    }

    if (r <= 1 && b != 1)
    {
       return 0;
    }

    return a == n_mulmod2_preinv(res, res, p, pinv);
}

int n_is_probabprime_sqrt(ulong n)
{
   ulong a, iters, maxiters;
   
   if (n < 11) {
      if (n == 2 || n == 3 || n == 5 || n == 7)   return 1;
      else                                        return 0;
   }
   
   if (!(n%2) || !(n%3) || !(n%5) || !(n%7))       return 0;
    
   if (n <  121) /* 11*11 */                       return 1;
    
   if (!(n%11) || !(n%13) || !(n%17) || !(n%19) ||
       !(n%23) || !(n%29) || !(n%31) || !(n%37) ||
       !(n%41) || !(n%43) || !(n%47) || !(n%53))   return 0;
   
   if (n < 3481) /* 59*59 */                       return 1;

   if (n_is_square(n))
      return 0;

   a = 2;
   
   if (n < 14981)
      maxiters = 1;
   else if (n < 486737)
      maxiters = 2;
   else if (n < 159874021)
      maxiters = 3;
   else if (n < 331658081)
      maxiters = 4;
   else if (n < 3056100623)
      maxiters = 5;
   else
      maxiters = 6;
    
   for (iters = 0; iters < maxiters; iters++)
   {
      for (; a < n; a++)
      {
         if (n_jacobi_unsigned(a, n) == 1 && !n_is_square(a))
            break;
      }

      if (n_sqrtmod_correct(a, n) == 0)
         return 0;

      a++;
   }

   return 1;
}

