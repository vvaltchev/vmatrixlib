
#pragma once

#include <cstdint>
#include <limits>
#include <algorithm>

#include <stdio.h>

namespace vmatrixlib {

template <class integer_type, class float_type>
inline bool num_out_of_range(float_type val)
{
   return val <= std::numeric_limits<integer_type>::min() ||
          val >= std::numeric_limits<integer_type>::max();
}


template <class integer_type>
integer_type gcd(integer_type a, integer_type b)
{
   int shift;
   integer_type u, v;
   integer_type diff;

   u = a >= 0 ? a : -a;
   v = b >= 0 ? b : -b;

   /* GCD(0,x) := x */
   if (u == 0 || v == 0)
      return u | v;

   /* Let shift := lg K, where K is the greatest power of 2
   dividing both u and v. */
   for (shift = 0; ((u | v) & 1) == 0; ++shift) {
      u >>= 1;
      v >>= 1;
   }

   while ((u & 1) == 0) {
      u >>= 1;
   }

   /* From here on, u is always odd. */
   do {
      while ((v & 1) == 0)  /* Loop X */
         v >>= 1;

      /* Now u and v are both odd, so diff(u, v) is even.
      Let u = min(u, v), v = diff(u, v)/2. */

      if (u < v) {

         v -= u;

      }
      else {

         diff = u - v;
         u = v;
         v = diff;
      }

      v >>= 1;

   } while (v != 0);


   return u << shift;
}


template <class integer_type, class float_type>
bool float_to_frac(float_type number,
                   integer_type& num, integer_type& den,
                   int precision = 6); // no generic body!


template <class integer_type>
bool float_to_frac(long double number,
                   integer_type& num,
                   integer_type& den,
                   int precision = 6)
{
   typedef integer_type inttype;
   constexpr const inttype lim = std::numeric_limits<inttype>::max();
   constexpr const int log10hi = std::numeric_limits<inttype>::digits10;
   static const long double log10floathi =
      log10l(std::numeric_limits<long double>::max());

   const int sign = number >= 0.0 ? 1 : -1;
   number = fabsl(number);

   long double int_part;
   const long double frac_part = modfl(number, &int_part);

   if (frac_part == 0.0) {

      if (int_part >= lim)
         return false;

      num = sign * static_cast<inttype>(int_part);
      den = 1;
      return true;
   }

   if (number == 0.0 || number == 1.0) {
      num = static_cast<inttype>(number);
      den = 1;
      return true;
   }

   const long double log10num = log10l(number);
   const int log10num_int =
      static_cast<int>(log10num >= 0.0 ? ceill(log10num) : floorl(log10num));

   long double lscale = log10hi - std::max(0, log10num_int);

   if (lscale < 0.0) {

      // The number is too big to fit in a 'integer_type', even without
      // considering its fractional part.

      return false;
   }

   // In the very *unlikely* case when we're using big integers which can
   // represent much bigger numbers when a long double, use the long double's
   // limit.

   if (lscale >= log10floathi) {
      lscale = log10floathi - 1;
   }

   // Artificially limit the scale to 10^6, in order to significantly reduce
   // the amount of floating point artifacts for 'reasonable numbers'.
   lscale = std::min(static_cast<long double>(precision), lscale);

   long double scale = powl(10.0, lscale);
   assert(number * scale < lim);

   num = sign * static_cast<inttype>(roundl(number * scale));
   den = static_cast<inttype>(roundl(scale));
   return true;
}


template <class integer_type>
bool float_to_frac(double number,
                   integer_type& num,
                   integer_type& den,
                   int precision = 6)
{
   return float_to_frac(static_cast<long double>(number), num, den, precision);
}


template <class T>
inline T frac_semplify(const T& t) {
   return t;
}

template <class T>
struct fp_type_of {
   typedef T type;
};

template <class T>
auto to_float(const T& t) {
   return static_cast<typename fp_type_of<T>::type>(t);
}

template <class T>
inline T to_frac_in_decimal_form(const T& t) {
   return t;
}

template <class T>
T numerator(const T& val) {
   return val;
}

template <class T>
T denominator(const T& val) {
   return 1;
}

} // namespace vmatrixlib
