
#pragma once

#include "util.h"
#include "fraction.h"

#include <cstdint>
#include <string>

namespace vmatrixlib {

template <class float_type>
std::string fpnum_to_string(float_type num, int p); // no generic body!

template <>
std::string fpnum_to_string(long double num, int p)
{
   char buf[256];
   char tbuf[16];

   // If the number is not an integer...
   if (static_cast<long double>(static_cast<std::int64_t>(num)) != num) {

      const long double absnum = fabsl(num);

      // and it is too big, use the standard notation.
      if (absnum <= 1e-6 || absnum >= 1e12) {
         sprintf(buf, "%.3LE", num);
         return buf;
      }
   }

   if (p == 0) {
      sprintf(buf, "%.0Lf", num);
      return buf;
   }

   sprintf(tbuf, "%%.%iLf", p);
   sprintf(buf, tbuf, num);

   char *ptr = buf;
   int count = 0;

   while (*ptr) {
      if (*ptr++ == '.')
         break;
   }

   assert(*ptr != 0);

   //I've found the dot.
   //ptr is in dot+1 position.

   char *ptr2 = buf + strlen(buf) - 1;

   while (ptr2 >= ptr) {

      if (*ptr2 != '0')
         break;

      ptr2--;
      count++;
   }

   //There are 'count' zeroes of p decimal digits
   //So, I have to show p-count decimal digits
   count = p - count;

   sprintf(tbuf, "%%.%iLf", count);
   sprintf(buf, tbuf, num);

   return buf;
}

template <>
std::string fpnum_to_string(double num, int p)
{
   return fpnum_to_string(static_cast<long double>(num), p);
}

template <class T>
inline std::string to_string(const T& t, int p) {
   /*
   * Generic to_string() with precision parameter 'p': fallback to
   * the standard to_string() since we don't know T.
   */
   return std::to_string(t);
}

template <>
inline std::string to_string(const long double& val, int p) {
   return fpnum_to_string(val, p);
}

template <>
inline std::string to_string(const double& val, int p) {
   return fpnum_to_string(val, p);
}


template <class float_type>
inline std::string to_string(const frac<long long, float_type>& f,
                             int precision = 6)
{
   char buf[256];

   if (f.is_using_fp()) {
      return fpnum_to_string(to_float(f), precision);
   }

   const long long num = f.int_numerator();
   const long long den = f.int_denominator();

   if (den == 1) {
      sprintf(buf, "%lld", num);
      return buf;
   }

   sprintf(buf, "%lld/%lld", den > 0 ? num : -num,
                             den > 0 ? den : -den);
   return buf;
}


template <class float_type>
inline std::string to_string(const frac<int32_t, float_type>& f,
                             int precision = 6)
{
   char buf[256];

   if (f.is_using_fp()) {
      return fpnum_to_string(to_float(f), precision);
   }

   const int32_t num = f.int_numerator();
   const int32_t den = f.int_denominator();

   if (den == 1) {
      sprintf(buf, "%d", num);
      return buf;
   }

   sprintf(buf, "%d/%d", den > 0 ? num : -num,
                         den > 0 ? den : -den);
   return buf;
}

} // namespace vmatrixlib
