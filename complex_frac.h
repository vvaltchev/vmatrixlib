
#pragma once

#include "fraction.h"
#include "to_string.h"

namespace vmatrixlib {

template <class frac_type>
class complex_frac {

public:

   // static methods

   static complex_frac zero() {
      return complex_frac();
   }

   static complex_frac one() {
      return complex_frac(1, 0);
   }

protected:

   frac_type re;
   frac_type im;

public:

   complex_frac() : re(0), im(0) { }

   complex_frac(const frac_type& r) : re(r), im(0.0) { }
   complex_frac(const frac_type& r, const frac_type& i) : re(r), im(i) { }

   complex_frac neg() { return complex_frac(-re, -im); }
   complex_frac conj() { return complex_frac(re, -im); }

   complex_frac operator+(const complex_frac& c2) const {
      return complex_frac(re + c2.re, im + c2.im);
   }

   complex_frac operator-(const complex_frac& c2) const {
      return complex_frac(re - c2.re, im - c2.im);
   }

   complex_frac operator*(const complex_frac& c2) const;
   complex_frac operator/(const complex_frac& c2) const;
   complex_frac operator-() { return neg(); }

   bool operator==(const complex_frac& c2) const {
      return re == c2.re && im == c2.im;
   }

   bool operator!=(const complex_frac& c2) const {
      return !operator==(c2);
   }

   bool operator==(const frac_type& c2) const {
      return operator==(complex_frac(c2));
   }

   bool operator!=(const frac_type& c2) const {
      return operator!=(complex_frac(c2));
   }

   complex_frac& operator=(const frac_type& n) {
      return *this = complex_frac(n);
   }

   complex_frac& operator+=(const complex_frac& c2) {
      return *this = operator+(c2);
   }

   complex_frac& operator-=(const complex_frac& c2) {
      return *this = operator+(c2);
   }

   complex_frac& operator*=(const complex_frac& c2) {
      return *this = operator*(c2);
   }

   complex_frac& operator/=(const complex_frac& c2) {
      return *this = operator/(c2);
   }


   frac_type real_part() const { return re; }
   frac_type imag_part() const { return im; }

   operator std::string() const;

   bool operator<(const complex_frac& c2) const {
      return re.fpval() < c2.re.fpval();
   }
   bool operator<=(const complex_frac& c2) const {
      return re.fpval() <= c2.re.fpval();
   }
   bool operator>(const complex_frac& c2) const {
      return re.fpval() > c2.re.fpval();
   }
   bool operator>=(const complex_frac& c2) const {
      return re.fpval() >= c2.re.fpval();
   }
};


template <class T>
inline complex_frac<T> to_frac_in_decimal_form(const complex_frac<T>& c) {
   return complex_frac<T>(to_frac_in_decimal_form(c.real_part()),
                          to_frac_in_decimal_form(c.imag_part()));
}

template <class frac_type>
inline std::ostream& operator<<(std::ostream& s,
                                const complex_frac<frac_type>& c)
{
   return s << to_string(c);
}

template <class frac_type>
complex_frac<frac_type>
complex_frac<frac_type>::operator*(const complex_frac& c2) const {

   frac_type _re, _im;

   if (im == frac_type(0) && c2.im == frac_type(0)) {
      return complex_frac(re*c2.re, 0);
   }

   _re = re*c2.re - im*c2.im;
   _im = im*c2.re + re*c2.im;

   return complex_frac(frac_semplify(_re), frac_semplify(_im));
}

template <class frac_type>
complex_frac<frac_type>
complex_frac<frac_type>::operator/(const complex_frac& c2) const {

   frac_type _re, _im, div;

   if (im == frac_type(0) && c2.im == frac_type(0)) {
      return complex_frac(re / c2.re, 0);
   }

   div = (c2.re) * (c2.re) + (c2.im) * (c2.im);

   _re = (re*c2.re + im*c2.im) / div;
   _im = (im*c2.re - re*c2.im) / div;

   return complex_frac(frac_semplify(_re), frac_semplify(_im));
}

/*
 * Converts a complex_frac to a (pretty) string.
 *
 * NOTE: it is assumed that:
 * - the float_type used by frac_type is convertible to long double
 * - OR fabsl(x) has an overload for float_type. 
 */

template <class frac_type>
std::string to_string(const complex_frac<frac_type>& val,
                      int precision = 6)
{
   typedef decltype(to_float(frac_type(0))) float_type;

   std::string realStr;
   std::string imgStr;
   std::string res;

   const frac_type re = val.real_part();
   frac_type im = val.imag_part();

   const float_type eps = powl(10, -precision);

   bool imNegative = false;

   if (to_float(re) != 0.0) {
      realStr = to_string(re, precision);
   }

   if (to_float(im) != 0.0) {

      if (fabsl(to_float(im) - 1.0) < eps) {

         imgStr = "i";

      } else if (fabsl(to_float(im) + 1.0) < eps) {


         if (realStr.size()) {

            imgStr = "i";
            imNegative = true;

         } else {

            imgStr = "-i";
         }

      } else {

         if (to_float(im) < 0.0 && realStr.size() != 0) {

            imNegative = true;
            im = -im;
         }

         imgStr = to_string(im, precision);

         if (denominator(im) == 1.0) {

            imgStr = imgStr + "i";

         } else {

            imgStr = "(" + imgStr + ")i";
         }

      }

   }

   if (!realStr.size() && !imgStr.size())
      return "0";

   if (realStr.size() && !imgStr.size())
      return realStr;

   if (!realStr.size() && imgStr.size())
      return imgStr;

   if (!imNegative)
      res = realStr + "+" + imgStr;
   else
      res = realStr + "-" + imgStr;

   return res;
}

template <class frac_type>
complex_frac<frac_type>::operator std::string() const {
   return to_string(*this);
}

} // namespace vmatrixlib
