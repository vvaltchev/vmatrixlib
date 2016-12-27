
#pragma once

#include <string>
#include <limits>

#include <cmath>
#include <cstdint>
#include <stdexcept>

#include "util.h"

namespace vmatrixlib {

template <class integerT, class floatT>
class frac {

public:

   typedef integerT integer_type;
   typedef floatT float_type;

protected:

   // members
   integer_type num = 0;
   integer_type den = 1;
   float_type fn = 0.0;


   // static functions

protected: // methods

   frac neg() const {

      if (!is_using_fp()) {
         return frac(-num, den);
      }

      frac res;
      res.fn = -fn;
      return res;
   }

   float_type fpval() const {
      return fn != 0.0 ? fn : (float_type)num / (float_type)den;
   }

public:

   frac() = default;

   frac(integer_type n, integer_type d)
      : num(n), den(d), fn(0.0)
   {
      if (d == 0) {
         throw std::domain_error("Division by zero!");
      }
   }

   frac(float_type number, int precision = 6) : frac()
   {
      if (float_to_frac(number, num, den, precision)) {
         semplify();
      } else {
         fn = number;
      }
   }

   static frac make_dec_frac(const float_type& n) {
      frac res;
      res.fn = n;
      return res;
   }

   frac& operator=(const frac& rhs) = default;

   frac to_frac_in_decimal_form() const {

      frac res;
      res.fn = fpval();
      return res;
   }

   bool is_using_fp() const {
      return fn != 0.0;
   }

   frac operator+(const frac& f2) const;
   frac operator-(const frac& f2) const { return operator+(f2.neg()); }
   frac operator*(const frac& f2) const;

   frac operator/(const frac& f2) const {

      frac f22 = !f2.is_using_fp()
                    ? frac(f2.den, f2.num)
                    : make_dec_frac(1.0 / f2.fn);

      return operator*(f22);
   }

   frac operator-() const { return neg(); }

   bool operator==(const frac& f2) const {
      float_type eps = 10 * std::numeric_limits<float_type>::epsilon();
      return std::abs(fpval() - f2.fpval()) <= eps;
   }

   bool operator!=(const frac& f2) const { return !operator==(f2); }

   frac& operator+=(const frac& f2) { return *this = operator+(f2); }
   frac& operator-=(const frac& f2) { return *this = operator-(f2); }
   frac& operator*=(const frac& f2) { return *this = operator*(f2); }
   frac& operator/=(const frac& f2) { return *this = operator/(f2); }

   float_type numerator() const { return fn != 0.0 ? fn : num; }
   float_type denominator() const { return fn != 0.0 ? 1.0 : den; }

   integer_type int_numerator() const {
      assert(!is_using_fp());
      return num;
   }

   integer_type int_denominator() const {
      assert(!is_using_fp());
      return den;
   }

   operator std::string() const;
   operator float_type() const { return fpval(); }
   operator integer_type() const {
      return !is_using_fp()
         ? static_cast<integer_type>(num / den)
         : static_cast<integer_type>(fn);
   }


   void semplify();
};

template <class integer_type, class float_type>
void frac<integer_type, float_type>::semplify()
{
   if (is_using_fp()) {
      return;
   }

   if (den < 0) {
      num = -num;
      den = -den;
   }

   integer_type d = gcd(num, den);
   num /= d;
   den /= d;
}


template <class integer_type, class float_type>
struct fp_type_of<frac<integer_type, float_type>> {
   typedef float_type type;
};


template <class integer_type, class float_type>
inline frac<integer_type, float_type>
to_frac_in_decimal_form(const frac<integer_type, float_type>& f) {
   return f.to_frac_in_decimal_form();
}

template <class T, class U>
auto numerator(const frac<T, U>& val) {
   return val.numerator();
}


template <class T, class U>
auto denominator(const frac<T, U>& val) {
   return val.denominator();
}


template <class integer_type, class float_type>
inline frac<integer_type, float_type>
operator+(float_type num, const frac<integer_type, float_type>& f) {
   return frac<integer_type, float_type>(num)+f;
}

template <class integer_type, class float_type>
inline frac<integer_type, float_type>
operator-(float_type num, const frac<integer_type, float_type>& f) {
   return frac<integer_type, float_type>(num)-f;
}

template <class integer_type, class float_type>
inline frac<integer_type, float_type>
operator*(float_type num, const frac<integer_type, float_type>& f) {
   return frac<integer_type, float_type>(num)*f;
}

template <class integer_type, class float_type>
inline frac<integer_type, float_type>
operator/(float_type num, const frac<integer_type, float_type&> f) {
   return frac<integer_type, float_type>(num)/f;
}


template <class integer_type, class float_type>
frac<integer_type, float_type>
frac_semplify(const frac<integer_type, float_type>& f)
{
   frac<integer_type, float_type> f2 = f;
   f2.semplify();
   return f2;
}

template <class integer_type, class float_type>
frac<integer_type, float_type>
frac<integer_type, float_type>::operator+(const frac& f2) const
{
   frac res;
   frac f11, f22;
   float_type sum;

   f11 = frac_semplify(*this);
   f22 = frac_semplify(f2);

   sum = to_float(f11) + to_float(f22);

   if (f11.is_using_fp() || f22.is_using_fp() ||
       num_out_of_range<integer_type>(sum) ||
       num_out_of_range<integer_type>(f11.numerator() * f22.denominator()) ||
       num_out_of_range<integer_type>(f11.numerator() * f22.denominator() +
       f22.numerator() * f11.denominator())) {
      return frac::make_dec_frac(sum);
   }

   const integer_type num =
      f11.int_numerator() * f22.int_denominator() +
      f22.int_numerator() * f11.int_denominator();

   const integer_type den = f11.int_denominator() * f22.int_denominator();

   return frac_semplify(frac(num, den));

}

template <class integer_type, class float_type>
frac<integer_type, float_type>
frac<integer_type, float_type>::operator*(const frac& f2) const
{
   frac res;
   frac f11, f22;
   float_type prod;

   f11 = frac_semplify(*this);
   f22 = frac_semplify(f2);

   prod = to_float(f11) * to_float(f22);

   if (f11.is_using_fp() || f22.is_using_fp() ||
       num_out_of_range<integer_type>(prod) ||
       num_out_of_range<integer_type>(f11.numerator() * f22.numerator()) ||
       num_out_of_range<integer_type>(f11.denominator() * f22.denominator())) {
      return frac::make_dec_frac(prod);
   }

   const integer_type num = f11.int_numerator() * f22.int_numerator();
   const integer_type den = f11.int_denominator() * f22.int_denominator();

   return frac_semplify(frac(num, den));
}


template <class integer_type, class float_type>
inline std::string to_string(const frac<integer_type, float_type>& f,
                             int precision = 6); // no generic body!



template <class integer_type, class float_type>
frac<integer_type, float_type>::operator std::string() const
{
   return to_string(*this);
}

} // namespace vmatrixlib
