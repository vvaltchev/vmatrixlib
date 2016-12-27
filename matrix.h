
#pragma once

#include <cassert>
#include <vector>
#include <random>
#include "complex_frac.h"

namespace vmatrixlib {

template <class T>
class matrix {

public:
   typedef T number_type;

protected:

   int _rows;
   int _cols;
   int _rowSwapsCount;
   std::vector<T> _data;

public:

   static matrix random(int rows, int cols, int min,
                        int max, int decimals, double zero_prob);

   matrix();
   matrix(int rows, int cols);
   matrix(int rows, int cols, T *arr);

   void load_data(T *arr);

   int rows() const { return _rows; }
   int cols() const { return _cols; }
   int size() const { return _rows*_cols; }
   bool is_square() const { return _rows == _cols; }


   void clear();
   void make_identity();

   const T& operator()(int r, int c) const;
   const T& operator()(int index) const;
   const T& get(int r, int c) const;
   const T& get(int n) const;

   T& operator()(int r, int c);
   T& operator()(int index);
   T& get(int r, int c);
   T& get(int n);

   void in_place_mul_by_constant(const T& n);
   void in_place_div_by_constant(const T& n);
   void in_place_mul(const matrix &m);

   void in_place_transpose();
   void in_place_sum(const matrix &m);
   void in_place_mul_row(int row, const T& k);
   void in_place_div_row(int row, const T& k);

   matrix transpose() const;
   matrix operator+(const matrix& m) const;
   matrix operator-(const matrix& m) const;
   matrix operator*(const T& n) const;
   matrix operator*(const matrix& m) const;
   bool operator==(const matrix& m) const;
   bool operator!=(const matrix& m) const {
      return !operator==(m);
   }

   matrix& operator*=(const T& n) {
      in_place_mul_by_constant(n);
      return *this;
   }

   matrix& operator/=(const T& n) {
      in_place_div_by_constant(n);
      return *this;
   }

   void swap(int i, int j, int x, int y) {
      const T tmp = get(i, j);
      get(i, j) = get(x, y);
      get(x, y) = tmp;
   }

   void swap_rows(int i, int j);

   bool is_lower_triangular() const;
   bool is_upper_triangular() const;

   bool is_triangular() const {
      return is_lower_triangular() || is_upper_triangular();
   }

   void add_row_mult_by_const_to_row(int srcRow, int destRow, T k);

   bool has_row_echelon_form() const;
   matrix make_triangular() const;
   T diagonal_product() const;

   int rank() const;
   T determinant() const;
   matrix compute_inverse() const;

   matrix sub_matrix_erasing_row_col(int r, int c) const;
   matrix row_reduce() const;
   matrix null_space() const;
   matrix col_space(matrix *cols = nullptr) const;

   matrix add_row(const matrix& row);
   matrix add_col(const matrix& col);

   void attach_sub_matrix(const matrix& m, int row, int col);
   void attach_col(const matrix& srcMatrix, int srcCol, int destCol);
   void attach_row(const matrix& srcMatrix, int srcRow, int destRow);

   matrix approx_matrix() const;

   bool is_row_null(int row) const;
   bool is_col_null(int col) const;

   int find_elem_in_col(int col, const T& elem) const;
   int find_elem_in_row(int row, const T& elem) const;
   int find_elem(const T& elem) const;

   void pretty_print(int precision = 6) const;
   void print_mathematica_style() const;
   void print_matlab_style() const;
};

// Fast matrix instantiation using double.
typedef matrix<double> fast_vmatrix;

// Slower, but more precise matrix instantiation.
typedef matrix<complex_frac<frac<long long, long double>>> vmatrix;


template <class T>
inline T& matrix<T>::get(int r, int c) {

   assert(r >= 0 && r < _rows && c >= 0 && c < _cols);
   return _data[r*_cols + c];
}

template <class T>
inline const T& matrix<T>::get(int r, int c) const {

   assert(r >= 0 && r < _rows && c >= 0 && c < _cols);
   return _data[r*_cols + c];
}

template <class T>
inline T& matrix<T>::get(int n) {

   assert(n >=0 && n <= size());
   return _data[n];
}

template <class T>
inline const T& matrix<T>::get(int n) const {

   assert(n >= 0 && n <= size());
   return _data[n];
}

template <class T>
inline const T& matrix<T>::operator()(int r, int c) const {

   assert(r >= 0 && r < _rows && c >= 0 && c < _cols);
   return _data[r*_cols + c];
}

template <class T>
inline T& matrix<T>::operator()(int r, int c) {

   assert(r >= 0 && r < _rows && c >= 0 && c < _cols);
   return _data[r*_cols + c];
}

template <class T>
inline T& matrix<T>::operator()(int index) {

   assert(index >= 0 && index < _rows*_cols);
   return _data[index];
}

template <class T>
inline const T& matrix<T>::operator()(int index) const {

   assert(index >= 0 && index < _rows*_cols);
   return _data[index];
}

template <class T>
inline matrix<T> operator*(T n, matrix<T> m) { return m*n; }

template <class T>
matrix<T>::matrix() : _rows(0), _cols(0), _rowSwapsCount(0) { }

template <class T>
matrix<T>::matrix(int r, int c)
   : _rows(r), _cols(c), _rowSwapsCount(0), _data(_rows * _cols)
{
   clear();
}

template <class T>
matrix<T>::matrix(int r, int c, T *arr)
   : _rows(r), _cols(c), _rowSwapsCount(0), _data(_rows, _cols)
{
   load_data(arr);
}

template <class T>
void matrix<T>::clear() {

   _rowSwapsCount=0;

   for (int i=0; i < _rows*_cols; i++)
      _data[i]=0;
}

template <class T>
void matrix<T>::make_identity() {

   if (!is_square())
      throw std::domain_error("The matrix isn't a square matrix");

   clear();

   for (int i=0; i < _rows; i++)
      get(i,i)=1;
}

template <class T>
void matrix<T>::in_place_mul_by_constant(const T& n) {

   for (int i=0; i < _rows*_cols; i++)
      _data[i] *= n;
}

template <class T>
void matrix<T>::in_place_div_by_constant(const T& n) {

   for (int i = 0; i < _rows*_cols; i++)
      _data[i] /= n;
}


template <class T>
void matrix<T>::in_place_mul(const matrix &m) {

   if (_rows != m._rows || _cols != m._cols)
      throw std::domain_error("Argument matrix MUST have the same size as object matrix");

   for (int i=0; i < _rows; i++)
      for (int j=0; j < _cols; j++)
         for (int k=0; k < _cols; k++)
            get(i,j)+=get(i,k)*m.get(k,j);

}

template <class T>
void matrix<T>::in_place_transpose() {

   if (!is_square())
      throw std::domain_error("in_place_transpose() can be used ONLY for square matrices");

   for (int i=0; i < _rows; i++)
      for (int j=0; j < _cols; j++)
         if (j < i)
            swap(i,j,j,i);
}

template <class T>
void matrix<T>::in_place_sum(const matrix &m) {

   if (_rows != m._rows || _cols != m._cols)
      throw std::domain_error("Argument matrix and object matrix MUST have the same size");

   for (int i=0; i < _rows; i++)
      for (int j=0; j < _cols; j++)
         get(i,j)+=m(i,j);
}

template <class T>
matrix<T> matrix<T>::operator+(const matrix<T>& m) const {

   matrix res(*this);
   res.in_place_sum(m);

   return res;
}

template <class T>
matrix<T> matrix<T>::operator-(const matrix<T>& m) const {

   matrix res = *this;
   matrix t = m * T(-1);
   res.in_place_sum(t);

   return res;
}

template <class T>
matrix<T> matrix<T>::operator*(const T& n) const {

   matrix res = *this;
   res.in_place_mul_by_constant(n);

   return res;
}

template <class T>
matrix<T> matrix<T>::operator*(const matrix<T>& m) const {

   if (m._rows != _cols)
      throw std::domain_error("Right matrix must have rows count equals to first matrix's columns count");

   int resR = _rows;
   int resC = m._cols;

   matrix res(resR,resC);

   for (int i=0; i < resR; i++)
      for (int j=0; j < resC; j++)
         for (int k=0; k < _cols; k++)
            res(i,j) += get(i,k)*m(k,j);

   return res;
}

template <class T>
bool matrix<T>::operator==(const matrix& m) const {

   if (_rows != m.rows() || _cols != m.cols())
      return false;

   for (int i=0; i < size(); i++)
      if (_data[i] != m._data[i])
         return false;

   return true;
}

template <class T>
matrix<T> matrix<T>::transpose() const {

   matrix res(_cols,_rows);

   for (int i=0; i < _rows; i++)
      for (int j=0; j < _cols; j++)
         res(j,i)=get(i,j);

   return res;
}

template <class T>
bool matrix<T>::is_lower_triangular() const {

   if (!is_square())
      return false;

   int i,j;

   for (i=0; i < _rows; i++)
      for (j=0; j < _cols; j++)
         if (i < j && get(i,j) != 0)
            return false;

   return true;
}

template <class T>
bool matrix<T>::is_upper_triangular() const {

   if (!is_square())
      return false;

   int i,j;

   for (i=0; i < _rows; i++)
      for (j=0; j < _cols; j++)
         if (i > j && get(i,j) != 0)
            return false;

   return true;
}

template <class T>
void matrix<T>::swap_rows(int i, int j) {

   for (int k=0; k < _cols; k++)
      swap(i,k,j,k);

   _rowSwapsCount++;
}

template <class T>
void matrix<T>::add_row_mult_by_const_to_row(int srcRow, int destRow, T k) {

   for (int i=0; i < _cols; i++)
      get(destRow,i)+=get(srcRow,i)*k;
}

template <class T>
matrix<T> matrix<T>::make_triangular() const {

   if (rows() == 1 || cols() == 1 || has_row_echelon_form())
      return *this;

   matrix res = *this;

   int i,j,k,u;

   i=0; j=0;

   while (i < res._rows && j < res._cols) {

      for (k=i; k < res._rows; k++)
         if (res(k,j) != 0)
            break;



      if (k == res._rows) {
         //we did not find a row 'k' with elem k,j != 0
         j++;
         continue;
      }


      if (k != i) {
         res.swap_rows(k, i);
      }

      T val = res(i,j);

      for (u=i+1; u < res._rows; u++) {

         if (res(u,j) == 0)
            continue;

         res.add_row_mult_by_const_to_row(i, u, -res(u,j)/val);

         //forced zero
         res(u,j)=0;
      }


      i++;
      j++;
   }


   if (!res.has_row_echelon_form()) {
      res.pretty_print();
      throw std::runtime_error("Calc error in make_triangular()");
   }

   return res;

}

template <class T>
bool matrix<T>::has_row_echelon_form() const {

   int i=0,j=0;
   bool canIncSteps=true;


   while (i < _rows) {

      if (j == _cols) {

         j=0;
         i++;
         continue;
      }

      if (get(i,j) == 0) {

         int k;
         for (k=i; k < _rows; k++)
            if (get(k,j) != 0)
               return false;

         j++;
         canIncSteps=true;

      } else {

         if (!canIncSteps)
            return false;

         canIncSteps=false;
         i++;
      }

   }

   return true;
}

template <class T>
T matrix<T>::diagonal_product() const {

   if (!is_square())
      throw std::domain_error("Diagonal product can be done only for square matrices");

   T res = T(1);

   for (int i=0; i < _rows; i++)
      res*=get(i,i);

   return res;
}

template <class T>
T matrix<T>::determinant() const {

   if (!is_square())
      throw std::domain_error("Determinant can be computed only for square matrices");

   if (_rows == 1)
      return get(0,0);

   if (_rows == 2)
      return get(0,0)*get(1,1) - get(0,1)*get(1,0);

   if (is_triangular()) {

      T det = diagonal_product();

      if (_rowSwapsCount == 0 || (_rowSwapsCount%2) == 0)
         return det;

      return -det;
   }

   matrix m = make_triangular();

   T det = m.diagonal_product();

   int rSwaps = m._rowSwapsCount;

   if (rSwaps == 0 || (rSwaps%2) == 0)
      return det;

   return -det;

}

template <class T>
int matrix<T>::rank() const {

   matrix m = make_triangular();

   int i=0,j=0;
   int steps=0;
   bool canIncSteps=true;

   while (i < m._rows) {

      if (j == m._cols) {

         j=0;
         i++;
         continue;
      }

      if (m.get(i,j) == 0) {

         j++;
         canIncSteps=true;

      } else {

         assert(canIncSteps);
         canIncSteps=false;
         steps++;
         i++;
      }
   }

   return steps;
}

template <class T>
matrix<T> matrix<T>::compute_inverse() const {

   T det = determinant();

   if (det == 0)
      throw std::runtime_error("Can't invert a singular matrix");

   matrix res(rows(),cols());

   int i,j;

   for (i=0; i < rows(); i++) {

      for (j=0; j < cols(); j++) {

         int sgn,sgnI=1,sgnJ=1;

         if (((i+1)%2) == 0) sgnI=-1;
         if (((j+1)%2) == 0) sgnJ=-1;

         //sgn = (-1)^(i+j)
         sgn=sgnI*sgnJ;

         matrix subM = sub_matrix_erasing_row_col(i, j);

         res(i,j) = T(sgn) * subM.determinant();
      }
   }

   res.in_place_transpose();
   res.in_place_div_by_constant(det);
   return res;
}

template <class T>
matrix<T> matrix<T>::sub_matrix_erasing_row_col(int row, int col) const {

   matrix res;

   if (row != -1 && col != -1)
      res = matrix(_rows-1,_cols-1);
   else if (row == -1 && col != -1)
      res = matrix(_rows,_cols-1);
   else if (row != -1 && col == -1)
      res = matrix(_rows-1,_cols);
   else
      return *this;

   int currR=0;
   int currC=0;
   int i,j;

   for (i=0; i < _rows; i++) {

      if (i == row)
         continue;

      currC=0;

      for (j=0; j < _cols; j++) {

         if (j == col)
            continue;

         res(currR,currC) = get(i,j);

         currC++;
      }

      currR++;
   }

   return res;

}

template <class T>
void matrix<T>::in_place_mul_row(int row, const T& k) {

   for (int i=0; i < _cols; i++)
      get(row,i) *= k;
}

template <class T>
void matrix<T>::in_place_div_row(int row, const T& k) {

   for (int i = 0; i < _cols; i++)
      get(row, i) /= k;
}


template <class T>
matrix<T> matrix<T>::row_reduce() const {


   if (_rows == 1 || _cols == 1) {
      return *this;
   }

   matrix t = make_triangular();

   int j=0;

   for (int i=_rows-1; i >= 0; i--) {

      j=i;

      if (j >= _cols)
         j=_cols-1;

      T f = t(i,j);


      if (f == 0) {

         for (j=i+1; j < _cols; j++) {

            if (t(i,j) != 0) {
               f=t(i,j);
               break;
            }
         }


         if (f == 0) {
            continue;
         }
      }

      t.in_place_div_row(i, f);

      //forced one (to overcome to approx errors)
      t(i,j)=1;

      if (i == 0)
         break;


      for (int k=i-1; k >= 0; k--) {

         if (t(k,j) == 0)
            continue;

         t.add_row_mult_by_const_to_row(i, k, -t(k,j));

         //forced zero (to overcome to approx errors)
         t(k,j)=0;

      }
   }

   return t;
}

template <class T>
void matrix<T>::attach_sub_matrix(const matrix<T>& m, int row, int col) {

   for (int i=0; i+row < _rows && i < m._rows; i++)
      for (int j=0; j+col < _cols && j < m._cols; j++)
         get(i+row,j+col) = m(i,j);

}

template <class T>
void matrix<T>::attach_col(const matrix<T>& srcMatrix, int srcCol, int destCol) {

   if (srcMatrix.rows() != rows())
      throw std::domain_error("srcMatrix must have the same number of rows as destination matrix");

   for (int i=0; i < rows(); i++)
      get(i,destCol) = srcMatrix(i,srcCol);
}

template <class T>
void matrix<T>::attach_row(const matrix<T>& srcMatrix, int srcRow, int destRow) {

   if (srcMatrix.cols() != cols())
      throw std::domain_error("srcMatrix must have the same number of cols as dest matrix");

   for (int i=0; i < cols(); i++)
      get(destRow,i) = srcMatrix(srcRow,i);
}

template <class T>
matrix<T> matrix<T>::add_row(const matrix<T>& row) {

   matrix res(_rows+1, _cols);

   if (row._rows != 1)
      throw std::domain_error("Row matrix MUST have only ONE row");

   if (row._cols != _cols)
      throw std::domain_error("The number of cols must be the same");

   res.attach_sub_matrix(*this, 0, 0);
   res.attach_sub_matrix(row, _rows, 0);

   return res;
}

template <class T>
matrix<T> matrix<T>::add_col(const matrix<T>& col) {

   matrix res(_rows, _cols+1);

   if (col._cols != 1)
      throw std::domain_error("Col matrix MUST have only ONE column");

   if (col._rows != _rows)
      throw std::domain_error("The number of rows must be the same");

   res.attach_sub_matrix(*this, 0, 0);
   res.attach_sub_matrix(col, 0, _cols);

   return res;
}

template <class T>
matrix<T> matrix<T>::null_space() const {

   matrix r = row_reduce();

   int kerDim=0;
   int depVarsCount=0;

   std::vector<int> depVarsRows(_rows);
   std::vector<int> indepVars(_cols);

   for (int i=0; i < _rows; i++) {

      int j;

      depVarsRows[i]=-1;

      for (j=0; j < _cols; j++) {

         if (r(i,j) == 1)
            break;
      }

      if (j != _cols) {
         depVarsRows[i]=j;
         depVarsCount++;
      }
   }


   for (int i=0; i < _cols; i++) {

      // check if the var 'i' is among depVarRows

      int j;
      for (j=0; j < _rows; j++)
         if (depVarsRows[j] == i)
            break;

      if (j == _rows) {
         // we found nothing: adding it to indepVar
         indepVars[kerDim++]=i;
      }
   }

   matrix res(_cols, kerDim);

   for (int k=0; k < kerDim; k++) {

      // building the coloumn 'indepVars[k]' for the independent variable 'k'.

      for (int i=0; i < _cols; i++) {

         if (i == indepVars[k]) {
            res(i,k)=1;
            continue;
         }

         int j;
         for (j=0; j < _rows; j++)
            if (depVarsRows[j] == i)
               break;

         res(i, k) = (j != _rows ? -r(j, indepVars[k]) : T(0));
      }

   }

   return res;
}

template <class T>
matrix<T> matrix<T>::col_space(matrix *cols) const {

   matrix r = row_reduce();

   int depVarsCount=0;
   std::vector<int> depVars(_cols);

   for (int i=0; i < _rows; i++) {

      int j;

      for (j=0; j < _cols; j++)
         if (r(i,j) == 1)
            break;


      if (j != _cols)
         depVars[depVarsCount++]=j;
   }

   matrix res(_rows, depVarsCount);

   for (int i=0; i < _rows; i++)
      for (int j=0; j < depVarsCount; j++)
         res(i,j) = get(i,depVars[j]);

   if (cols) {

      *cols = matrix(depVarsCount,1);

      for (int i=0; i < depVarsCount; i++)
         (*cols)(i,0)=depVars[i];
   }

   return res;
}

template <class T>
matrix<T> matrix<T>::random(int rows, int cols, int min,
                            int max, int decimals, double zero_prob)
{
   using namespace std;

   assert(zero_prob >= 0.0f && zero_prob <= 1.0f);
   assert(decimals <= 6);

   zero_prob *= 1000;

   long double den = 1;

   for (int k = 0; k < decimals; k++)
      den *= 10.0;

   random_device rdev;
   default_random_engine e(rdev());

   uniform_int_distribution<> dist(min * (int)den, max * (int)den);
   uniform_int_distribution<> zerodist(0, 999);

   matrix res(rows,cols);

   for (int i=0; i < rows; i++) {
      for (int j=0; j < cols; j++) {

         if (zerodist(e) < zero_prob) {
            res(i,j) = 0;
            continue;
         }

         res(i, j) = T( static_cast<long double>(dist(e)) / den );
      }
   }

   return res;
}

template <class T>
int matrix<T>::find_elem_in_col(int col, const T& elem) const {

   int i;

   for (i=0; i < rows(); i++)
      if (get(i,col) == elem)
         return i;

   return -1;
}

template <class T>
int matrix<T>::find_elem_in_row(int row, const T& elem) const {

   int i;

   for (i=0; i < cols(); i++)
      if (get(row,i) == elem)
         return i;

   return -1;
}

template <class T>
int matrix<T>::find_elem(const T& elem) const {

   int i;

   for (i=0; i < size(); i++)
      if (_data[i] == elem)
         return i;

   return -1;
}

template <class T>
void matrix<T>::pretty_print(int precision) const {


   int i,j,k,maxlen=0;
   printf("matrix (%i x %i)\n", _rows, _cols);

   if (!_rows || !_cols)
      return;

   std::vector<std::string> strmatrix(_rows * _cols);

#define getelem(r,c) (strmatrix[((r) * _cols + (c))])

   for (i=0; i < _rows; i++) {
      for (j=0; j < _cols; j++) {

         T num = get(i,j);

         const std::string& s = to_string(num, precision);
         getelem(i,j) = s + " | ";

         if ((int)s.size() > maxlen)
            maxlen = s.size();
      }
   }

   printf("\n");

   for (i=0; i < _rows; i++) {

      for (j=0; j < _cols; j++) {

         if (j==0)
            printf("| ");

         int diff = maxlen - (getelem(i,j).size() - 3);

         for (k=0; k < diff; k++)
            printf(" ");

         printf("%s", getelem(i,j).c_str());
      }

      printf("\n");
   }

   printf("\n");

#undef getelem
}

template <class T>
matrix<T> matrix<T>::approx_matrix() const {

   matrix r(_rows,_cols);

   for (int i=0; i < _rows; i++)
      for (int j=0; j < _cols; j++)
         r(i,j) = to_frac_in_decimal_form(get(i,j));

   return r;
}

template <class T>
void matrix<T>::load_data(T *arr) {

   T *ptr = &_data[0];

   for (int i=0; i < _rows*_cols; i++)
      *ptr++ = *arr++;
}

template <class T>
bool matrix<T>::is_row_null(int row) const {

   for (int i=0; i < _cols; i++)
      if (get(row,i) != 0)
         return false;

   return true;
}

template <class T>
bool matrix<T>::is_col_null(int col) const {

   for (int i=0; i < _rows; i++)
      if (get(i,col) != 0)
         return false;

   return true;
}

template <class T>
void matrix<T>::print_mathematica_style() const {

   printf("{");

   for (int i=0; i < _rows; i++) {

      printf("{");

      for (int j=0; j < _cols; j++) {

         printf("%s", to_string(get(i,j)).c_str());

         if (j < _cols-1)
            printf(",");
      }

      printf("}");

      if (i < _rows-1)
         printf(",\n");
   }

   printf("}\n");
}

template <class T>
void matrix<T>::print_matlab_style() const {


   printf("[");

   for (int i=0; i < _rows; i++) {


      for (int j=0; j < _cols; j++) {

         printf("%s", to_string(to_frac_in_decimal_form(get(i,j))).c_str());

         if (j < _cols-1)
            printf(" ");
      }


      if (i < _rows-1)
         printf(";\n");
   }

   printf("]\n");
}

} // namespace vmatrixlib
