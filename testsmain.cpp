
#include <cstring>
#include <iostream>
#include <random>

#include "matrix.h"

using namespace std;
using namespace vmatrixlib;

void testing_float_to_frac()
{
   random_device rdev;
   default_random_engine e(rdev());

   lognormal_distribution<> dist(5.0, 3.5);

   cout << "Testing float to frac... ";
   cout.flush();

   char strRepr[256];
   char newStrRepr[256];

   for (int i = 0; i < 100000; i++) {

      long double rval = dist(e);

      sprintf(strRepr, "%.6Lf", rval);
      rval = atof(strRepr); 

      int64_t num, den;
      bool r = float_to_frac(rval, num, den);

      long double fpval =
         static_cast<long double>(num) / static_cast<long double>(den);

      if (!r || fpval != rval) {

         sprintf(newStrRepr, "%.6Lf", fpval);

         if (!strcmp(strRepr, newStrRepr)) continue;

         printf("[FAIL]\n");
         printf("Failed for %.16Lf:\n", rval);
         printf("orig strRepr: %s\n", strRepr);
         printf("new strRepr:  %s\n", newStrRepr);

         cout << "num = " << num << endl;
         cout << "den = " << den << endl;
         return;
      }
   }

   cout << "[PASS]\n";
}

void testing_triang_matrix()
{
   cout << "Making many random matrixes triangular... ";
   cout.flush();

   for (int i = 0; i < 1000; i++) {
      vmatrix A = vmatrix::random(3, 3, -20, 20, 1, 0.5);
      vmatrix ar = A.make_triangular();
   }

   cout << "[PASS]\n";
}

void testing_inv_matrix()
{
   cout << "Inverting matrixes... ";
   cout.flush();

   for (int i = 0; i < 10000; i++) {

      vmatrix A = vmatrix::random(4, 4, -4, 4, 1, 0.35);

      if (A.determinant() == 0)
         continue;

      vmatrix inv = A.compute_inverse();
      vmatrix invinv = inv.compute_inverse();

      if (invinv != A) {

         cout << "[FAIL]\n";
         cout << "A:\n";
         A.pretty_print();

         cout << endl << endl;
         cout << "A^-1:\n";
         inv.pretty_print();

         cout << "(A^-1)^-1:\n";
         invinv.pretty_print();
         return;
      }

   }

   cout << "[PASS]\n";
}

int main(int argc, char ** argv) {

   cout << "sizeof long double: " << sizeof(long double) << endl;
   cout << "eps of long double: "
        << numeric_limits<long double>::epsilon() << endl << endl;

   cout << "One random matrix...\n";

   vmatrix m = vmatrix::random(6, 6, -20, 20, 1, 0.35);

   m.pretty_print();
   //m.print_mathematica_style();
   //m.print_matlab_style();

   testing_float_to_frac();
   testing_triang_matrix();
   testing_inv_matrix();

   //getchar();
   return 0;
}