#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

double error(const vector<double> &sum_prog,const vector<double> &sum_prog2, int n){
   double e;
   if (n==0) e=0;
   else e=sqrt((sum_prog2[n]-pow(sum_prog[n],2))/n);  
   
   return e;
}

double chisquare(const vector<int> &counts, double exp_counts){
   double chi=0.;
   double dim_counts = counts.size();
   for (int i=0; i<dim_counts; i++){
      chi += pow((counts[i]- exp_counts),2)/exp_counts;
   }
   return chi;
}

int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   // Constants and needed vectors
   /****************************************************************/
   int M=10000; // # of throws
   int N=100; // # of blocks
   int L=M/N; // # of throws/block

   vector<double> ave1, ave2, var1, var2; 
   vector<double> sum_prog, sum_prog2, error_prog;
   vector<double> sum_var_prog, sum_var_prog2, error_var_prog;
   /****************************************************************/

   // Cycling -- Integral evaluation
   /****************************************************************/
   for(int i=0; i<N; i++){
      double sum=0;
      double sum_var=0.;
      for(int j=0; j<L; j++){
         sum+=rnd.Rannyu();
         sum_var+=pow(rnd.Rannyu()-0.5,2); 
      }
      ave1.push_back(sum/L);
      ave2.push_back(ave1[i]*ave1[i]);          
      var1.push_back(sum_var/L);
      var2.push_back(var1[i]*var1[i]);    
   }

   for(int i=0; i<N; i++){
     sum_prog.push_back(0);
     sum_prog2.push_back(0);
     sum_var_prog.push_back(0);
     sum_var_prog2.push_back(0);

     for(int j=0; j<i+1; j++){
         sum_prog[i] += ave1[j];
         sum_prog2[i] += ave2[j];
         sum_var_prog[i] += var1[j];
         sum_var_prog2[i] += var2[j];
      }
     sum_prog[i] /= (i+1); 
     sum_prog2[i] /= (i+1);
     sum_var_prog[i] /= (i+1); 
     sum_var_prog2[i] /= (i+1);
     error_prog.push_back(error(sum_prog,sum_prog2,i));
     error_var_prog.push_back(error(sum_var_prog,sum_var_prog2,i));
   }
   /****************************************************************/

   // Writing on output file
   /****************************************************************/
   ofstream output_average;
   ofstream output_variance;
   output_average.open("output_average.data");
   output_variance.open("output_variance.data");
   for(int i=0; i<N; i++){
      output_average <<(i+1)*L << "\t" << sum_prog[i] << "\t" << error_prog[i] << endl;
      output_variance <<(i+1)*L << "\t" << sum_var_prog[i] << "\t" << error_var_prog[i] << endl;
   }
   output_average.close();
   output_variance.close();
   /****************************************************************/

   // Terminal
   /****************************************************************/
   cout << "\n\nNumber of blocks = " << N << endl;
   cout << "Number of throws/block = " << L << endl;
   cout << "Number of total throws = " << M << endl;
   cout << "\n\nBest integral value = " << sum_prog[N-1] << " +- " << error_prog[N-1] << "\n\n" << endl;
   /****************************************************************/

   // Chi-square test
   /****************************************************************/
   int n=10000;              // Number of throws
   int m=100;                // Number of intervals in [0,1)
   int exp_counts = n / m;   // Number of expected counts in each interval
   vector<int> test_cases = {100, 1000, 10000};

   // Calculation
   /****************************************************************/
   for (int numtest : test_cases) {
      vector<double> chi_squared_values;
      string output_counts_file = "output_counts_" + to_string(numtest) + ".data";
      string output_chi_file = "output_chi_" + to_string(numtest) + ".data";

      ofstream output_counts(output_counts_file);
        
      for (int i = 0; i < numtest; i++) {
         vector<int> counts(m, 0.);
         for (int j = 0; j < n; j++) {
            int index = static_cast<int>(rnd.Rannyu() * m);   // Index of the interval containing rnd.Rannyu()
            counts[index]++;
         }
         double chi_squared = chisquare(counts, exp_counts);
         chi_squared_values.push_back(chi_squared);
         for (int j = 0; j < m; j++) {
            output_counts << counts[j] << endl;
         }
      }
      output_counts.close();
      ofstream output_chi(output_chi_file);
      for (double chi_value : chi_squared_values) {
          output_chi << chi_value << endl;
      }
      output_chi.close();
   }
rnd.SaveSeed(); 
return 0;
}
