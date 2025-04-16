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

double function_cosine(double x){
   return (M_PI/2.) * cos((M_PI * x)/2.);
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
   /************************************************************************************************/
   int M=100000;   // # of throws
   int N=100;      // # of blocks
   int L=M/N;      // # of throws/block

   vector<double> ave1, ave2, var1, var2, ave1_impsampling, ave2_impsampling; 
   vector<double> sum_prog, sum_prog2, error_prog;
   vector<double> sum_var_prog, sum_var_prog2, error_var_prog;
   vector<double> sum_prog_impsampling, sum_prog2_impsampling, error_prog_impsampling;
   /************************************************************************************************/
   // Cycling -- Integral evaluation
   /************************************************************************************************/
   for(int i=0; i<N; i++){
      double sum=0;
      double sum_var=0.;
      double sum_impsampling = 0.;
      for(int j=0; j<L; j++){
         double x_uniform = rnd.Rannyu();
         sum += function_cosine(x_uniform);
         sum_var += pow(x_uniform - 1, 2);
         
         double x_imp = rnd.Retta();
         sum_impsampling += function_cosine(x_imp) / (-2.0 * x_imp + 2.0);       
   }
      ave1.push_back(sum/L);
      ave2.push_back(ave1[i]*ave1[i]);          
      var1.push_back(sum_var/L);
      var2.push_back(var1[i]*var1[i]);
      ave1_impsampling.push_back(sum_impsampling/L);
      ave2_impsampling.push_back(ave1_impsampling[i]*ave1_impsampling[i]);
   }

   for(int i=0; i<N; i++){
   // Not calculating at the 1st step
     sum_prog.push_back(0);
     sum_prog2.push_back(0);
     
     sum_var_prog.push_back(0);
     sum_var_prog2.push_back(0);
     sum_prog_impsampling.push_back(0);
     sum_prog2_impsampling.push_back(0);

     for(int j=0; j<i+1; j++){
         sum_prog[i] += ave1[j];
         sum_prog2[i] += ave2[j];
         sum_var_prog[i] += var1[j];
         sum_var_prog2[i] += var2[j];
         sum_prog_impsampling[i] += ave1_impsampling[j];
         sum_prog2_impsampling[i] += ave2_impsampling[j];
      }
     sum_prog[i] /= (i+1); 
     sum_prog2[i] /= (i+1);
     sum_var_prog[i] /= (i+1); 
     sum_var_prog2[i] /= (i+1);
     sum_prog_impsampling[i] /= (i+1);
     sum_prog2_impsampling[i] /= (i+1);

     error_prog.push_back(error(sum_prog,sum_prog2,i));
     error_var_prog.push_back(error(sum_var_prog,sum_var_prog2,i));
     error_prog_impsampling.push_back(error(sum_prog_impsampling,sum_prog2_impsampling,i));
   }
   /************************************************************************************************/
   // Writing on output file
   /************************************************************************************************/
   ofstream output_average;
   ofstream output_variance;
   ofstream output_impsampling;
   output_average.open("output_average.data");
   output_variance.open("output_variance.data");
   output_impsampling.open("output_impsampling.data");
   for(int i=0; i<N; i++){
      output_average <<(i+1)*L << "\t" << sum_prog[i] << "\t" << error_prog[i] << endl;
      output_variance <<(i+1)*L << "\t" << sum_var_prog[i] << "\t" << error_var_prog[i] << endl;
      output_impsampling << (i+1)*L << "\t" << sum_prog_impsampling[i] << "\t" << error_prog_impsampling[i] << endl;
   }
   output_average.close();
   output_variance.close();
   output_impsampling.close();
   /************************************************************************************************/
   // Terminal
   /************************************************************************************************/
   cout << "\n\nNumber of blocks = " << N << endl;
   cout << "Number of throws/block = " << L << endl;
   cout << "Number of total throws = " << M << endl;
   cout << "\n\nBest integral value= " << sum_prog[N-1] << " +- " << error_prog[N-1] << "\n\n" << endl;
   cout << "\n\nBest integral value with importance sampling = " << sum_prog_impsampling[N-1] << " +- " << error_prog_impsampling[N-1] << "\n\n" << endl;

   /************************************************************************************************/
rnd.SaveSeed(); 
return 0;
}
