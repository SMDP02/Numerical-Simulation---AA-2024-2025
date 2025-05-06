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
   int M=10000;     // # of throws
   int N=100;       // # of blocks
   int L=M/N;       // # of throws/block
   double l = 0.8;  // legth of the needle
   double d = 1.;   // grid spacing
   double alpha, y_center;

   vector<double> ave1, ave2, var1, var2; 
   vector<double> sum_prog, sum_prog2, error_prog;
   vector<double> sum_var_prog, sum_var_prog2, error_var_prog;
   /****************************************************************/
   // Cycling -- Buffon experiment
   /****************************************************************/
   for(int i=0; i<N; i++){
   double N_hit = 0;                                      
      for(int j=0; j<L; j++){
         alpha = rnd.Angle();
         y_center = rnd.Rannyu(0,1);

         // Hitting the grid
         double y_edge_sup = y_center + (l / 2.) * sin(alpha);
         double y_edge_inf = y_center - (l / 2.) * sin(alpha);

         if (y_edge_sup > d || y_edge_inf < 0) N_hit++;
      }
   ave1.push_back((2 * l * L )/ (d * N_hit));
   ave2.push_back(ave1[i]*ave1[i]); 
   }

   for(int i=0; i<N; i++){
    sum_prog.push_back(0);
    sum_prog2.push_back(0);

    for(int j=0; j<i+1; j++){
        sum_prog[i] += ave1[j];
        sum_prog2[i] += ave2[j];
    }
    sum_prog[i] /= (i+1); 
    sum_prog2[i] /= (i+1);
    error_prog.push_back(error(sum_prog,sum_prog2,i)); 
   }
   /****************************************************************/
   // Writing on output file
   /****************************************************************/
   ofstream output_average;
   output_average.open("output_average.data");
   for(int i=0; i<N; i++){
      output_average <<(i+1)*L << "\t" << sum_prog[i] << "\t" << error_prog[i] << endl;
   }
   output_average.close();
   /****************************************************************/

   // Terminal
   /****************************************************************/
   cout << "\n\nNumber of blocks = " << N << endl;
   cout << "Number of throws/block = " << L << endl;
   cout << "Number of total throws = " << M << endl;
   cout << "\n\nBest integral value = " << sum_prog[N-1] << " +- " << error_prog[N-1] << "\n\n" << endl;
   /****************************************************************/
rnd.SaveSeed(); 
return 0;
}
