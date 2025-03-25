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

double S(double S0, double r, double sigma, double Z, double t){
   return S0 * exp((r - 1/2. * sigma * sigma)*t + sigma * Z *sqrt(t));
}

double call(double r, double t, double S, double K){
   return exp(-r*t) * max(0., S-K);
}

double put(double r, double t, double S, double K){
   return exp(-r*t) * max(0., K-S);
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
   // Data blocking
   int M=100000; // # of throws
   int N=100; // # of blocks
   int L=M/N; // # of throws/block
   // Option pricing
   double S0 = 100.;
   double K = 100.;
   double T = 1.;
   double r = 0.1;
   double sigma = 0.25;

   vector<double> ave1_call, ave1_put, ave2_call, ave2_put; 
   vector<double> sum_prog_call, sum_prog2_call, error_prog_call;
   vector<double> sum_prog_put, sum_prog2_put, error_prog_put;
   /****************************************************************/
   // 
   /****************************************************************/
   
   for(int i=0; i<N; i++){
      double st=0;
      double c=0;
      double p=0;
      
      for(int j=0; j<L; j++){
         double Z = rnd.Gauss(0,1);
         st = S(S0, r, sigma, Z, T);
         c+=call(r,T,st,K);
         p+=put(r,T,st,K);
      }
   
      ave1_call.push_back(c/L);
      ave2_call.push_back(ave1_call[i]*ave1_call[i]);    
      ave1_put.push_back(p/L);
      ave2_put.push_back(ave1_put[i]*ave1_put[i]);
   }
   for(int i=0; i<N; i++){
   // Not calculating at the 1st step
     sum_prog_call.push_back(0);
     sum_prog2_call.push_back(0);
     sum_prog_put.push_back(0);
     sum_prog2_put.push_back(0);

     for(int j=0; j<i+1; j++){
         sum_prog_call[i] += ave1_call[j];
         sum_prog2_call[i] += ave2_call[j];
         sum_prog_put[i] += ave1_put[j];
         sum_prog2_put[i] += ave2_put[j];
      }
     sum_prog_call[i] /= (i+1); 
     sum_prog2_call[i] /= (i+1);
     sum_prog_put[i] /= (i+1); 
     sum_prog2_put[i] /= (i+1);
     error_prog_call.push_back(error(sum_prog_call,sum_prog2_call,i));
     error_prog_put.push_back(error(sum_prog_put,sum_prog2_put,i));
   }
   /****************************************************************/
   // Writing on output file
   /****************************************************************/

   ofstream output_call;
   ofstream output_put;
   output_put.open("output_put.data");
   output_call.open("output_call.data");
   for(int i=0; i<N; i++){
      output_call <<(i+1)*L << "\t" << sum_prog_call[i] << "\t" << error_prog_call[i] << endl;
      output_put <<(i+1)*L << "\t" << sum_prog_put[i] << "\t" << error_prog_put[i] << endl;
   
   }
   output_call.close();
   output_put.close();

   vector<double> ave1_call_discretized, ave1_put_discretized, ave2_call_discretized, ave2_put_discretized; 
   vector<double> sum_prog_call_discretized, sum_prog2_call_discretized, error_prog_call_discretized;
   vector<double> sum_prog_put_discretized, sum_prog2_put_discretized, error_prog_put_discretized;
   /****************************************************************/
   // 
   /****************************************************************/
   
   for(int i=0; i<N; i++){
      double c_discretized=0;
      double p_discretized=0;
      
      for(int j=0; j<L; j++){
         double st_discretized=S0;
         for(int k=0; k<100; k++){
            double Z_discretized = rnd.Gauss(0,1);
            st_discretized = S(st_discretized, r, sigma, Z_discretized, T/100.);
         }
         c_discretized += call(r,T,st_discretized,K);
         p_discretized += put(r,T,st_discretized,K);
      }
   
      ave1_call_discretized.push_back(c_discretized/L);
      ave2_call_discretized.push_back(ave1_call_discretized[i]*ave1_call_discretized[i]);    
      ave1_put_discretized.push_back(p_discretized/L);
      ave2_put_discretized.push_back(ave1_put_discretized[i]*ave1_put_discretized[i]);
   }
   for(int i=0; i<N; i++){
   // Not calculating at the 1st step
     sum_prog_call_discretized.push_back(0);
     sum_prog2_call_discretized.push_back(0);
     sum_prog_put_discretized.push_back(0);
     sum_prog2_put_discretized.push_back(0);

     for(int j=0; j<i+1; j++){
         sum_prog_call_discretized[i] += ave1_call_discretized[j];
         sum_prog2_call_discretized[i] += ave2_call_discretized[j];
         sum_prog_put_discretized[i] += ave1_put_discretized[j];
         sum_prog2_put_discretized[i] += ave2_put_discretized[j];
      }
     sum_prog_call_discretized[i] /= (i+1); 
     sum_prog2_call_discretized[i] /= (i+1);
     sum_prog_put_discretized[i] /= (i+1); 
     sum_prog2_put_discretized[i] /= (i+1);
     error_prog_call_discretized.push_back(error(sum_prog_call_discretized,sum_prog2_call_discretized,i));
     error_prog_put_discretized.push_back(error(sum_prog_put_discretized,sum_prog2_put_discretized,i));
   }
   /****************************************************************/
   // Writing on output file
   /****************************************************************/

   ofstream output_call_discretized;
   ofstream output_put_discretized;
   output_put_discretized.open("output_put_discretized.data");
   output_call_discretized.open("output_call_discretized.data");
   for(int i=0; i<N; i++){
      output_call_discretized <<(i+1)*L << "\t" << sum_prog_call_discretized[i] << "\t" << error_prog_call_discretized[i] << endl;
      output_put_discretized <<(i+1)*L << "\t" << sum_prog_put_discretized[i] << "\t" << error_prog_put_discretized[i] << endl;
   }
   output_call_discretized.close();
   output_put_discretized.close();

rnd.SaveSeed(); 
return 0;
}
