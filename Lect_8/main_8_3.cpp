#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include "random.h"

using namespace std;
/************************************************************************************************/
// Functions
/************************************************************************************************/

double error(const vector<double> &sum_prog, const vector<double> &sum_prog2, int n) {
    if (n == 0) return 0;
    return sqrt((sum_prog2[n] - pow(sum_prog[n], 2)) / n);
}

double acceptance_ratio(double prob_old, double prob_new) {
   return min(1., prob_new / prob_old);
}

double wf_square_modulus(double x, double sigma, double mu){
    return exp(-pow(x - mu, 2) / (pow(sigma, 2))) + exp(-pow(x + mu, 2) / (pow(sigma, 2))) + 2 * exp( - (pow(x + mu, 2) + pow(x - mu, 2))/ (2 * pow(sigma, 2)));
}

double wf(double x, double sigma, double mu){
    return exp(-pow(x - mu, 2) / (2*pow(sigma, 2))) + exp(-pow(x + mu, 2) / (2*pow(sigma, 2)));
}

double laplace(double x, double sigma, double mu){
    double factor = -1./(sigma*sigma);
    double frac_minus = pow((x + mu)/sigma,2);
    double frac_plus = pow((x - mu)/sigma,2);
    double exp_minus= exp(-1. * frac_minus/2.);
    double exp_plus = exp(-1. * frac_plus/2.);
    return factor * ((1. - frac_minus) * exp_minus + (1. - frac_plus) * exp_plus);
}

double kinetic_energy(double x, double sigma, double mu){
    return - 1./2. * laplace(x, sigma, mu)/wf(x, sigma, mu); // Natural units
}

double potential_energy(double x){
    return pow(x, 4) - 2.5 * x * x;
}

double ave_hamiltonian(double x, double sigma, double mu){
    return kinetic_energy(x, sigma, mu) + potential_energy(x);
}


double Metropolis_aveEnergy(double mu, double sigma, Random &rnd) {
   int M=100000;
   double ave_energy=0.;
   double x=0.;
   double passo = 2.;

   for (int i=0; i<M; i++) {

      double x1=rnd.Rannyu(-passo, passo);
      double alpha = min(1., wf_square_modulus(x+x1, sigma, mu)/wf_square_modulus(x, sigma, mu));
      
      if (rnd.Rannyu()<=alpha) x=x+x1;
      ave_energy += kinetic_energy(x, sigma, mu)+potential_energy(x);
   }     
   return ave_energy/M;
}

double Temperature (int step) {
   return 3./(1.+4.*step);
}

double Acceptance_Ratio_SA (double energy_old, double energy_new, double T) {
   return min(1., exp(-1./T*(energy_new-energy_old)));
}

void DataBlocking_wf_square_modulus(int M, int N, double sigma, double mu, Random &rnd){
    int count_accept = 0;
    double x_start = 0.;
    int L = M/N;
    double step = 3.;
    vector<double> ave1, ave2;
    vector<double> sum_prog, sum_prog2, error_prog;

    ofstream sampling_data("sampling.data");
    ofstream ave_hamiltonian_data("ave_hamiltonian.data");
    for (int i = 0; i < N; i++) {
        double sum = 0;
        for (int j = 0; j < L; j++) {
            double x = rnd.Rannyu(-step, step);
            double acceptance = acceptance_ratio(wf_square_modulus(x_start, sigma, mu), wf_square_modulus(x_start + x, sigma, mu));
        
            if (rnd.Rannyu() <= acceptance) {
               x_start = x_start + x;
               count_accept++;
            }
            sampling_data << x_start << endl;

            sum += ave_hamiltonian(x_start, sigma, mu);
        }
        ave1.push_back(sum / L);
        ave2.push_back(ave1[i] * ave1[i]);
    }

    // Data blocking
    for (int i = 0; i < N; i++) {
        sum_prog.push_back(0);
        sum_prog2.push_back(0);
        for (int j = 0; j < i + 1; j++) {
            sum_prog[i] += ave1[j];
            sum_prog2[i] += ave2[j];
        }
        sum_prog[i] /= (i + 1);
        sum_prog2[i] /= (i + 1);
        error_prog.push_back(error(sum_prog, sum_prog2, i));
        ave_hamiltonian_data << (i + 1) * L << "\t" << sum_prog[i] << "\t" << error_prog[i] << "\t" << endl;
    }
}
// **********************************************************************************************/
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
   
   DataBlocking_wf_square_modulus(1000000, 100, 0.618567, 0.803221, rnd);
   rnd.SaveSeed();
   return 0;
}

