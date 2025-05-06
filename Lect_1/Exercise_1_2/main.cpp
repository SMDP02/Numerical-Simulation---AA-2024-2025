#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

int main(int argc, char *argv[]) {
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
            rnd.SetRandom(seed,p1,p2);   // Fixing random seed for reproducibility
            }
        } 
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    /****************************************************************/
    // Constants and required vectors
    /****************************************************************/
    const int M = 10000;
    const double lambda = 1;
    const double gamma = 1;
    const double mu = 0; 
    
    vector<int> N = {1, 2, 10, 100};
    
    vector<ofstream> outputFiles;
    vector<string> filenames = {"output_N1.txt", "output_N2.txt", "output_N10.txt", "output_N100.txt"};
    
    for (size_t j = 0; j < N.size(); j++) {
        outputFiles.emplace_back(filenames[j]);
        if (!outputFiles[j].is_open()) {
            cerr << "PROBLEM: Unable to open " << filenames[j] << endl;
            return 1;
        }
    }
    
    /****************************************************************/
    // Calculation
    /****************************************************************/
    for (size_t j = 0; j < N.size(); j++) {
        for (int i = 0; i < M; i++) {
            double sum_unif = 0.0, sum_exp = 0.0, sum_cauchylorentz = 0.0;
            
            // Generate values and compute averages
            for (int k = 0; k < N[j]; k++) {
                sum_unif += rnd.Rannyu();
                sum_exp += rnd.Exponential(lambda);
                sum_cauchylorentz += rnd.CauchyLorentz(gamma, mu);
            }
            sum_unif /= N[j];
            sum_exp /= N[j];
            sum_cauchylorentz /= N[j];
            
            // Write values to respective file
            outputFiles[j] << sum_unif << " " << sum_exp << " " << sum_cauchylorentz << "\n";
        }
        outputFiles[j].close();
    }
    
    rnd.SaveSeed();
    return 0;
}
