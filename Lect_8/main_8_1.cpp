#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include "random.h"

using namespace std;

/************************************************************************************************/
// Funzioni
/************************************************************************************************/
double error(const vector<double> &sum_prog, const vector<double> &sum_prog2, int n) {
    if (n == 0) return 0;
    return sqrt((sum_prog2[n] - pow(sum_prog[n], 2)) / n);
}

double acceptance_ratio(double prob_old, double prob_new) {
    return min(1., prob_new / prob_old);
}

double wf_square_modulus(double x, double sigma, double mu){
    return exp(-pow(x - mu, 2) / pow(sigma, 2)) + exp(-pow(x + mu, 2) / pow(sigma, 2)) + 
           2 * exp(-(pow(x + mu, 2) + pow(x - mu, 2)) / (2 * pow(sigma, 2)));
}

double wf(double x, double sigma, double mu){
    return exp(-pow(x - mu, 2) / (2 * pow(sigma, 2))) + exp(-pow(x + mu, 2) / (2 * pow(sigma, 2)));
}

double laplace(double x, double sigma, double mu){
    double factor = -1. / (sigma * sigma);
    double frac_minus = pow((x + mu) / sigma, 2);
    double frac_plus = pow((x - mu) / sigma, 2);
    double exp_minus = exp(-frac_minus / 2.);
    double exp_plus = exp(-frac_plus / 2.);
    return factor * ((1. - frac_minus) * exp_minus + (1. - frac_plus) * exp_plus);
}

double kinetic_energy(double x, double sigma, double mu){
    return -0.5 * laplace(x, sigma, mu) / wf(x, sigma, mu);
}

double potential_energy(double x){
    return pow(x, 4) - 2.5 * x * x;
}

double ave_hamiltonian(double x, double sigma, double mu){
    return kinetic_energy(x, sigma, mu) + potential_energy(x);
}

/************************************************************************************************/
// Main
/************************************************************************************************/
int main() {
    Random rnd;
    int seed[4], p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) Primes >> p1 >> p2;
    else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();

    ifstream input("seed.in");
    string property;
    if (input.is_open()) {
        while (!input.eof()) {
            input >> property;
            if (property == "RANDOMSEED") {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else {
        cerr << "PROBLEM: Unable to open seed.in" << endl;
        return 1;
    }

    ifstream params("output_SA_accepted.txt");
    if (!params.is_open()) {
        cerr << "PROBLEM: Unable to open output_SA_accepted.txt" << endl;
        return 1;
    }

    ofstream output_accepted("ave_hamiltonian_sa_steps.data");
    if (!output_accepted.is_open()) {
        cerr << "PROBLEM: Unable to open ave_hamiltonian_sa_steps.data" << endl;
        return 1;
    }
    output_accepted << "# mu\t\tsigma\t\tenergy\t\terror\n";

    string line;
    int count_line = 0;
    while (getline(params, line)) {
        istringstream iss(line);
        double dummy, mu, sigma;
        if (!(iss >> dummy >> mu >> sigma)) continue;

        /************************************************************************************************/
        // Simulation Parameters
        /************************************************************************************************/
        int M = 10000, N = 100, L = M / N;
        int count_accept = 0;
        double step = 3.;
        double x_start = 0.;

        vector<double> ave1, ave2;
        vector<double> sum_prog, sum_prog2, error_prog;

        for (int i = 0; i < N; i++) {
            double sum = 0;
            for (int j = 0; j < L; j++) {
                double x = rnd.Rannyu(-step, step);
                double acceptance = acceptance_ratio(wf_square_modulus(x_start, sigma, mu), wf_square_modulus(x_start + x, sigma, mu));

                if (rnd.Rannyu() <= acceptance) {
                    x_start += x;
                    count_accept++;
                }
                sum += ave_hamiltonian(x_start, sigma, mu);
            }
            ave1.push_back(sum / L);
            ave2.push_back(ave1[i] * ave1[i]);
        }

        for (int i = 0; i < N; i++) {
            double s1 = 0, s2 = 0;
            for (int j = 0; j <= i; j++) {
                s1 += ave1[j];
                s2 += ave2[j];
            }
            sum_prog.push_back(s1 / (i + 1));
            sum_prog2.push_back(s2 / (i + 1));
            error_prog.push_back(error(sum_prog, sum_prog2, i));
        }

        count_line++;
        double final_energy = sum_prog[N - 1];
        double final_error = error_prog[N - 1];
        output_accepted << count_line << "\t" << mu << "\t" << sigma << "\t" << final_energy << "\t" << final_error << endl;
    }

    params.close();
    output_accepted.close();
    rnd.SaveSeed();

    cout << "=================================================================================\n";
    cout << "The program has finished successfully.\n";
    cout << "=================================================================================\n";

    return 0;
}
