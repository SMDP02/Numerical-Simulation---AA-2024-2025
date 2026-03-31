#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

/************************************************************************************************/
// Funzioni
/************************************************************************************************/
double error(const vector<double> &sum_prog, const vector<double> &sum_prog2, int n) {
    if (n == 0) return 0;
    return sqrt((sum_prog2[n] - pow(sum_prog[n], 2)) / n);
}

double wf_1s(double x, double y, double z) {
   return 1 / sqrt(M_PI) * exp(-sqrt(x * x + y * y + z * z));
}

double wf_2p(double x, double y, double z) {
   return 1 / 8. * sqrt(2 / M_PI) * exp(-sqrt(x * x + y * y + z * z) / 2.) * z;
}

double acceptance_ratio(double prob_old, double prob_new) {
   return min(1., prob_new / prob_old);
}

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
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

   /************************************************************************************************/
   // Parameters, Vectors, Counters
   /************************************************************************************************/
    int M = 1000000, N = 100, L = M / N;
    // step[0] = half the step for 1S with uniform distribution
    // step[1] = standard deviation for 1S with Gaussian distribution
    // step[2] = half the step for 2P with uniform distribution
    // step[3] = standard deviation for 2P with Gaussian distribution
    vector<double> step = {1.22, 0.7, 2.8, 1.6};
    
    double x_1s, y_1s, z_1s, x_1s_far, y_1s_far, z_1s_far,
           x_1s_start, y_1s_start, z_1s_start,
           x_1s_far_start, y_1s_far_start, z_1s_far_start;
    double x_2p, y_2p, z_2p, x_2p_far, y_2p_far, z_2p_far,
           x_2p_start, y_2p_start, z_2p_start,
           x_2p_far_start, y_2p_far_start, z_2p_far_start;
    int count_accept_1s = 0, count_accept_2p = 0, count_accept_far_1s = 0, count_accept_far_2p = 0;

    // Scelta della distribuzione
    int distribution_choice;
    cout << "Choose a distribution for the random walk:" << endl;
    cout << "1. Uniform" << endl;
    cout << "2. Gaussian" << endl;
    cout << "Your choice (1 or 2): ";
    cin >> distribution_choice;

    cout << "Insert the initial coordinates for the 1s orbital (near (0,0,0)), format: x y z → ";
    cin >> x_1s >> y_1s >> z_1s;

    cout << "Insert the initial coordinates for the 1s orbital (near from (0,0,0)), format: x y z → ";
    cin >> x_1s_far >> y_1s_far >> z_1s_far;

    cout << "Insert the initial coordinates for the 2p orbital (near (0,0,1)), format: x y z → ";
    cin >> x_2p >> y_2p >> z_2p;

    cout << "Insert the initial coordinates for the 2p orbital (far from (0,0,1)), format: x y z → ";
    cin >> x_2p_far >> y_2p_far >> z_2p_far;
    cout << "\n\n=================================================================================\n";
    cout << "\nCALCULATION ..." << endl;

    // Copy the initial positions for output
    x_1s_start = x_1s; y_1s_start = y_1s; z_1s_start = z_1s;
    x_1s_far_start = x_1s_far; y_1s_far_start = y_1s_far; z_1s_far_start = z_1s_far;
    x_2p_start = x_2p; y_2p_start = y_2p; z_2p_start = z_2p;
    x_2p_far_start = x_2p_far; y_2p_far_start = y_2p_far; z_2p_far_start = z_2p_far;    

    // Output files
    string suffix = (distribution_choice == 1) ? "_uniform" : "_gauss";

    ofstream pos_output_1s("positions_1s" + suffix + ".data");
    ofstream pos_output_2p("positions_2p" + suffix + ".data");
    ofstream pos_output_1s_far("positions_1s_far" + suffix + ".data");
    ofstream pos_output_2p_far("positions_2p_far" + suffix + ".data");
    ofstream avg_output_1s("output_average_1s" + suffix + ".data");
    ofstream avg_output_2p("output_average_2p" + suffix + ".data");
    ofstream avg_output_1s_far("output_average_1s_far" + suffix + ".data");
    ofstream avg_output_2p_far("output_average_2p_far" + suffix + ".data");
    

    vector<double> ave1_1s, ave2_1s, sum_prog_1s, sum_prog2_1s, error_prog_1s;
    vector<double> ave1_1s_far, ave2_1s_far, sum_prog_1s_far, sum_prog2_1s_far, error_prog_1s_far;
    vector<double> ave1_2p, ave2_2p, sum_prog_2p, sum_prog2_2p, error_prog_2p;
    vector<double> ave1_2p_far, ave2_2p_far, sum_prog_2p_far, sum_prog2_2p_far, error_prog_2p_far;

   /************************************************************************************************/
   // 1S Orbital
   /************************************************************************************************/
    for (int i = 0; i < N; i++) {
        double sum_1s = 0;
        double sum_1s_far = 0;
        
        for (int j = 0; j < L; j++) {

            double x_new, y_new, z_new, x_new_far, y_new_far, z_new_far;
            if (distribution_choice == 1) {
                x_new = x_1s + rnd.Rannyu(-step[0], step[0]);
                y_new = y_1s + rnd.Rannyu(-step[0], step[0]);
                z_new = z_1s + rnd.Rannyu(-step[0], step[0]);
            
                x_new_far = x_1s_far + rnd.Rannyu(-step[0], step[0]);
                y_new_far = y_1s_far + rnd.Rannyu(-step[0], step[0]);
                z_new_far = z_1s_far + rnd.Rannyu(-step[0], step[0]);
            }
            else if (distribution_choice == 2) {
                x_new = x_1s + rnd.Gauss(0, step[1]);
                y_new = y_1s + rnd.Gauss(0, step[1]);
                z_new = z_1s + rnd.Gauss(0, step[1]);
            
                x_new_far = x_1s_far + rnd.Gauss(0, step[1]);
                y_new_far = y_1s_far + rnd.Gauss(0, step[1]);
                z_new_far = z_1s_far + rnd.Gauss(0, step[1]);
            }
            
            double acceptance = acceptance_ratio(pow(wf_1s(x_1s, y_1s, z_1s), 2), pow(wf_1s(x_new, y_new, z_new), 2));
            double acceptance_far = acceptance_ratio(pow(wf_1s(x_1s_far, y_1s_far, z_1s_far), 2), pow(wf_1s(x_new_far, y_new_far, z_new_far), 2));

            if (rnd.Rannyu() <= acceptance) {
               x_1s = x_new; y_1s = y_new; z_1s = z_new;
               count_accept_1s++;
            }
            if (rnd.Rannyu() <= acceptance_far) {
                x_1s_far = x_new_far; y_1s_far = y_new_far; z_1s_far = z_new_far;
                count_accept_far_1s++;
            }
            pos_output_1s << x_1s << "\t" << y_1s << "\t" << z_1s << endl;
            pos_output_1s_far << x_1s_far << "\t" << y_1s_far << "\t" << z_1s_far << endl;

            sum_1s += sqrt(x_1s * x_1s + y_1s * y_1s + z_1s * z_1s);
            sum_1s_far += sqrt(x_1s_far * x_1s_far + y_1s_far * y_1s_far + z_1s_far * z_1s_far);
        }
        ave1_1s.push_back(sum_1s / L);
        ave2_1s.push_back(ave1_1s[i] * ave1_1s[i]);
        ave1_1s_far.push_back(sum_1s_far / L);
        ave2_1s_far.push_back(ave1_1s_far[i] * ave1_1s_far[i]);
    }

    // Data blocking for 1S
    for (int i = 0; i < N; i++) {
        sum_prog_1s.push_back(0);
        sum_prog2_1s.push_back(0);
        sum_prog_1s_far.push_back(0);
        sum_prog2_1s_far.push_back(0);
        for (int j = 0; j < i + 1; j++) {
            sum_prog_1s[i] += ave1_1s[j];
            sum_prog2_1s[i] += ave2_1s[j];
            sum_prog_1s_far[i] += ave1_1s_far[j];
            sum_prog2_1s_far[i] += ave2_1s_far[j];
        }
        sum_prog_1s[i] /= (i + 1);
        sum_prog2_1s[i] /= (i + 1);
        error_prog_1s.push_back(error(sum_prog_1s, sum_prog2_1s, i));
        sum_prog_1s_far[i] /= (i + 1);
        sum_prog2_1s_far[i] /= (i + 1);
        error_prog_1s_far.push_back(error(sum_prog_1s_far, sum_prog2_1s_far, i));
        avg_output_1s << (i + 1) * L << "\t" << sum_prog_1s[i] << "\t" << error_prog_1s[i] << "\t" << endl;
        avg_output_1s_far << (i + 1) * L << "\t" << sum_prog_1s_far[i] << "\t" << error_prog_1s_far[i] << "\t" << endl;
    }

   /************************************************************************************************/
   // 2P Orbital
   /************************************************************************************************/
    for (int i = 0; i < N; i++) {
        double sum_2p = 0;
        double sum_2p_far = 0;

        for (int j = 0; j < L; j++) {

            double x_new_2p, y_new_2p, z_new_2p, x_new_2p_far, y_new_2p_far, z_new_2p_far;

            if (distribution_choice == 1) {
                x_new_2p = x_2p + rnd.Rannyu(-step[2], step[2]);
                y_new_2p = y_2p + rnd.Rannyu(-step[2], step[2]);
                z_new_2p = z_2p + rnd.Rannyu(-step[2], step[2]);

                x_new_2p_far = x_2p_far + rnd.Rannyu(-step[2], step[2]);
                y_new_2p_far = y_2p_far + rnd.Rannyu(-step[2], step[2]);
                z_new_2p_far = z_2p_far + rnd.Rannyu(-step[2], step[2]);
            }

            else if (distribution_choice == 2) {
                x_new_2p = x_2p + rnd.Gauss(0, step[3]);
                y_new_2p = y_2p + rnd.Gauss(0, step[3]);
                z_new_2p = z_2p + rnd.Gauss(0, step[3]);

                x_new_2p_far = x_2p_far + rnd.Gauss(0, step[3]);
                y_new_2p_far = y_2p_far + rnd.Gauss(0, step[3]);
                z_new_2p_far = z_2p_far + rnd.Gauss(0, step[3]);
            }

            double acceptance_2p = acceptance_ratio(pow(wf_2p(x_2p, y_2p, z_2p), 2), pow(wf_2p(x_new_2p, y_new_2p, z_new_2p), 2));
            double acceptance_2p_far = acceptance_ratio(pow(wf_2p(x_2p_far, y_2p_far, z_2p_far), 2), pow(wf_2p(x_new_2p_far, y_new_2p_far, z_new_2p_far), 2));

            if (rnd.Rannyu() <= acceptance_2p) {
               x_2p = x_new_2p; y_2p = y_new_2p; z_2p = z_new_2p;
               count_accept_2p++;
            }
            if (rnd.Rannyu() <= acceptance_2p_far) {
                x_2p_far = x_new_2p_far; y_2p_far = y_new_2p_far; z_2p_far = z_new_2p_far;
                count_accept_far_2p++;
            }

            pos_output_2p << x_2p << "\t" << y_2p << "\t" << z_2p << endl;
            pos_output_2p_far << x_2p_far << "\t" << y_2p_far << "\t" << z_2p_far << endl;

            sum_2p += sqrt(x_2p * x_2p + y_2p * y_2p + z_2p * z_2p);
            sum_2p_far += sqrt(x_2p_far * x_2p_far + y_2p_far * y_2p_far + z_2p_far * z_2p_far);
        }

        ave1_2p.push_back(sum_2p / L);
        ave2_2p.push_back(ave1_2p[i] * ave1_2p[i]);
        ave1_2p_far.push_back(sum_2p_far / L);
        ave2_2p_far.push_back(ave1_2p_far[i] * ave1_2p_far[i]);
    }

    // Data blocking for 2P
    for (int i = 0; i < N; i++) {
        sum_prog_2p.push_back(0);
        sum_prog2_2p.push_back(0);
        sum_prog_2p_far.push_back(0);
        sum_prog2_2p_far.push_back(0);
        for (int j = 0; j < i + 1; j++) {
            sum_prog_2p[i] += ave1_2p[j];
            sum_prog2_2p[i] += ave2_2p[j];
            sum_prog_2p_far[i] += ave1_2p_far[j];
            sum_prog2_2p_far[i] += ave2_2p_far[j];
        }
        sum_prog_2p[i] /= (i + 1);
        sum_prog2_2p[i] /= (i + 1);
        error_prog_2p.push_back(error(sum_prog_2p, sum_prog2_2p, i));
        sum_prog_2p_far[i] /= (i + 1);
        sum_prog2_2p_far[i] /= (i + 1);
        error_prog_2p_far.push_back(error(sum_prog_2p_far, sum_prog2_2p_far, i));
        avg_output_2p << (i + 1) * L << "\t" << sum_prog_2p[i] << "\t" << error_prog_2p[i] << "\t" << endl;
        avg_output_2p_far << (i + 1) * L << "\t" << sum_prog_2p_far[i] << "\t" << error_prog_2p_far[i] << "\t" << endl;
    }

    int last_block = N - 1;
    cout << "=================================================================================\n";
    
    cout << "Orbital 1s (starting from ("<< x_1s_start << ", " << y_1s_start << ", " << z_1s_start << ")):\n";
    cout << " <r> = " << sum_prog_1s[last_block] << " ± " << error_prog_1s[last_block] << " a₀" << endl;
    cout << "Acceptance rate 1s: " << (double)count_accept_1s / M << endl;

    cout << "\n\nOrbital 1s (starting from ("<< x_1s_far_start << ", " << y_1s_far_start << ", " << z_1s_far_start << ")):\n";
    cout << " <r> = " << sum_prog_1s_far[last_block] << " ± " << error_prog_1s_far[last_block] << " a₀" << endl;
    cout << "Acceptance rate 1s: " << (double)count_accept_far_1s / M << endl;

    cout << "\n\nOrbital 2p (starting from ("<< x_2p_start << ", " << y_2p_start << ", " << z_2p_start << ")):\n";
    cout << " <r> = " << sum_prog_2p[last_block] << " ± " << error_prog_2p[last_block] << " a₀" << endl;
    cout << "Acceptance rate 2p: " << (double)count_accept_2p / M << endl;

    cout << "\n\nOrbital 2p (starting from ("<< x_2p_far_start << ", " << y_2p_far_start  << ", " << z_2p_far_start << ")):\n";
    cout << " <r> = " << sum_prog_2p_far[last_block] << " ± " << error_prog_2p_far[last_block] << " a₀" << endl;
    cout << "Acceptance rate 2p: " << (double)count_accept_far_2p / M << endl;

    cout << "=================================================================================\n";
    cout << "The program has finished successfully." << endl;

    rnd.SaveSeed();
    return 0;
}
