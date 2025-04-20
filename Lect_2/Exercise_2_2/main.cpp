#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

double error(const vector<double> &avg, const vector<double> &avg2, int n) {
    if (n == 0) return 0;
    return sqrt((avg2[n] - pow(avg[n], 2)) / n);
}

int main(int argc, char *argv[]) {

    Random rnd;
    int seed[4], p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        Primes >> p1 >> p2;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
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

/***************************************************************************************************/ 
    int RW_choice;
    double a;

    cout << "Choose a type of random walk:" << endl;
    cout << "1. Discrete" << endl;
    cout << "2. Continuous" << endl;
    cout << "Your choice (1 or 2): ";
    cin >> RW_choice;
    if (RW_choice != 1 && RW_choice != 2) {
        cerr << "Invalid choice. Please run the program again and choose either 1 or 2." << endl;
        return 1;
    }
/***************************************************************************************************/ 
    if (RW_choice == 1) {
        cout << "Set the step length (a): ";
        cin >> a;
        cout << "\nCALCULATION ..." << endl;

        int M = 10000;     // Total random walks
        int N = 100;       // Number of blocks
        int L = M / N;     // Walks per block
        int Nsteps = 100;  // Max number of steps

        vector<double> global_avg(Nsteps, 0.), global_avg2(Nsteps, 0.), prog_avg(Nsteps, 0.), prog_avg2(Nsteps, 0.), prog_err(Nsteps, 0.);
        ofstream traj_discrete("trajectory_discrete.dat");

        for (int i = 0; i < N; i++) { 
            vector<double> block_avg(Nsteps, 0.);

            for (int j = 0; j < L; j++) {
                double x = 0., y = 0., z = 0.;

                for (int step = 0; step < Nsteps; step++) {
                    int direction = int(rnd.Rannyu(0, 3)); 
                    int step_sign = (rnd.Rannyu(-1, 1) < 0) ? 1 : -1; 

                    if (direction == 0) x += a * step_sign;
                    else if (direction == 1) y += a * step_sign;
                    else z += a * step_sign;

                    double r2 = x * x + y * y + z * z;
                    block_avg[step] += r2;

                    // Saving only one trajectory 
                    if (j == 0 && i == 0)
                        traj_discrete << step + 1 << "\t" << x << "\t" << y << "\t" << z << endl;
                }
            }

            for (int step = 0; step < Nsteps; step++) {
                block_avg[step] /= L;
                global_avg[step] += block_avg[step];
                global_avg2[step] += block_avg[step] * block_avg[step];
            }
        }
        traj_discrete.close();

        string string_a = to_string(a);
        string filename = "output_RW_discrete_" + string_a + ".data";
        ofstream output_discrete(filename);

        for (int step = 0; step < Nsteps; step++) {
            prog_avg[step] = global_avg[step] / N;
            prog_avg2[step] = global_avg2[step] / N;
            prog_err[step] = error(prog_avg, prog_avg2, step);

            output_discrete << step + 1 << "\t" << sqrt(prog_avg[step]) << "\t" << prog_err[step] / (2 * sqrt(prog_avg[step])) << endl;
        }
        output_discrete.close();
    }
/***************************************************************************************************/ 
    if (RW_choice == 2) {
        cout << "Set the step length (a): ";
        cin >> a;
        cout << "\nCALCULATION ..." << endl;

        int M = 10000;     // Total random walks
        int N = 100;       // Number of blocks
        int L = M / N;     // Walks per block
        int Nsteps = 100;  // Max number of steps

        vector<double> global_avg(Nsteps, 0.), global_avg2(Nsteps, 0.), prog_avg(Nsteps, 0.), prog_avg2(Nsteps, 0.), prog_err(Nsteps, 0.);
        ofstream traj_continuous("trajectory_continuous.dat");

        for (int i = 0; i < N; i++) { 
            vector<double> block_avg(Nsteps, 0.);

            for (int j = 0; j < L; j++) {
                double x = 0., y = 0., z = 0.;

                for (int step = 0; step < Nsteps; step++) {
                    double costheta = rnd.Rannyu(-1, 1);
                    double theta_unif = acos(costheta);

                    double phi = rnd.Rannyu(0, 2 * M_PI);

                    x += a * sin(theta_unif) * cos(phi);
                    y += a * sin(theta_unif) * sin(phi);
                    z += a * cos(theta_unif);

                    double r2 = x * x + y * y + z * z;
                    block_avg[step] += r2;

                    // Saving only one trajectory 
                    if (j == 0 && i == 0)
                        traj_continuous << step + 1 << "\t" << x << "\t" << y << "\t" << z << endl;
                }
            }

            for (int step = 0; step < Nsteps; step++) {
                block_avg[step] /= L;
                global_avg[step] += block_avg[step];
                global_avg2[step] += block_avg[step] * block_avg[step];
            }
        }
        traj_continuous.close();

        string string_a = to_string(a);
        string filename = "output_RW_continuous_" + string_a + ".data";
        ofstream output_continuous(filename);

        for (int step = 0; step < Nsteps; step++) {
            prog_avg[step] = global_avg[step] / N;
            prog_avg2[step] = global_avg2[step] / N;
            prog_err[step] = error(prog_avg, prog_avg2, step);

            output_continuous << step + 1 << "\t" << sqrt(prog_avg[step]) << "\t" << prog_err[step] / (2 * sqrt(prog_avg[step])) << endl;
        }
        output_continuous.close();
    }
    // Uniformly distributed random numbers on the sphere
    // Typical error
    /***************************************************************************************************/ 
    ofstream sphere_wrong("sphere_wrong.dat");
    ofstream sphere_correct("sphere_correct.dat");
    for (int i=0; i<3000; i++) {
        double x = 0., y = 0., z = 0.;
        double x_err = 0., y_err = 0., z_err = 0.;

        double theta = rnd.Rannyu(0, M_PI);
        double costheta = rnd.Rannyu(-1, 1);
        double theta_unif = acos(costheta);

        double phi = rnd.Rannyu(0, 2 * M_PI);

        x += sin(theta_unif) * cos(phi);
        y += sin(theta_unif) * sin(phi);
        z += cos(theta_unif);

        sphere_correct << x << "\t" << y << "\t" << z << endl;
        
        x_err += sin(theta) * cos(phi);
        y_err += sin(theta) * sin(phi);
        z_err += cos(theta);

        sphere_wrong << x_err << "\t" << y_err << "\t" << z_err << endl;

    }
    sphere_wrong.close();
    sphere_correct.close();
    
    cout << "Simulation completed.\n";
    rnd.SaveSeed();
    return 0;
}
