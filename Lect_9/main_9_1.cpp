#include <iostream>
#include <fstream>
#include "system.h"
#include "random.h"

using namespace std;

int main() {
    Random rnd;
    int seed[4], p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) Primes >> p1 >> p2;
    else { cerr << "Error: Primes file missing" << endl; return 1; }
    Primes.close();

    ifstream input("seed.in");
    if (input.is_open()) {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed, p1, p2);
        input.close();
    } else { cerr << "Error: seed.in missing" << endl; return 1; }

    // Parameters
    int n_cities = 34;
    int n_pop = 200;
    int n_gen = 600;
    double mut_prob = 0.2;
    double p_exponent = 3.;

    // Circle
    System sys_circle;
    sys_circle.initialize(rnd, "circle", n_cities, n_pop, n_gen, mut_prob, p_exponent);
    
    for(int i = 0; i < sys_circle.get_ngenerations(); i++) {
        sys_circle.new_generation();
    }
    sys_circle.finalize();

    // Square
    System sys_square;
    sys_square.initialize(rnd, "square", n_cities, n_pop, n_gen, mut_prob, p_exponent);
    
    for(int i = 0; i < sys_square.get_ngenerations(); i++) {
        sys_square.new_generation();
    }
    sys_square.finalize();

    cout << "All simulations complete." << endl;
    rnd.SaveSeed();
    
    return 0;
}