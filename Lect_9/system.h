#ifndef __SYSTEM_H__
#define __SYSTEM_H__

#include <vector>
#include <string>
#include <fstream>
#include <armadillo>
#include <algorithm>
#include "individual.h"
#include "random.h"

using namespace std;
using namespace arma;

class System {
private:
    int _ncities;
    int _npop;
    int _ngen;
    int _generation_index;
    double _mut_prob;
    double _p_exp;
    string _type;
    
    mat _positions;
    mat _distance_matrix;
    vector<Individual> _individuals;
    Random* _rnd;

    ofstream _out_loss;

    void initialize_circle();
    void initialize_square();
    void compute_distance_matrix();
    
    void fitness_sorting();
    Individual rank_selection();
    void crossover(Individual &p1, Individual &p2, Individual &c1, Individual &c2);
    double get_average_best_half();

public:
    System();
    ~System();
    
    void initialize(Random& rnd, const string& type, int ncities, int npop, int ngen, double mut_prob, double p_exp);
    void new_generation();
    int get_ngenerations();
    void finalize();
};

#endif