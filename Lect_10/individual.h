#ifndef __INDIVIDUAL_H__
#define __INDIVIDUAL_H__

#include <armadillo>
#include "random.h"

using namespace std;
using namespace arma;

class Individual {
private:
    uvec _x;
    double _fitness;
    int _ncities;
    mat _distance_matrix;

public:
    Individual(int ncities, const mat &dist_matrix);
    
    // Core Functions
    void initialize(Random& rnd);
    void check();
    void compute_fitness();
    void swap(int i, int j);

    // Mutations
    void pair_permutation(Random& rnd);
    void shift(Random& rnd);
    void contiguous_permutation(Random& rnd);
    void inversion(Random& rnd);

    // Getters and Setters
    double get_fitness() const { return _fitness; }
    uvec get_x() const { return _x; }
    void set_x(const uvec &x);
};

#endif