#include "individual.h"
#include <iostream>

using namespace std;
using namespace arma;

// ================================================================================================
// INDIVIDUAL 
// ================================================================================================

Individual::Individual(int ncities, const mat &dist_matrix) {
    _ncities = ncities;
    _distance_matrix = dist_matrix;
    _x.set_size(_ncities);
}

// Sets initial path and changes it 
void Individual::initialize(Random& rnd) {
    _x.set_size(_ncities); 
    // we exclude the 1st city, which must be fixed
    for (int i = 0; i < _ncities; i++) _x(i) = i;
    for (int i = 0; i < _ncities; i++) {
        int j = static_cast<int>(rnd.Rannyu(1, _ncities));
        int k = static_cast<int>(rnd.Rannyu(1, _ncities));
        this->swap(j, k);
    }
    this->check();
    this->compute_fitness();
}

// Checking that
// 1. The first city is 0 
// 2. All cities appear exactly once
void Individual::check() {
    if (_x(0) != 0) { 
        cerr << "Error: first city is not 0!" << endl;
        exit(1); 
    }
    uvec v_unique = unique(_x); 
    if (v_unique.n_elem != (unsigned int)_ncities) {
        cerr << "Error: duplicate cities found!" << endl;
        exit(1);
    }
}

// Total path length --> our fitness 
void Individual::compute_fitness() {
    double distance = 0.;
    for (int i = 0; i < (_ncities - 1); i++) {
        distance += _distance_matrix(_x(i), _x(i + 1));
    }
    distance += _distance_matrix(_x(_ncities - 1), _x(0)); // Periodic boundary --> we include the comeback to the 1st city
    _fitness = distance;
}

// Useful function that swaps two elements i and j --> see the mutation below
void Individual::swap(int i, int j) {
    unsigned int temp = _x(i);
    _x(i) = _x(j);
    _x(j) = temp;
}

// Mutation 1: Swap two random cities
void Individual::pair_permutation(Random& rnd) {
    // 1 in order to exclude the 1st
    int i = static_cast<int>(rnd.Rannyu(1, _ncities));
    int j = static_cast<int>(rnd.Rannyu(1, _ncities));
    while (i == j) j = static_cast<int>(rnd.Rannyu(1, _ncities));
    this->swap(i, j);
    this->check();
}

// Mutation 2: Shifts a block of m cities by n positions
void Individual::shift(Random& rnd) {
    // m --> the block must to be 
    // 1. of 2+ cities 
    // 2. shorter than the whole sequence
    int m = static_cast<int>(rnd.Rannyu(2, _ncities - 2)); 
    // n
    // 1. has to be 1+ --> unless it will stay in the same position
    // 2. at most (_ncities - m) because is the missing length of the vector
    int n = static_cast<int>(rnd.Rannyu(1, _ncities - m)); 
    // the last city has to be at most (_ncities - n - m) because within the mutation 
    // the last city will be (i + n + m - 1) as i = starting point and the last position of the change will be (n+m)
    int i = static_cast<int>(rnd.Rannyu(1, _ncities - n - m)); 
    uvec y = _x;
    for (int k = 0; k < m; k++) _x(i + n + k) = y(i + k);
    for (int k = 0; k < n; k++) _x(i + k) = y(i + k + m);
    this->check();
}

// Mutation 3: Swap two contiguous blocks of size m
void Individual::contiguous_permutation(Random& rnd) {
    // the minimum length is 2 cities and the maximum must be half of the sequence
    int m = static_cast<int>(rnd.Rannyu(2, _ncities / 2));
    // starting point --> after the 1st fixed city and it must be at most (_ncities - 2 * m) 
    // because 2m is the number of vector positions that we need
    int i = static_cast<int>(rnd.Rannyu(1, _ncities - 2 * m));
    // the second block starts after the first and it must be at most (_ncities - m + 1)
    int j = static_cast<int>(rnd.Rannyu(i + m, _ncities - m + 1));
    for (int k = 0; k < m; k++) {
        this->swap(i + k, j + k);
    }
    this->check();
}

// Mutation 4: Reverse the order of cities within a random block
void Individual::inversion(Random& rnd) {
    int i = static_cast<int>(rnd.Rannyu(1, _ncities - 3)); 
    int l = static_cast<int>(rnd.Rannyu(2, _ncities - i)); 
    int j = i + l - 1;
    for (int k = 0; k < l / 2; k++) {
        this->swap(i + k, j - k);
    }
    this->check();
}

void Individual::set_x(const uvec &x) {
    _x = x;
    this->check();
    this->compute_fitness();
}