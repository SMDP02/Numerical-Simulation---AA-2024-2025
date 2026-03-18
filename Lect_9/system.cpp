#include "system.h" 
#include <iostream>
#include <cmath> 
#include <string> 
#include <algorithm>

using namespace std;
using namespace arma;

// ================================================================================================ 
// POPULATION 
// ================================================================================================

System::System() {} // Default constructor 

System::~System() { // Destructor
    if(_out_loss.is_open()) _out_loss.close();
} 

void System::initialize(Random& rnd, const string& type, int ncities, int npop, int ngen, double mut_prob, double p_exp) { 
    _rnd = &rnd; 
    _type = type;  
    _ncities = ncities;  
    _npop = npop;  
    _ngen = ngen; 
    _mut_prob = mut_prob;  
    _p_exp = p_exp;  
    _generation_index = 0; 

    _positions.set_size(_ncities, 2); // Allocate a matrix (ncities x 2)
    
    if (_type == "circle") initialize_circle(); 
    else if (_type == "square") initialize_square(); 

    // Calculate the distances between all pairs of cities
    compute_distance_matrix(); 

    // Loop to create the initial population of npop individuals
    for(int i=0; i<_npop; i++) { 
        Individual indiv(_ncities, _distance_matrix); 
        indiv.initialize(*_rnd); // Random initialization of individual's genome --> order of the visited cities
        // Add the new individual to the population
        _individuals.push_back(indiv); 
    } 
    fitness_sorting(); // Sort the new population

    // Output file
    _out_loss.open(_type + "_losses.dat");
    cout << "\nStarting Optimization: " << _type << "..." << endl;
} // End of initialize method

void System::initialize_circle() { 
    for(int i=0; i<_ncities; ++i) {
        double theta = _rnd->Rannyu(0, 2 * M_PI); // Random angle theta in [0 and 2*PI)
        _positions(i, 0) = cos(theta); // X coordinate
        _positions(i, 1) = sin(theta); // Y coordinate 
    }
}

void System::initialize_square() { 
    for(int i=0; i<_ncities; ++i) { 
        _positions(i, 0) = _rnd->Rannyu(-1.0, 1.0); // Random X coordinate uniformly distributed in (-1,1)
        _positions(i, 1) = _rnd->Rannyu(-1.0, 1.0); // Random Y coordinate '''''''''''''''''''''''''''''''
    } 
}

void System::compute_distance_matrix() { // Pre-calculation of all distances
    // Allocatation of a square matrix (ncities x ncities)
    _distance_matrix.set_size(_ncities, _ncities); 
    // Loop over the row index corresponding to the first city
    for(int i=0; i<_ncities; ++i) {
        // Loop over the column index for the second city
        for(int j=i+1; j<_ncities; ++j) { // from the (i+1)-th to avoid double calculations
            // Distance between city i and j
            _distance_matrix(i,j) = norm(_positions.row(i) - _positions.row(j)); 
            // Matrix symmetry
            _distance_matrix(j,i) = _distance_matrix(i,j);
        }
    }
} 

// Method to sort the population from the best to the worst
void System::fitness_sorting() { 
    sort(_individuals.begin(), _individuals.end(), [](const Individual &a, const Individual &b) { 
         // Return true if individual 'a' has a lower fitness than 'b'
        return a.get_fitness() < b.get_fitness();
    }); 
}
// Method to select an individual for reproduction 
Individual System::rank_selection() { 
    // Calculate a biased index using a random number raised to the power of _p_exp
    int index = static_cast<int>(_npop * pow(_rnd->Rannyu(), _p_exp)); 
    // Return a copy of the selected individual from the sorted population
    return _individuals[index];
}

void System::crossover(Individual &p1, Individual &p2, Individual &c1, Individual &c2) {
    // Randomly select a crossover cut point (not the first or last element)
    int cut = static_cast<int>(_rnd->Rannyu(1, _ncities - 1)); 
    uvec x1 = p1.get_x();
    uvec x2 = p2.get_x();
    // Allocate uninitialized vectors for the children's genomes
    uvec s1(_ncities), s2(_ncities); 
    // Loop to copy the exact sequence from parents to children up to the cut point
    for(int i=0; i<=cut; i++) { 
        s1(i) = x1(i); // Copy allele from parent 1 to child 1
        s2(i) = x2(i); // ''''''''''''''''''''''  2 ''  ''   2
    } 
    // Insertion indices for the remaining parts of the children's genomes
    int i1 = cut+1, i2 = cut+1; 
    // Iterate through the entire sequence of the partner parent
    for(int i=0; i<_ncities; i++) { 
        // Boolean flags to track if a city is already present in the child's inherited sequence
        bool f1 = false, f2 = false; 
        // Inner loop to check against the already inherited part of the child's genome
            for(int j=0; j<=cut; j++) { 
            if(x2(i) == s1(j)) f1 = true; // Set flag true if the city from parent 2 is already in child 1
            if(x1(i) == s2(j)) f2 = true; // ''''''''''''''''''''''''''''''''''''' 1 ''''''''''''''''''' 2
        } // End of inner loop
        if(!f1) s1(i1++) = x2(i); // If the city is not in child 1, append it and increment the child 1 index
        if(!f2) s2(i2++) = x1(i); // ''''''''''''''''''''''''''' 2 ''''''''''''''''''''''''''''''''   2  ''
    } // End of missing city filling loop
    c1.set_x(s1); // Assign the new genome to child 1
    c2.set_x(s2); // '''''''''''''''''''''''''''''' 2
}

void System::new_generation() {
    // Empty temporary vector
    vector<Individual> next_gen; 
    // Pre-allocate memory for the new population
    next_gen.reserve(_npop); 
    // Loop until the new generation reaches the required population size
    while(next_gen.size() < static_cast<size_t>(_npop)) {
        Individual p1 = rank_selection(); // Select the first parent using the biased rank selection operator
        Individual p2 = rank_selection(); // Select the second parent   ''''''''''''''''''''''''''''''''''
        // Initialize two empty children individuals
        Individual c1(_ncities, _distance_matrix), c2(_ncities, _distance_matrix);
        
        // 70% of performing crossover
        if(_rnd->Rannyu() < 0.7) crossover(p1, p2, c1, c2); 
        // Otherwise --> clone parents' genomes
        else { c1.set_x(p1.get_x()); c2.set_x(p2.get_x()); } 

        if(_rnd->Rannyu() < _mut_prob) c1.pair_permutation(*_rnd);
        if(_rnd->Rannyu() < _mut_prob) c1.shift(*_rnd);
        if(_rnd->Rannyu() < _mut_prob) c1.contiguous_permutation(*_rnd);
        if(_rnd->Rannyu() < _mut_prob) c1.inversion(*_rnd);

        if(_rnd->Rannyu() < _mut_prob) c2.pair_permutation(*_rnd);
        if(_rnd->Rannyu() < _mut_prob) c2.shift(*_rnd);
        if(_rnd->Rannyu() < _mut_prob) c2.contiguous_permutation(*_rnd);
        if(_rnd->Rannyu() < _mut_prob) c2.inversion(*_rnd);

        c1.compute_fitness();
        c2.compute_fitness();
        // Add the second child if the population limit hasn't been reached
        next_gen.push_back(c1);
        if(next_gen.size() < static_cast<size_t>(_npop)) next_gen.push_back(c2);
    }
    _individuals = next_gen; // Replace old population with the new one
    fitness_sorting();
    // Compute the loss of the best individual in this generation
    double best_l = _individuals[0].get_fitness(); 
    // Calculate the average loss of the top 50% of the population to monitor convergence
    double ave_l = get_average_best_half(); 
    _out_loss << _generation_index << " " << best_l << " " << ave_l << endl; 

    if(_generation_index % 100 == 0) cout << _type << " - generation " << _generation_index << " | best: " << best_l << endl; 
    _generation_index++; 
}

double System::get_average_best_half() { 
    double sum = 0;
    for(int i=0; i<_npop/2; i++) sum += _individuals[i].get_fitness(); // Loop through the first half + add their fitness
    return sum / (double(_npop)/2.0); // Divide the sum by half the population size to return the average
} 

int System::get_ngenerations() { // Get method for the total number of generations to run
    return _ngen; 
} 

void System::finalize() { 
    if(_out_loss.is_open()) _out_loss.close(); // Check --> loss tracking file open and close it

    uvec best_path = _individuals[0].get_x(); // Extract the genome of the best individual
    ofstream out_path(_type + "_best_path.dat"); 

    for(int i=0; i<_ncities; ++i) { 
        int idx = best_path(i); 
        out_path << _positions(idx, 0) << " " << _positions(idx, 1) << endl;
    } 
    // Write the first city's coordinates again to close the loop of the path
    out_path << _positions(best_path(0), 0) << " " << _positions(best_path(0), 1) << endl; 
    out_path.close(); 
    
    cout << "Done!" << endl; 
}