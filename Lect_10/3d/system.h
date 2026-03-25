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
    // Parametri
    int _ncities;
    int _npop;
    int _ngen;
    int _generation_index;
    double _mut_prob;
    double _p_exp;
    string _type;
    
    // Oggetti Armadillo e Membri
    mat _positions;
    mat _distance_matrix;
    vector<Individual> _individuals;
    Random* _rnd;

    ofstream _out_loss;

    // Metodi privati di supporto
    void initialize_circle();
    void initialize_square();
    void compute_distance_matrix();
    void fitness_sorting();
    Individual rank_selection();
    void crossover(Individual &p1, Individual &p2, Individual &c1, Individual &c2);
    double get_average_best_half();

public:
    // Costruttore e Distruttore
    System();
    ~System();
    
    // Inizializzazione
    void initialize(Random& rnd, const string& type, int ncities, int npop, int ngen, double mut_prob, double p_exp);
    void initialize_from_file(Random& rnd, const string& filename, int ncities, int npop, int ngen, double mut_prob, double p_exp);
    
    // Evoluzione
    void new_generation();
    int get_ngenerations();

    // Metodi per la PARALLELIZZAZIONE (quelli che causavano l'errore)
    // Ritorna il percorso del miglior individuo attuale
    arma::uvec get_best_individual_path() { return _individuals[0].get_x(); }
    // Accetta un individuo da un altro rank
    void accept_immigrant(const arma::uvec &path);

    // Finalizzazione
    void finalize();            // Per esecuzioni seriali
    void finalize(int rank);    // Per esecuzioni parallele (aggiunge il rank al nome file)
};

#endif