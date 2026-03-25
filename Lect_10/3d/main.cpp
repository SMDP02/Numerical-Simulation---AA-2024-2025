#include <iostream>
#include <fstream>
#include <vector>
#include <numeric> 
#include <algorithm>
#include <random> 
#include "mpi.h"
#include "system.h"
#include "random.h"

using namespace std;

int main(int argc, char* argv[]) {
    // 1. MPI INITIALIZATION
    int rank, size;
    MPI_Init(&argc, &argv);              // Initialize the MPI environment
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the ID of the current process (0, 1, 2...)
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get total number of processes running

    // 2. PARALLEL RANDOM SEEDING
    // For a GA to be effective in parallel, every rank MUST explore different paths.
    // We read different lines of the "Primes" file based on the rank ID.
    Random rnd;
    int seed[4], p1, p2;
    ifstream Primes("Primes");
    if (Primes.is_open()) {
        for(int i = 0; i <= rank; i++) Primes >> p1 >> p2; // Skip lines until we reach our rank's unique primes
    } else { 
        if(rank == 0) cerr << "Error: Primes file missing" << endl; 
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    Primes.close();

    ifstream input("seed.in");
    if (input.is_open()) {
        input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
        rnd.SetRandom(seed, p1, p2); // Every rank now has a unique random sequence
        input.close();
    }
    int n_cities = 186;    
    int n_pop = 600;       
    int n_gen = 5000;     
    int n_migr = 300;      
    double mut_prob = 0.1; 
    double p_exponent = 3.5; 

    System sys;
    sys.initialize_from_file(rnd, "concap_filtered_un.csv", n_cities, n_pop, n_gen, mut_prob, p_exponent);

    for(int g = 1; g <= n_gen; g++) {
        sys.new_generation();
        if(g % n_migr == 0) {
            vector<int> destinations(size);
            iota(destinations.begin(), destinations.end(), 0); // Fills 0, 1, 2, ..., size-1

            if(rank == 0) {
                std::random_device rd;
                std::mt19937 g_engine(rd());
                std::shuffle(destinations.begin(), destinations.end(), g_engine); // Randomly shuffle destinations
            }
            
            // Rank 0 broadcasts the shuffle map so everyone knows who to send their best individual to
            MPI_Bcast(destinations.data(), size, MPI_INT, 0, MPI_COMM_WORLD);

            // PREPARING DATA FOR TRANSMISSION
            // We take the best path (uvec) and get its raw memory pointer for MPI
            arma::uvec best_path = sys.get_best_individual_path();
            arma::uword* send_buf = best_path.memptr();
            arma::uword* recv_buf = new arma::uword[n_cities];

            // Determine who I am sending to and who I am receiving from
            int target = destinations[rank];
            int source = -1;
            for(int i=0; i<size; i++) {
                if(destinations[i] == rank) source = i;
            }

            // MPI_Sendrecv: A combined operation that prevents deadlocks.
            // Rank X sends its best path to 'target' and simultaneously receives a path from 'source'.
            MPI_Status stat;
            MPI_Sendrecv(send_buf, n_cities, MPI_UNSIGNED_LONG_LONG, target, 1,
                         recv_buf, n_cities, MPI_UNSIGNED_LONG_LONG, source, 1,
                         MPI_COMM_WORLD, &stat);

            // Convert raw buffer back into an Armadillo vector and inject it
            arma::uvec immigrant_path(n_cities);
            for(int i=0; i<n_cities; i++) immigrant_path(i) = recv_buf[i];
            
            sys.accept_immigrant(immigrant_path);

            delete[] recv_buf; // Clean up temporary heap memory
        }
    }
    
    sys.finalize(rank); 
    MPI_Finalize();
    return 0;
}