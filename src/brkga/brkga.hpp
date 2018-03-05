/*
 * BRKGA.h
 *
 * This class encapsulates a Biased Random-key Genetic Algorithm (for minimization problems) with K
 * independent Populations stored in two vectors of Population, current and previous. It supports
 * multi-threading via OpenMP, and implements the following key methods:
 *
 * - BRKGA() constructor: initializes the populations with parameters described below.
 * - evolve() operator: evolve each Population following the BRKGA methodology. This method
 *                      supports OpenMP to evolve up to K independent Populations in parallel.
 *                      Please note that double Decoder::decode(...) MUST be thread-safe.
 *
 * Required hyperparameters:
 * - n: number of genes in each chromosome
 * - p: number of elements in each population
 * - pe: pct of elite items into each population
 * - pm: pct of mutants introduced at each generation into the population
 * - rhoe: probability that an offspring inherits the allele of its elite parent
 *
 * Optional parameters:
 * - K: number of independent Populations
 * - MAX_THREADS: number of threads to perform parallel decoding -- WARNING: Decoder::decode() MUST
 *                be thread-safe!
 *
 * Required templates are:
 * RNG: random number generator that implements the methods below.
 *     - RNG(unsigned long seed) to initialize a new RNG with 'seed'
 *     - double rand() to return a double precision random deviate in range [0,1)
 *     - unsigned long randInt() to return a >=32-bit unsigned random deviate in range [0,2^32-1)
 *     - unsigned long randInt(N) to return a unsigned random deviate in range [0, N] with N < 2^32
 *
 * Decoder: problem-specific decoder that implements any of the decode methods outlined below. When
 *          compiling and linking BRKGA with -fopenmp (i.e., with multithreading support via
 *          OpenMP), the method must be thread-safe.
 *     - double decode(const vector< double >& chromosome) const, if you don't want to change
 *       chromosomes inside the framework, or
 *     - double decode(vector< double >& chromosome) const, if you'd like to update a chromosome
 *
 *  Created on : Jun 22, 2010 by rtoso
 *  Last update: Jul 03, 2012 by andrade
 *      Authors: Rodrigo Franco Toso <rtoso@cs.rutgers.edu>
 *      Collaborator: Carlos Eduardo de Andrade <andrade@ic.unicamp.br>
 *
 * The MIT License (MIT)
 *
 * Copyright (c) 2018
 * Rodrigo Franco Toso (rfrancotoso@gmail.com) and
 * Mauricio G.C. Resende
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef BRKGA_H
#define BRKGA_H

#include <omp.h>
#include <algorithm>
#include <exception>
#include <stdexcept>
#include <limits>
#include "population.hpp"

template< class Decoder, class RNG, template<class T> class Compare >
class BRKGA {
public:
	/*
	 * Default constructor
	 * Required hyperparameters:
	 * - n: number of genes in each chromosome
	 * - p: number of elements in each population
	 * - pe: pct of elite items into each population
	 * - pm: pct of mutants introduced at each generation into the population
	 * - rhoe: probability that an offspring inherits the allele of its elite parent
	 *
	 * Optional parameters:
	 * - K: number of independent Populations
	 * - MAX_THREADS: number of threads to perform parallel decoding
	 *                WARNING: Decoder::decode() MUST be thread-safe; safe if implemented as
	 *                + double Decoder::decode(Population::Chromosome& chromosome) const
	 */
	BRKGA(unsigned n, unsigned p, double pe, double pm, double rhoe,
		  Decoder& refDecoder, RNG& refRNG, unsigned K = 1, unsigned MAX_THREADS = 1,
            unsigned cut_point = 0,
            unsigned left_lb = numeric_limits< Population::Allele >::min(),
            unsigned left_ub = numeric_limits< Population::Allele >::max(),
            unsigned right_lb = numeric_limits< Population::Allele >::min(),
            unsigned right_ub = numeric_limits< Population::Allele >::max());

	/**
	 * Destructor
	 */
	~BRKGA();

    /**
     * Set individuals to initial population (only one population in case of multiple ones).
     * @param chromosomes a set of individuals described as double vectors
     *        between 0 and 1.
     */
    void setInitialPopulation(const std::vector< Population::Chromosome >& chromosomes);


    /**
     * Initialize the populations and others parameters of the algorithm. If a
     * initial population is supplied, this method completes the remain
     * individuals, if they not exist.
     * @param true_init indicates a true initialization or other from reset phase.
     */
	void initialize(bool true_init = true);

	/**
	 * Resets all populations with brand new keys
	 * @param partial_reset if true indicates partial reset considering the initial population.
	 */
	void reset(bool partial_reset = false);

	/**
	 * Evolve the current populations following the guidelines of BRKGAs
	 * @param generations number of generations (must be even and nonzero)
	 */
	void evolve(unsigned generations = 1);

	/**
	 * Exchange elite-solutions between the populations
	 * @param M number of elite chromosomes to select from each population
	 */
	void exchangeElite(unsigned M);

	/**
	 * Returns the current population
	 */
	const Population& getCurrentPopulation(unsigned k = 0) const;

	/**
	 * Returns the chromosome with best fitness so far among all populations
	 */
	const Population::Chromosome& getBestChromosome() const;

	/**
	 * Returns the best fitness found so far among all populations
	 */
	double getBestFitness() const;

private:
	// Hyperparameters:
	const unsigned n;	// number of genes in the chromosome
	const unsigned p;	// number of elements in the population
	const unsigned pe;	// number of elite items in the population
	const unsigned pm;	// number of mutants introduced at each generation into the population
	const double rhoe;	// probability that an offspring inherits the allele of its elite parent

	// Templates:
	RNG& refRNG;			// reference to the random number generator
	Decoder& refDecoder;	// reference to the problem-dependent Decoder

	// Parallel populations parameters:
	const unsigned K;				// number of independent parallel populations
	const unsigned MAX_THREADS;		// number of threads for parallel decoding

	// Parameters to limit the allele generation
	const unsigned cut_point;   // Separate the chromosome in two sections
	const unsigned left_lb;     // Mininum allele value on left section
	const unsigned left_ub;     // Maximum allele value on left section
	const unsigned right_lb;    // Mininum allele value on right section
	const unsigned right_ub;    // Maximum allele value on right section

	// Data:
	std::vector< Population* > previous;	// previous populations
	std::vector< Population* > current;		// current populations
	bool initialPopulation;                 // indicate if a initial population is set
	bool initialized;                       // indicate if the algorithm was proper initialized
	bool reset_phase;                       // indicate if the algorithm have been reset

    /**
     * Binary function object (functor) used to compare the chromosome fitness.
     * It's important to supply a Strict Weak Ordering Function like
     * less<>, greater<> and other compatible functions.
     */
    const Compare< double > betterThan;

	// Local operations:
	void evolution(Population& curr, Population& next);
	bool isRepeated(const Population::Chromosome& chrA, const Population::Chromosome& chrB) const;
};

template< class Decoder, class RNG, template<class T> class Compare >
BRKGA< Decoder, RNG, Compare >::BRKGA(unsigned _n, unsigned _p, double _pe, double _pm, double _rhoe,
		Decoder& decoder, RNG& rng, unsigned _K, unsigned MAX, unsigned _cut_point, unsigned _left_lb,
		unsigned _left_ub, unsigned _right_lb, unsigned _right_ub):

		n(_n), p(_p),
		pe(unsigned(_pe * p)), pm(unsigned(_pm * p)), rhoe(_rhoe),
		refRNG(rng), refDecoder(decoder), K(_K), MAX_THREADS(MAX),
        cut_point(_cut_point), left_lb(_left_lb), left_ub(_left_ub),
        right_lb(_right_lb), right_ub(_right_ub), previous(K, 0),
		current(K, 0), initialPopulation(false), initialized(false),
		reset_phase(false), betterThan()
{
	// Error check:
	using std::range_error;
	if(n == 0) { throw range_error("Chromosome size equals zero."); }
	if(p == 0) { throw range_error("Population size equals zero."); }
	if(pe == 0) { throw range_error("Elite-set size equals zero."); }
	if(pe > p) { throw range_error("Elite-set size greater than population size (pe > p)."); }
	if(pm > p) { throw range_error("Mutant-set size (pm) greater than population size (p)."); }
	if(pe + pm > p) { throw range_error("elite + mutant sets greater than population size (p)."); }
	if(K == 0) { throw range_error("Number of parallel populations cannot be zero."); }
}

template< class Decoder, class RNG, template<class T> class Compare >
BRKGA< Decoder, RNG, Compare >::~BRKGA() {
	for(unsigned i = 0; i < K; ++i) { delete current[i]; delete previous[i]; }
}

template< class Decoder, class RNG, template<class T> class Compare >
const Population& BRKGA< Decoder, RNG, Compare >::getCurrentPopulation(unsigned k) const {
	return (*current[k]);
}

template< class Decoder, class RNG, template<class T> class Compare >
double BRKGA< Decoder, RNG, Compare >::getBestFitness() const {
	double best = current[0]->fitness[0].first;
	for(unsigned i = 1; i < K; ++i) {
		if(betterThan(current[i]->fitness[0].first, best)) { best = current[i]->fitness[0].first; }
	}

	return best;
}

template< class Decoder, class RNG, template<class T> class Compare >
const Population::Chromosome& BRKGA< Decoder, RNG, Compare >::getBestChromosome() const {
	unsigned bestK = 0;
	for(unsigned i = 1; i < K; ++i)
		if(betterThan(current[i]->getBestFitness(), current[bestK]->getBestFitness()) ) { bestK = i; }
	return current[bestK]->getChromosome(0);	// The top one :-)
}

template< class Decoder, class RNG, template<class T> class Compare >
void BRKGA< Decoder, RNG, Compare >::reset(bool partial_reset) {
    if(!initialized) {
        throw std::runtime_error("The algorithm hasn't been initialized. Don't forget to call initialize() method");
    }
    reset_phase = true;
    initialize(partial_reset);
}

template< class Decoder, class RNG, template<class T> class Compare >
void BRKGA< Decoder, RNG, Compare >::evolve(unsigned generations) {
    if(!initialized) {
        throw std::runtime_error("The algorithm hasn't been initialized. Don't forget to call initialize() method");
    }

    if(generations == 0) { throw std::range_error("Cannot evolve for 0 generations."); }

	for(unsigned i = 0; i < generations; ++i) {
		for(unsigned j = 0; j < K; ++j) {
			evolution(*current[j], *previous[j]);	// First evolve the population (curr, next)
			std::swap(current[j], previous[j]);		// Update (prev = curr; curr = prev == next)
		}
	}
}

template< class Decoder, class RNG, template<class T> class Compare >
void BRKGA< Decoder, RNG, Compare >::exchangeElite(unsigned M) {
	if(M == 0 || M >= p) { throw std::range_error("M cannot be zero or >= p."); }

	#ifdef _OPENMP
		#pragma omp parallel for num_threads(MAX_THREADS)
	#endif
	for(int i = 0; i < int(K); ++i) {
		// Population i will receive some elite members from each Population j below:
		unsigned dest = p - 1;	// Last chromosome of i (will be overwritten below)
		for(unsigned j = 0; j < K; ++j) {
			if(j == unsigned(i)) { continue; }

			// Copy the M best of Population j into Population i:
			for(unsigned m = 0; m < M; ++m) {
				// Copy the m-th best of Population j into the 'dest'-th position of Population i:
				const Population::Chromosome& bestOfJ = current[j]->getChromosome(m);

				std::copy(bestOfJ.begin(), bestOfJ.end(), current[i]->getChromosome(dest).begin());
				current[i]->fitness[dest].first = current[j]->fitness[m].first;

				--dest;
			}
		}
	}

	// PARALLEL REGION: re-sort each population since they were modified
	#ifdef _OPENMP
		#pragma omp parallel for num_threads(MAX_THREADS)
	#endif
	for(int i = 0; i < int(K); ++i) { current[i]->sortFitness< Compare >(); }
}

template< class Decoder, class RNG, template<class T> class Compare >
void BRKGA< Decoder, RNG, Compare >::setInitialPopulation(const std::vector< Population::Chromosome >& chromosomes) {
//    if(initialPopulation) {
//        throw std::runtime_error("You cannot set initial population twice!");
//    }

    if(initialPopulation)
        delete current[0];

    current[0] = new Population(n, chromosomes.size());
    unsigned i = 0;

    for(vector< Population::Chromosome >::const_iterator it_chrom = chromosomes.begin();
        it_chrom != chromosomes.end(); ++it_chrom, ++i) {

        if(it_chrom->size() != n) {
            throw std::runtime_error("Error on setting initial population: number of genes isn't equal!");
        }
        std::copy(it_chrom->begin(), it_chrom->end(), current[0]->population[i].begin());
    }
    initialPopulation = true;
}

template< class Decoder, class RNG, template<class T> class Compare >
void BRKGA< Decoder, RNG, Compare >::initialize(bool true_init) {
    unsigned start = 0;

    // Verify the initial population and complete or prune it!
    if(initialPopulation && true_init) {
        if(current[0]->population.size() < p) {
            Population* pop = current[0];
            Population::Chromosome chromosome(n);
            unsigned last_chromosome = pop->population.size();

            pop->population.resize(p);
            pop->fitness.resize(p);

            for(; last_chromosome < p; ++last_chromosome) {
                for(unsigned k = 0; k < n; ++k) {
                    //TODO: fix this.
                    //chromosome[k] = refRNG.rand();         // for doubles
                    //chromosome[k] = refRNG.randInt(n, unsigned(-(n + 1)));   // for ints
                    if(k < cut_point)
                        chromosome[k] = refRNG.randInt(left_lb, left_ub);
                    else
                        chromosome[k] = refRNG.randInt(right_lb, right_ub);
                }
                pop->population[last_chromosome] = chromosome;
            }
        }
        // Prune some additional chromosomes
        else if(current[0]->population.size() > p) {
            current[0]->population.resize(p);
            current[0]->fitness.resize(p);
        }
        start = 1;
    }

    // Initialize each chromosome of the current population
    for(; start < K; ++start) {
        // Allocate:
        if(!reset_phase)
            current[start] = new Population(n, p);

        for(unsigned j = 0; j < p; ++j) {
            for(unsigned k = 0; k < n; ++k) {
                //TODO: fix this.
                //(*current[start])(j, k) = refRNG.randInt(n, unsigned(-(n + 1)));   //refRNG.rand();
                if(k < cut_point)
                    (*current[start])(j, k) = refRNG.randInt(left_lb, left_ub);
                else
                    (*current[start])(j, k) = refRNG.randInt(right_lb, right_ub);
            }
        }
    }

    // Initialize and decode each chromosome of the current population, then copy to previous:
    for(unsigned i = 0; i < K; ++i) {
        // Decode:
        #ifdef _OPENMP
            #pragma omp parallel for num_threads(MAX_THREADS)
        #endif
        for(int j = 0; j < int(p); ++j) {
            current[i]->setFitness(j, refDecoder.decode((*current[i])(j)) );
        }

        // Sort:
        current[i]->sortFitness< Compare >();

        // Then just copy to previous:
        if(!reset_phase)
            previous[i] = new Population(*current[i]);

        if(reset_phase && true_init) {
            delete previous[i];
            previous[i] = new Population(*current[i]);
        }
    }

    initialized = true;
}

template< class Decoder, class RNG, template<class T> class Compare >
inline void BRKGA< Decoder, RNG, Compare >::evolution(Population& curr, Population& next) {
	// We now will set every chromosome of 'current', iterating with 'i':
	unsigned i = 0;	// Iterate chromosome by chromosome
	unsigned j = 0;	// Iterate allele by allele

	// 2. The 'pe' best chromosomes are maintained, so we just copy these into 'current':
	while(i < pe) {
		for(j = 0 ; j < n; ++j) { next(i,j) = curr(curr.fitness[i].second, j); }

		next.fitness[i].first = curr.fitness[i].first;
		next.fitness[i].second = i;
		++i;
	}

	// 3. We'll mate 'p - pe - pm' pairs; initially, i = pe, so we need to iterate until i < p - pm:
	while(i < p - pm) {
		// Select an elite parent:
		const unsigned eliteParent = (refRNG.randInt(pe - 1));

		// Select a non-elite parent:
		const unsigned noneliteParent = pe + (refRNG.randInt(p - pe - 1));

		// Mate:
		for(j = 0; j < n; ++j) {
			const unsigned sourceParent = ((refRNG.rand() < rhoe) ? eliteParent : noneliteParent);

			next(i, j) = curr(curr.fitness[sourceParent].second, j);

			//next(i, j) = (refRNG.rand() < rhoe) ? curr(curr.fitness[eliteParent].second, j) :
			//		                              curr(curr.fitness[noneliteParent].second, j);
		}

		++i;
	}

	// We'll introduce 'pm' mutants:
	while(i < p) {
	    //TODO: fix this.
	    for(j = 0; j < n; ++j) {
	        //next(i, j) = refRNG.randInt(n, unsigned(-(n + 1))); //refRNG.rand();

            if(j < cut_point)
                next(i, j) = refRNG.randInt(left_lb, left_ub);
            else
                next(i, j) = refRNG.randInt(right_lb, right_ub);
	    }
		++i;
	}

	// Time to compute fitness, in parallel:
	#ifdef _OPENMP
		#pragma omp parallel for num_threads(MAX_THREADS)
	#endif
	for(int i = int(pe); i < int(p); ++i) {
		next.setFitness(i, refDecoder.decode(next.population[i]) );
	}

	// Now we must sort 'current' by fitness, since things might have changed:
	next.sortFitness< Compare >();
}

#endif
