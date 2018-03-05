/*
 * Population.h
 *
 * Encapsulates a population of chromosomes represented by a vector of doubles. We don't decode
 * nor deal with random numbers here; instead, we provide private support methods to set the
 * fitness of a specific chromosome as well as access methods to each allele. Note that the BRKGA
 * class must have access to such methods and thus is a friend.
 *
 *  Created on : Jun 21, 2010 by rtoso
 *  Last update: Nov 08, 2011 by andrade
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

#ifndef POPULATION_HPP
#define POPULATION_HPP

//#include "population.hpp"

#include <vector>

class Population {
	template< class Decoder, class RNG, template<class T> class Compare >
	friend class BRKGA;

public:
	typedef unsigned Allele;
	typedef std::vector< Allele > Chromosome;

	Population(const Population& other);
	Population(unsigned n, unsigned p);
	~Population();

	unsigned getN() const;	// Size of each chromosome
	unsigned getP() const;	// Size of population

	Allele operator()(unsigned i, unsigned j) const;	// Direct access to allele j of chromosome i

	// These methods REQUIRE fitness to be sorted, and thus a call to sortFitness() beforehand
	// (and I won't implement a lazy scheme to automatically update sortFitness() here).
	double getBestFitness() const;	// Returns the best fitness in this population
	double getFitness(unsigned i) const;	// Returns the fitness of chromosome i
	const Chromosome& getChromosome(unsigned i) const;	// Returns i-th best chromosome

private:
	std::vector< Chromosome > population;		// Population as vectors of prob.
	std::vector< std::pair< double, unsigned > > fitness;	// Fitness (double) of a each chromosome

	template< template<class T> class Compare >
	void sortFitness(); // Sorts 'fitness' by its first parameter by predicate template class.
	void setFitness(unsigned i, double f);				// Sets the fitness of chromosome i
	Chromosome& getChromosome(unsigned i);	// Returns a chromosome

	Allele& operator()(unsigned i, unsigned j);		// Direct access to allele j of chromosome i
	Chromosome& operator()(unsigned i);	// Direct access to chromosome i
};

template< template<class T> class Compare >
void Population::sortFitness() {
    sort(fitness.begin(), fitness.end(), Compare< std::pair< double, unsigned > >());
}

#endif //POPULATION
