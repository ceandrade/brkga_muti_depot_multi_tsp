/*
 * Population.cpp
 *
 * Encapsulates a population of chromosomes represented by a vector of doubles. We don't decode
 * nor deal with random numbers here; instead, we provide private support methods to set the
 * fitness of a specific chromosome as well as access methods to each allele. Note that the BRKGA
 * class must have access to such methods and thus is a friend.
 *
 *  Created on : Nov 08, 2010 by andrade
 *  Last update: Nov 08, 2011 by andrade
 *      Authors: Rodrigo Franco Toso <rtoso@cs.rutgers.edu>
 *               Carlos Eduardo de Andrade <andrade@ic.unicamp.br>
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

#ifndef POPULATION_H
#define POPULATION_H

#include <algorithm>
#include <exception>
#include <stdexcept>

#include "population.hpp"

Population::Population(const Population& pop) :
		population(pop.population),
		fitness(pop.fitness) {
}

Population::Population(const unsigned n, const unsigned p) :
		population(p, Chromosome(n, 0.0)), fitness(p) {
	if(p == 0) { throw std::range_error("Population size p cannot be zero."); }
	if(n == 0) { throw std::range_error("Chromosome size n cannot be zero."); }
}

Population::~Population() {
}

unsigned Population::getN() const {
	return population[0].size();
}

unsigned Population::getP() const {
	return population.size();
}

double Population::getBestFitness() const {
	return getFitness(0);
}

double Population::getFitness(unsigned i) const {
	return fitness[i].first;
}

const Population::Chromosome& Population::getChromosome(unsigned i) const {
	return population[ fitness[i].second ];
}

Population::Chromosome& Population::getChromosome(unsigned i) {
	return population[ fitness[i].second ];
}

void Population::setFitness(unsigned i, double f) {
	fitness[i].first = f;
	fitness[i].second = i;
}

Population::Allele Population::operator()(unsigned chromosome, unsigned allele) const {
	return population[chromosome][allele];
}

Population::Allele& Population::operator()(unsigned chromosome, unsigned allele) {
	return population[chromosome][allele];
}

Population::Chromosome& Population::operator()(unsigned chromosome) {
	return population[chromosome];
}

#endif
