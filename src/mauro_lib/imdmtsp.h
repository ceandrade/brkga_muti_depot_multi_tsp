/******************************************************************************
 * Authors: Mauro C. Lopes <maurolopes@gmail.com>
 *
 * (c) Copyright 2012-2018, Mauro C. Lopes, Carlos Eduardo de Andrade.
 *     All Rights Reserved.
 *
 * Last update: Jul 11, 2014 by maurolopes
 *
 * This code is released under LICENSE.md.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *****************************************************************************/

#ifndef IMDMTSP_H
#define IMDMTSP_H

#include <lemon/core.h>
//#include <lemon/lp_base.h>
//#include <lemon/lp.h>
#include <lemon/full_graph.h>
#include <lemon/graph_to_eps.h>
#include <lemon/adaptors.h>
#include <lemon/gomory_hu.h>
//#include <lemon/time_measure.h>
#include <lemon/path.h>
#include <lemon/arg_parser.h>
#include <lemon/dim2.h>
#include <lemon/random.h>

#include <cstdio>
#include <ctime>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <string>
#include <set>
#include <limits>
#include <algorithm>

/*
 #ifdef DEBUG
 void assert_failed(string file, int line, string function, string condition);
 #define	ASSERT(e) ((e) ? (void)0 : assert_failed(__FILE__, __LINE__, __func__, #e))
 #else
 #define ASSERT(ignore)((void) 0)
 #endif
 */

//#include "gurobi_c++.h"
#include "../include/concorde_wrapper.hpp"

using namespace std;
using namespace lemon;

//types
//typedef unsigned long long int ull;
typedef dim2::Point<double> PointType;

//#ifdef INT_DIST
typedef int DistanceType;
//#else
//typedef double DistanceType;
//#endif

class Instance
{
public:
	FullGraph &graph;
	FullGraph::NodeMap<PointType> &coords;
	FullGraph::EdgeMap<DistanceType> &dist;
	int K, n;
	string name;
	bool has_coords;
	string output_name;

	Instance(FullGraph &_graph, FullGraph::NodeMap<PointType> &_coords,
	         FullGraph::EdgeMap<DistanceType> &_dist);

private:
	void load(string filename);
};

class Solution
{
public:
	double value;
	Instance *instance;
	FullGraph::NodeMap<int> cycle_with_node;
	FullGraph::EdgeMap<int> cycle_with_edge;
	vector<StaticPath<FullGraph> > cycles;
	vector<vector<FullGraph::Node> > subcycles;
	ConcordeWrapper *concorde;

	void draw(string filename) const;
//	string to_string() const;
	string check_validity() const;
//	void set_as_mip_start(FullGraph::EdgeMap<GRBVar>& x,
//			FullGraph::NodeMap<GRBVar>& y, FullGraph::EdgeMap<GRBVar>& z,
//			FullGraph::NodeMap<GRBVar*>& t);
	Solution(Instance *instance, ConcordeWrapper *concorde);
	Solution(Instance *instance, vector<vector<FullGraph::Node> > &subcycles);
	Solution(Instance *instance, vector<vector<FullGraph::Node> > &subcycles, ConcordeWrapper *concorde);
	~Solution();
//	Solution(Instance &instance, FullGraph::EdgeMap<GRBVar> &x,
//			FullGraph::NodeMap<GRBVar> &y, FullGraph::EdgeMap<GRBVar> &z);
	void set_subcycles(const vector<vector<FullGraph::Node> > &subcycles);
	void lin_kernighan();
	void choose_depots_heuristic();
//    void node_exchange_heuristic();
    void node_moving_heuristic();
    void improve();
};

Solution* kmeans_heuristic(Instance *instance);

#endif

