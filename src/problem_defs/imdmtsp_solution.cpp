/******************************************************************************
 * imdmtsp_solution.cpp: Implementation for IMDMTSP Solution.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2012-2018, Carlos Eduardo de Andrade. All Rights Reserved.
 *
 * Created on : Jun 21, 2012 by andrade
 * Last update: Jul 31, 2012 by andrade
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

#include <iostream>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <lemon/adaptors.h>
#include <lemon/graph_to_eps.h>

using namespace std;
using namespace lemon;

#include "imdmtsp_solution.hpp"

//--------------------------[ Default Constructor ]---------------------------//

IMDMTSP_Solution::IMDMTSP_Solution(const IMDMTSP_Instance &_instance):
    K(0), cycles(), value(std::numeric_limits< double >::max()), instance(_instance)
{}

//-----------------------------[ Destructor ]---------------------------------//

IMDMTSP_Solution::~IMDMTSP_Solution() {}

//----------------------------[ Draw Instance ]-------------------------------//

void IMDMTSP_Solution::draw(const string &file_name) {
    typedef FullGraph::EdgeMap<bool> EdgeFilter;
    FullGraph::NodeMap<Color> node_colors(instance.graph, WHITE);
    FullGraph::EdgeMap<double> widths(instance.graph, 1);
    FullGraph::EdgeMap<Color> edge_colors(instance.graph);

    EdgeFilter edge_filter(instance.graph, false);
    IdMap<FullGraph, FullGraph::Node> id_map(instance.graph);

    vector< deque< FullGraph::Node > >::iterator cy_it1, cy_it2;
    deque< FullGraph::Node >::iterator nd_it1, nd_it2;
    FullGraph::Edge edge;

    for(cy_it2 = cycles.begin(); cy_it2 != cycles.end(); ) {
        // Color the external cycles
        for(nd_it2 = cy_it2->begin(); nd_it2 != cy_it2->end(); ) {
            nd_it1 = nd_it2++;
            if(nd_it2 != cy_it2->end()) {
                edge = instance.graph.edge(*nd_it1, *nd_it2);
                node_colors[instance.graph.u(edge)] = GREEN;
                node_colors[instance.graph.v(edge)] = GREEN;
                edge_colors[edge] = GREEN;
                edge_filter[edge] = true;
            }
        }

        if(cy_it2->size() > 1) {
            edge = instance.graph.edge(*nd_it1, *(cy_it2->begin()));
            node_colors[instance.graph.u(edge)] = GREEN;
            node_colors[instance.graph.v(edge)] = GREEN;
            edge_colors[edge] = GREEN;
            edge_filter[edge] = true;
        }

        // Color the internal cycle
        cy_it1 = cy_it2++;
        if(cy_it2 != cycles.end()) {
            edge = instance.graph.edge(*(cy_it1->begin()), *(cy_it2->begin()));
            node_colors[instance.graph.u(edge)] = RED;
            node_colors[instance.graph.v(edge)] = RED;
            edge_colors[edge] = RED;
            edge_filter[edge] = true;
        }
    }
    edge = instance.graph.edge(*(cy_it1->begin()), *(cycles.begin()->begin()));
    node_colors[instance.graph.u(edge)] = RED;
    node_colors[instance.graph.v(edge)] = RED;
    edge_colors[edge] = RED;
    edge_filter[edge] = true;

    typedef dim2::Point<double> PointType;
    PointType min_coords = dim2::makePoint(12345678, 12345678);
    PointType max_coords = dim2::makePoint(-12345678, -12345678);

    // Code from Mauro Lopes.
    for (FullGraph::NodeIt v(instance.graph); v != INVALID; ++v) {
        max_coords.x = max(max_coords.x, instance.coords[v].x);
        max_coords.y = max(max_coords.y, instance.coords[v].y);
        min_coords.x = min(min_coords.x, instance.coords[v].x);
        min_coords.y = min(min_coords.y, instance.coords[v].y);
    }
    double scale = max(max_coords.x - min_coords.x, max_coords.y - min_coords.y);

    FilterEdges<const FullGraph, const EdgeFilter> graph_without_edges(instance.graph, edge_filter);

    stringstream title;
    title << file_name << " | Value: " << value;

    string eps_filename = file_name + ".eps";
    string dot_filename = file_name + ".dot";

    graphToEps(graph_without_edges, eps_filename.c_str())
        .coords(instance.coords)
        .nodeScale(0.01)
        .nodeColors(node_colors)
        .edgeColors(edge_colors)
        .nodeTexts(id_map).nodeTextSize(scale/30)
        .title(title.str())
        .run();

    // Write dot file.
    ofstream dot_file(dot_filename.c_str(), ios::out);
    if(!dot_file) {
        stringstream error_msg;
        error_msg << "\n**** Cannot open file \"" << file_name << "\" "
                  << "to write the graph in \"dot\" format";
        throw runtime_error(error_msg.str());
    }

    dot_file << "/* To better visualization use \"dot -Kfdp -s20.0  -Tpdf\" */\n"
             << "graph \"" << title.str() << "\" {\n"
             << "\tnode [shape=circle, fixedsize=\"true\"];\n";

    for(FullGraph::NodeIt v(instance.graph); v != INVALID; ++v) {
        dot_file << "\t" << instance.graph.id(v)
                 << " [pos = \""  << instance.coords[v].x << "," << instance.coords[v].y
                 << "!\"];\n";
    }

    for(FullGraph::EdgeIt e(instance.graph); e != INVALID; ++e) {
        if(edge_filter[e]) {
            dot_file << "\t"
                     << instance.graph.id(instance.graph.u(e))
                     << " -- "
                     << instance.graph.id(instance.graph.v(e))
                     << "";

            if(edge_colors[e].red() == RED.red()) {
                dot_file << " [style=\"dashed\", weight=1, penwidth=3]";
            }
            dot_file << ";\n";
        }
    }

    dot_file << "}\n";
    dot_file.close();
}

//---------------------------[ Print Instance ]-------------------------------//

ostream& operator<<(ostream &os, const IMDMTSP_Solution &solution) {
    for(unsigned i = 0; i < solution.K; ++i) {
        os << "\nCycle " << i << ": ";

        deque< FullGraph::Node >::const_iterator it = solution.cycles[i].begin();
        os << "(" << solution.instance.graph.id(*it) << ") ";

        while(++it != solution.cycles[i].end())
            os << solution.instance.graph.id(*it) << " ";
    }

    os << "\nValue: " << solution.value;

    return os;
}

