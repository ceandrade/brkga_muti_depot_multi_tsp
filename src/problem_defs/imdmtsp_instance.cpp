/******************************************************************************
 * imdmtsp_instance.cpp: Implementation for IMDMTSP Instance class.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2012-2018, Carlos Eduardo de Andrade. All Rights Reserved.
 *
 * Created on : Jun 11, 2012 by andrade
 * Last update: Jun 29, 2012 by andrade
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

#ifdef DEBUG
#include <iostream>
#include <iomanip>
#include <iterator>
#include <cassert>
using namespace std;
#endif

#include "imdmtsp_instance.hpp"
#include <lemon/graph_to_eps.h>
#include <lemon/adaptors.h>
using namespace lemon;

//--------------------------[ Default Constructor ]---------------------------//

IMDMTSP_Instance::IMDMTSP_Instance():
    name("unknown"), graph(), dist(graph), has_coords(false), coords(graph),
    alpha(1.0)
{}

//------------------------------[ Destructor ]--------------------------------//

IMDMTSP_Instance::~IMDMTSP_Instance() {}

//----------------------------[ Draw Instance ]-------------------------------//

void IMDMTSP_Instance::draw() {
    typedef FullGraph::EdgeMap<bool> EdgeFilter;
    EdgeFilter edge_filter(graph, false);
    IdMap<FullGraph, FullGraph::Node> id_map(graph);

    FilterEdges<FullGraph, EdgeFilter> graph_without_edges(graph, edge_filter);

    graphToEps(graph_without_edges, (name + ".eps"))
        .coords(coords)
//        .scale(0.1)
        .nodeScale(0.01)
        .nodeTexts(id_map).nodeTextSize(2)
        .run();
}

