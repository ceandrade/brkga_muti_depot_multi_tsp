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

#include "imdmtsp.h"

//extern "C++" string str(Instance &instance, FullGraph::Arc a);

void Solution::set_subcycles(const vector<vector<FullGraph::Node> > &subcycles)
{
	FullGraph &graph = instance->graph;
	int K = instance->K;

	assert(int(subcycles.size())==K+1);
	assert(int(subcycles[K].size())==K);

	for(FullGraph::NodeIt v(graph); v != INVALID; ++v) {
	    cycle_with_node[v] = -1;
	}

	this->subcycles = subcycles;

	cycles.clear();
	cycles = vector<StaticPath<FullGraph> >(K + 1);

	value = 0;

	for (int i = 0; i <= K; i++)
	{
		SimplePath<FullGraph> path;
		const vector<FullGraph::Node>& cycle = subcycles[i];
		int cycle_size = int(cycle.size());

		if(cycle_size > 1)
            for (int j = 0; j < cycle_size; j++)
            {
                if (i == K)
                    cycle_with_node[cycle[j]] += K;
                else
                    cycle_with_node[cycle[j]] = i;
                FullGraph::Arc arc = graph.arc(cycle[j],
                        cycle[(j + 1) % cycle_size]);

                value += instance->dist[arc];
                cycle_with_edge[arc] = i;
                path.addBack(arc);
            }

		cycles[i] = path;
		//cycles.push_back(StaticPath<FullGraph>(path));
		//path.clear();
	}
}

Solution::Solution(Instance *_instance, ConcordeWrapper *_concorde):
		value(-1), instance(_instance), cycle_with_node(instance->graph), cycle_with_edge(
				instance->graph, -1), cycles(), subcycles(), concorde(_concorde)
{
}

Solution::Solution(Instance *_instance,
        vector<vector<FullGraph::Node> > &subcycles):
        value(-1), instance(_instance), cycle_with_node(instance->graph), cycle_with_edge(
                instance->graph, -1), cycles(), subcycles(), concorde(NULL)
{
    this->set_subcycles(subcycles);
}


Solution::Solution(Instance *instance,
		vector<vector<FullGraph::Node> > &subcycles, ConcordeWrapper *_concorde):
		value(-1), instance(instance), cycle_with_node(instance->graph), cycle_with_edge(
				instance->graph, -1), cycles(), subcycles(), concorde(_concorde)
{
	this->set_subcycles(subcycles);
}

Solution::~Solution() {
}

//Solution::Solution(Instance &instance, FullGraph::EdgeMap<GRBVar> &x,
//		FullGraph::NodeMap<GRBVar> &y, FullGraph::EdgeMap<GRBVar> &z) :
//		value(-1), instance(instance), cycle_with_node(instance.graph), cycle_with_edge(
//				instance.graph, -1), cycles(), subcycles(instance.K + 1)
//{
//	FullGraph &graph = instance.graph;
//	int K = instance.K;
//	cerr << "13 " << y[graph(13)].get(GRB_DoubleAttr_X) << endl;
//	cerr << "12 " << y[graph(12)].get(GRB_DoubleAttr_X) << endl;
//	cerr << "6 " << y[graph(6)].get(GRB_DoubleAttr_X) << endl;
//	//exit(0);
//	FullGraph::NodeMap<bool> added(graph, false);
//
//	FullGraph::Node first_depot = graph(0);
//	//find some depot node
//	for (FullGraph::NodeIt v(graph); v != INVALID; ++v)
//	{
//		if (y[v].get(GRB_DoubleAttr_X) >= 0.999)
//		{
//			first_depot = v;
//			break;
//		}
//	}
//
//	//find other depot nodes (in order)
//	for (FullGraph::Node v = first_depot; !added[v];)
//	{
//		added[v] = true;
//		subcycles[K].push_back(v);
//		//find some new neighbor depot
//		for (FullGraph::IncEdgeIt e(graph, v); e != INVALID; ++e)
//		{
//			FullGraph::Node other = graph.oppositeNode(v, e);
//			if (!added[other] && x[e].get(GRB_DoubleAttr_X) >= 0.999)
//			{
//				v = other;
//				break;
//			}
//		}
//	}
//
//	for (int k = 0; k < K; k++)
//	{
//		FullGraph::Node depot_k = subcycles[K][k];
//		added[depot_k] = false;
//		for (FullGraph::Node v = depot_k; !added[v];)
//		{
//			cerr << graph.id(v) << " ";
//
//			added[v] = true;
//			subcycles[k].push_back(v);
//			//find some new neighbor node
//			for (FullGraph::IncEdgeIt e(graph, v); e != INVALID; ++e)
//			{
//				FullGraph::Node other = graph.oppositeNode(v, e);
//				if (!added[other] && z[e].get(GRB_DoubleAttr_X) >= 0.999)
//				{
//					v = other;
//					break;
//				}
//			}
//		}
//
//		set_subcycles(subcycles);
//	}
//}

void Solution::draw(string filename) const
{
	FullGraph::NodeMap<Color> node_colors(instance->graph, WHITE);
	FullGraph::EdgeMap<double> widths(instance->graph, 1);
	FullGraph::EdgeMap<Color> edge_colors(instance->graph);

	typedef FullGraph::EdgeMap<bool> EdgeFilter;
	EdgeFilter edge_filter(instance->graph, false);

	for (StaticPath<FullGraph>::ArcIt a(cycles[instance->K]); a != INVALID; ++a)
	{
		FullGraph::Arc arc = a;
		node_colors[instance->graph.source(arc)] = RED;
		edge_colors[arc] = RED;
		edge_filter[arc] = true;
	}

	for (int i = 0; i < instance->K; i++)
		for (StaticPath<FullGraph>::ArcIt a(cycles[i]); a != INVALID; ++a)
		{
			FullGraph::Arc arc = a;
			edge_colors[arc] = GREEN;
			edge_filter[arc] = true;
		}

	FilterEdges<FullGraph, EdgeFilter> new_graph(instance->graph, edge_filter);

	IdMap<FullGraph, FullGraph::Node> id_map(instance->graph);

	ostringstream title;
	title << filename << " | Value: " << value;

	graphToEps(new_graph, filename.c_str()).coords(instance->coords).scale(0.1).nodeScale(
			0.005).nodeColors(node_colors).edgeColors(edge_colors).nodeTexts(
			id_map).nodeTextSize(.5).title(title.str()).run();
}

//string Solution::to_string() const
//{
//	ostringstream oss;
//	int cycles_count = int(cycles.size());
//	for (int i = 0; i < cycles_count; i++)
//	{
//		for (StaticPath<FullGraph>::ArcIt a(cycles[i]); a != INVALID; ++a)
//			oss << str(instance, a) << " ";
//		oss << endl;
//	}
//	return oss.str();
//}

string Solution::check_validity() const
{
	ostringstream message;
	int cycles_count = int(cycles.size()) - 1;

	if (cycles_count != instance->K)
		message << "INVALID: solution has " << cycles_count
				<< " external cycles (should be " << instance->K << ")" << endl;

	if (int(cycles[instance->K].length()) != instance->K)
		message << "INVALID: internal cycle has length"
				<< int(cycles[instance->K].length()) << " (should be "
				<< instance->K << ")" << endl;

	for (int i = 0; i < cycles_count + 1; i++)
	{
		if (!checkPath(instance->graph, cycles[i]))
			message << "INVALID: cycle " << i << " is invalid" << endl;
		//each cycle must start and finish in the same node
		if (pathSource(instance->graph, cycles[i])
				!= pathTarget(instance->graph, cycles[i]))
			message << "INVALID: \"cycle\" " << i
					<< " is not closed (i.e, it is a path)" << endl;
	}

	int max_cycle_size = instance->n - 3 * (instance->K - 1);
	for (int i = 0; i < cycles_count; i++)
	{
		int cycle_size = int(cycles[i].length());
		if (cycle_size < 3 || cycle_size > max_cycle_size)
			message << "INVALID: cycle " << i << " has size " << cycle_size
					<< " (should be in range [3, " << max_cycle_size << "])"
					<< endl;
	}

	FullGraph::NodeMap<int> count(instance->graph, 0);
	for (int i = 0; i < cycles_count; i++)
		for (StaticPath<FullGraph>::ArcIt a(cycles[i]); a != INVALID; ++a)
		{
			FullGraph::Arc arc = a;
			++count[instance->graph.source(arc)];
			++count[instance->graph.target(arc)];
		}
	for (FullGraph::NodeIt v(instance->graph); v != INVALID; ++v)
		if (count[v] != 2)
			message << "INVALID: vertex " << instance->graph.id(v)
					<< " has degree " << count[v] << " (should be 2)" << endl;

	string ret = message.str();
	if (ret == "")
	{
		return "VALID";
	}
	else
	{
		return ret;
	}
}

//void Solution::set_as_mip_start(FullGraph::EdgeMap<GRBVar>& x,
//		FullGraph::NodeMap<GRBVar>& y, FullGraph::EdgeMap<GRBVar>& z,
//		FullGraph::NodeMap<GRBVar*>& t)
//{
//	int K = instance.K;
//	FullGraph &graph = instance.graph;
//
//	for (FullGraph::EdgeIt e(graph); e != INVALID; ++e)
//	{
//		int cycle = cycle_with_edge[e];
//		if (cycle == K) //internal cycle
//		{
//			x[e].set(GRB_DoubleAttr_Start, 1.0);
//			z[e].set(GRB_DoubleAttr_Start, 0.0);
//		}
//		else if (cycle != -1) //external cycle
//		{
//			x[e].set(GRB_DoubleAttr_Start, 0.0);
//			z[e].set(GRB_DoubleAttr_Start, 1.0);
//		}
//	}
//
//	for (FullGraph::NodeIt v(graph); v != INVALID; ++v)
//	{
//		for (int k = 0; k < K; k++)
//			t[v][k].set(GRB_DoubleAttr_Start, 0.0);
//
//		int cycle = cycle_with_node[v];
//		if (cycle >= K) //internal cycle
//		{
//			y[v].set(GRB_DoubleAttr_Start, 1.0);
//			t[v][cycle - K].set(GRB_DoubleAttr_Start, 1.0);
//		}
//		else //external cycle
//		{
//			y[v].set(GRB_DoubleAttr_Start, 0.0);
//			t[v][cycle].set(GRB_DoubleAttr_Start, 1.0);
//		}
//	}
//}

void Solution::lin_kernighan()
{
	vector<vector<FullGraph::Node> > subcycles(instance->K + 1);

	for (int k = 0; k <= instance->K; k++)
	{
		concorde->getLinKernTour(this->subcycles[k], &subcycles[k]);
	}

	set_subcycles(subcycles);
}

void Solution::choose_depots_heuristic()
{
	int K = instance->K;
	FullGraph &graph = instance->graph;
	for (bool improved = true; improved;)
	{
		improved = false;
		for (int k = 0; k < K; k++)
		{
			vector<FullGraph::Node> &subcycle = subcycles[k];
			FullGraph::Node pred_depot = subcycles[K][(k - 1 + K) % K];
			FullGraph::Node succ_depot = subcycles[K][(k + 1) % K];

			DistanceType min_dist_sum = numeric_limits<DistanceType>::max();
			int best_depot = 0;
			int cycle_size = int(subcycle.size());

			for (int i = 0; i < cycle_size; i++)
			{
				FullGraph::Node depot = subcycle[i];

				DistanceType dist_sum = instance->dist[graph.edge(depot,
						pred_depot)]
						+ instance->dist[graph.edge(depot, succ_depot)];

				if (dist_sum < min_dist_sum)
				{
					min_dist_sum = dist_sum;
					best_depot = i;
				}
			}

			if (best_depot != 0)
			{
				//cerr << "improved" << endl;
				improved = true;
				rotate(subcycle.begin(), subcycle.begin() + best_depot,
						subcycle.end());
				subcycles[K][k] = subcycle[0];
			}
		}
	}
	set_subcycles(subcycles);
}

//void Solution::node_exchange_heuristic()
//{
//    int K = instance.K;
//    FullGraph &graph = instance.graph;
//    FullGraph::EdgeMap<DistanceType> &dist = instance.dist;
//
//    bool improved;
//    do
//    {
//        improved = false;
//        for (int k1 = 0; k1 < K; k1++)
//            for (int k2 = k1 + 1; k2 < K; k2++)
//            {
//                vector<FullGraph::Node> &subcycle1 = subcycles[k1], &subcycle2 =
//                        subcycles[k2];
//                int cycle_size1 = int(subcycle1.size());
//                int cycle_size2 = int(subcycle2.size());
//
//                //do not exchange depots
//                for (int i = 1; i < cycle_size1; i++)
//                {
//                    int pre_i = (i == 0) ? (cycle_size1 - 1) : (i - 1);
//                    int suc_i = (i == cycle_size1 - 1) ? 0 : (i + 1);
//                    FullGraph::Node node1 = subcycle1[i];
//                    FullGraph::Node pre_node1 = subcycle1[pre_i];
//                    FullGraph::Node suc_node1 = subcycle1[suc_i];
//
//                    //do not exchange depots
//                    for (int j = 1; j < cycle_size2; j++)
//                    {
//                        int pre_j = (j == 0) ? (cycle_size2 - 1) : (j - 1);
//                        int suc_j = (j == cycle_size2 - 1) ? 0 : (j + 1);
//                        FullGraph::Node node2 = subcycle2[j];
//                        FullGraph::Node pre_node2 = subcycle2[pre_j];
//                        FullGraph::Node suc_node2 = subcycle2[suc_j];
//
//                        //(pre_i, i), (i, suc_i), (pre_j, j), (j, suc_j)
//                        //(pre_i, j), (j, suc_i), (pre_j, i), (i, suc_j)
//
//                        DistanceType dist1 = dist[graph.edge(pre_node1, node1)];
//                        DistanceType dist2 = dist[graph.edge(node1, suc_node1)];
//                        DistanceType dist3 = dist[graph.edge(pre_node2, node2)];
//                        DistanceType dist4 = dist[graph.edge(node2, suc_node2)];
//                        DistanceType dist5 = dist[graph.edge(pre_node1, node2)];
//                        DistanceType dist6 = dist[graph.edge(node2, suc_node1)];
//                        DistanceType dist7 = dist[graph.edge(pre_node2, node1)];
//                        DistanceType dist8 = dist[graph.edge(node1, suc_node2)];
//
//                        DistanceType old_dist = dist1 + dist2 + dist3 + dist4;
//                        DistanceType new_dist = dist5 + dist6 + dist7 + dist8;
//                        /*
//                         cerr << graph.id(node1) << " " << graph.id(node2)
//                         << " dist1=" << dist1 << " dist2=" << dist2
//                         << " dist3=" << dist3 << " dist4=" << dist4
//                         << " dist5=" << dist5 << " dist6=" << dist6
//                         << " dist7=" << dist7 << " dist8=" << dist8
//                         << " old_dist=" << old_dist << ", new_dist="
//                         << new_dist << endl;
//                         */
//                        if (new_dist < old_dist)
//                        {
//                            cerr << graph.id(subcycles[k1][i]) << " "
//                                    << graph.id(subcycles[k2][j]) << "    ";
//                            swap(subcycle1[i], subcycle2[j]); //CHECK THIS!
//                            cerr << graph.id(subcycles[k1][i]) << " "
//                                    << graph.id(subcycles[k2][j]) << endl;
//                            //cerr << "improved" << endl;
//                            improved = true;
//                            //cerr << new_dist << " para " << old_dist << endl;
//
//                            set_subcycles(subcycles); //maybe unnecessary
//                        }
//                    }
//                }
//            }
//    } while (improved);
//    set_subcycles(subcycles);
//}
//
void Solution::node_moving_heuristic()
{
    int K = instance->K;
    FullGraph &graph = instance->graph;
    FullGraph::EdgeMap<DistanceType> &dist = instance->dist;

    improved:

    //cerr << "Current value: " << value << endl;

    for (int k1 = 0; k1 < K; k1++)
    {
        vector<FullGraph::Node> &subcycle1 = subcycles[k1];
        int cycle_size1 = int(subcycle1.size());

        if (cycle_size1 < 3)
            continue;

        for (int k2 = k1 + 1; k2 < K; k2++)
        {
            vector<FullGraph::Node> &subcycle2 = subcycles[k2];
            int cycle_size2 = int(subcycle2.size());

            //do not move depots
            for (int i = 1; i < cycle_size1; i++)
            {
                int pre_i = (i == 0) ? (cycle_size1 - 1) : (i - 1);
                int suc_i = (i == cycle_size1 - 1) ? 0 : (i + 1);
                FullGraph::Node node1 = subcycle1[i];
                FullGraph::Node pre_node1 = subcycle1[pre_i];
                FullGraph::Node suc_node1 = subcycle1[suc_i];

                //do not move depots
                for (int j = 1; j < cycle_size2; j++)
                {
                    int pre_j = (j == 0) ? (cycle_size2 - 1) : (j - 1);
                    FullGraph::Node node2 = subcycle2[j];
                    FullGraph::Node pre_node2 = subcycle2[pre_j];

                    //(pre_i, i), (i, suc_i), (pred_j, j)
                    //(pre_i, suc_i), (pre_j, i), (i, j)

                    DistanceType dist1 = dist[graph.edge(pre_node1, node1)];
                    DistanceType dist2 = dist[graph.edge(node1, suc_node1)];
                    DistanceType dist3 = dist[graph.edge(pre_node2, node2)];

                    DistanceType dist4 = dist[graph.edge(pre_node1, suc_node1)];
                    DistanceType dist5 = dist[graph.edge(pre_node2, node1)];
                    DistanceType dist6 = dist[graph.edge(node1, node2)];

                    DistanceType old_dist = dist1 + dist2 + dist3;
                    DistanceType new_dist = dist4 + dist5 + dist6;

                    if (new_dist < old_dist)
                    {
                        subcycle2.insert(subcycle2.begin() + j, node1);
                        subcycle1.erase(subcycle1.begin() + i);

                        set_subcycles(subcycles); //maybe unnecessary
                        goto improved;
                    }
                }
            }
        }
    }
}

void Solution::improve()
{
    //cerr << "Improved a solution from " << value;

    double old_value = numeric_limits<double>::max();
    int i = 0;
//    draw(string("ze_") + instance->output_name + "-" + char(0 + '0') + ".eps");

    while((old_value - value) > 0.0001)
    {
        old_value = value;

        lin_kernighan();
        //choose_depots_heuristic();
        //node_exchange_heuristic();
        node_moving_heuristic();

//        draw(string("ze_") + instance->output_name + "-" + char(i + '0') + ".eps");
        ++i;
    }

//  node_exchange_heuristic();
//  set_subcycles(subcycles);
//  draw(instance.output_name + "-" + char(i + '0') + "D.eps");

    //cerr << " to " << value << " (after " << i << " improvements)" << endl;
}

