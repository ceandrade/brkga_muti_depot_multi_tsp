/******************************************************************************
 * imdmtsp_decoder.hpp: Implementation for K-IMDMTSP Decoder class.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2012-2018, Carlos Eduardo de Andrade. All Rights Reserved.
 *
 * Created on : Jun 06, 2012 by andrade
 * Last update: Sep 19, 2012 by andrade
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
#include <stdexcept>
#include <sstream>
#endif

#include <algorithm>
#include <functional>
#include <limits>
#include <omp.h>
using namespace std;

#include "k_imdmtsp_decoder.hpp"

//------------------[ Default Constructor and Destructor ]--------------------//

K_IMDMTSP_Decoder::K_IMDMTSP_Decoder(const K_IMDMTSP_Instance &_instance,
                                     MTRand &_rng, ConcordeWrapper &_concorde):
        instance(_instance), rng(_rng), concorde(_concorde),
        K(instance.K), max_cycle_size(instance.max_cycle_size),
        sorted_nodes_per_thread(_concorde.getThreadsNum(),
                                PairVector(_instance.graph.nodeNum(), KeyIndexPair()))
{
}

K_IMDMTSP_Decoder::~K_IMDMTSP_Decoder() {}

//--------------------------------[ Decode ]----------------------------------//

double K_IMDMTSP_Decoder::decode(Population::Chromosome& chromosome) {
    #ifdef DEBUG
    cout << setiosflags(ios::fixed) << setprecision(6)
         << "\n------------------------------------------------------\n"
         << "> Decoding chromosome"
         << endl;

    copy(chromosome.begin(), chromosome.end(), ostream_iterator< Population::Allele >(cout, " "));

    Population::Chromosome original_chr(chromosome.size());
    copy(chromosome.begin(), chromosome.end(), original_chr.begin());

    cout << "\n\n> External cycles sizes (Max. cycle size: "
         << max_cycle_size << "):\n";
    copy(chromosome.begin(), chromosome.begin() + instance.K, ostream_iterator< Population::Allele >(cout, " "));
    cout << endl;

    double original_solution_value = 0.0;
    #endif

    //////////////////////////////////////////////////////////////////////
    // First, we need to verify and adjust the external cycles size.
    //////////////////////////////////////////////////////////////////////

    int total = 0;
    for(unsigned i = 0; i < K; ++i) {
//        // If less than 3, put some more nodes in the cycle.
//        if(chromosome[i] < 3) {
//            chromosome[i] = 3;
//        }
//        else
//        // If more then max_size then remove some nodes.
//        if(chromosome[i] > max_cycle_size) {
//            chromosome[i] -= chromosome[i] - max_cycle_size;
//        }
        total += chromosome[i];
    }
//    #ifdef FULLDEBUG
//    cout << "\n> Adjusted cycles sizes (phase 1, total: " << total
//         << ", needed: " << instance.graph.nodeNum() << "):\n";
//    copy(chromosome.begin(), chromosome.begin() + K, ostream_iterator< Population::Allele >(cout, " "));
//    cout << endl;
//    #endif

    // If we have cycles such the sum of number of nodes is less than the
    // the total number of nodes of the graph, we must augment the cycles size.
    if(total != instance.graph.nodeNum()) {
        deque< unsigned > indexes;
        for(unsigned i = 0; i < K; ++i)
            indexes.push_back(i);

        int diff = total - instance.graph.nodeNum();

        // If all the cycles don't cover all nodes, we augment them. If we have
        // cycles such the sum of number of nodes is greater than the total
        // number of nodes of the graph, we must reduce the cycles size.
        int sign = 1;
        if(diff > 0)
            sign = -1;

        while(diff != 0) {
            // Choose a cycle randomly.
            int rn = rng.randInt(indexes.size()-1);
            int index = indexes[rn];

            if((chromosome[index] > 1 && sign < 0) ||
               (chromosome[index] < max_cycle_size && sign > 0)) {

                chromosome[index] += sign;
                diff += sign;
            }
            else // Prevent this cycle to be chosen again.
                indexes.erase(indexes.begin() + rn);
        }
    }

    #ifdef FULLDEBUG
    total = 0;
    for(unsigned i = 0; i < K; ++i)
        total += chromosome[i];
    cout << "\n> Adjusted cycles sizes (phase 2, total: " << total
         << ", needed: " << instance.graph.nodeNum() << "):\n";
    copy(chromosome.begin(), chromosome.begin() + K, ostream_iterator< Population::Allele >(cout, " "));
    cout << endl;
    #endif

    //////////////////////////////////////////////////////////////////////////
    // Now, we proceed with the optimization of the cycles. Each cycle is
    // a sequence of indexes in non-increasing order of key/allele value.
    //////////////////////////////////////////////////////////////////////////
     double solution_value = 0.0;

    // So first, we sort the keys in non-decreasing order by key value;
    PairVector &sorted_nodes = sorted_nodes_per_thread[omp_get_thread_num()];

    for(unsigned i = K; i < chromosome.size(); ++i) {
        KeyIndexPair &index_pair = sorted_nodes[i - K];
        index_pair.index = i;
        index_pair.key = chromosome[i];
    }

    sort(sorted_nodes.begin(), sorted_nodes.end(), greater<KeyIndexPair>());

    #ifdef FULLDEBUG
    cout << "\n> Raw Sorted Chromosome:\n";
    copy(sorted_nodes.begin(), sorted_nodes.end(), ostream_iterator< KeyIndexPair >(cout, " "));
    #endif

    // If we have some identical keys, modify their values.
    for(unsigned i = 0, j = 0, key = numeric_limits<unsigned>::max(); i < sorted_nodes.size(); ++i) {
        if(sorted_nodes[i].key == key)
            sorted_nodes[i].key -= ++j;
        else if(sorted_nodes[i].key <= (key - j)) {
            key = --sorted_nodes[i].key;
            j = 0;
        }
        else if(sorted_nodes[i].key > (key - j)) {
            key -= ++j;
            j = 0;
            sorted_nodes[i].key = key;
        }
    }

    #ifdef FULLDEBUG
    cout << "\n> Adjusted Sorted Chromosome:\n";
    copy(sorted_nodes.begin(), sorted_nodes.end(), ostream_iterator< KeyIndexPair >(cout, " "));
    cout << "\n\n> Cycles:" << endl;
    #endif

    // Vector used in local search.
    vector<FullGraph::Node> initial, opt;
    initial.reserve(max_cycle_size);
    opt.reserve(max_cycle_size);

    // WARNING: The Lemon maps aren't thread safe, unfortunately.
    FullGraph::NodeMap< deque< FullGraph::Node > > *cycles;
    FullGraph::NodeMap< unsigned > *cycle_sizes;

    #ifdef _OPENMP
    #pragma omp critical
    {
    #endif
    cycles = new FullGraph::NodeMap< deque< FullGraph::Node > >(instance.graph);
    cycle_sizes = new FullGraph::NodeMap< unsigned >(instance.graph);
    #ifdef _OPENMP
    }
    #endif

    //////////////////////////////////////////////////////////////////////////
    // Lin Kernighan Neighborhood
    //////////////////////////////////////////////////////////////////////////
    #ifdef FULLDEBUG
    cout << "\n\n|------- Applying Lin Kernighan in the cycles ----------";
    cout.flush();
    #endif

    // Traverse each external cycle and try optimize them at first.
    unsigned offset = K;
    double cycle_value = 0.0;
    for(unsigned cycle_id = 0; cycle_id < K; ++cycle_id) {
        // First, take the depot/gateway node.
        FullGraph::Node node = instance.graph.nodeFromId(sorted_nodes[cycle_id].index - K);
        deque< FullGraph::Node > &ll = (*cycles)[node];

        // First, treated the common cases (cycles >= 3)
        if(chromosome[cycle_id] > 2) {
            const unsigned cycle_size = chromosome[cycle_id] - 1;
            (*cycle_sizes)[node] = cycle_size + 1;

            initial.clear();
            opt.clear();

            // Push the depot.
            initial.push_back(node);

            // Copy the remain nodes
            for(unsigned i = 0; i < cycle_size; ++i)
                initial.push_back(instance.graph.nodeFromId(sorted_nodes[i + offset].index - K));

            // Try improve using the Lin Kernighan heuristic.
            cycle_value = concorde.getLinKernTour(initial, &opt);
            solution_value += cycle_value;

            // Copy the nodes of new cycle except the depot.
            for(vector<FullGraph::Node>::iterator it = opt.begin(); *it != node; ++it)
                ll.push_back(*it);

            for(vector<FullGraph::Node>::reverse_iterator it = opt.rbegin(); *it != node; ++it)
                ll.push_front(*it);

            // Adjust the offset to next cycle.
            offset += cycle_size;
        }
        // Degenerated case: cycle with size 2 (antenna)
        else if(chromosome[cycle_id] == 2) {
            // Put the only the neighbor node in "cycle".
            ll.push_back(instance.graph.nodeFromId(sorted_nodes[offset].index - K));

            // Double the edge value because the round trip.
            solution_value += 2 * instance.dist[instance.graph.edge(node, *ll.begin())];

            (*cycle_sizes)[node] = 2;

            // Adjust the offset to next cycle.
            ++offset;

            #ifdef FULLDEBUG
            cycle_value = 2 * instance.dist[instance.graph.edge(node, *ll.begin())];
            original_solution_value += cycle_value;
            #endif
        }
        // Degenerated case: cycle with size 1 (empty). Only we need is to set the cycle size.
        else if(chromosome[cycle_id] == 1) {
            (*cycle_sizes)[node] = 1;
            #ifdef FULLDEBUG
            cycle_value = 0;
            #endif
        }

        #ifdef DEBUG
        double previous_value = 0.0;
        if(chromosome[cycle_id] > 2) {
            vector<FullGraph::Node>::const_iterator it1, it2;
            for(it2 = initial.begin(); it2 != initial.end(); ) {
                it1 = it2++;
                if(it2 != initial.end())
                    previous_value += instance.dist[instance.graph.edge(*it1,*it2)];
            }
            previous_value += instance.dist[instance.graph.edge(*it1,*(initial.begin()))];
            original_solution_value += previous_value;
            #ifdef FULLDEBUG
            cout << "\n- Init. Cycle[" << cycle_id << "]: " << previous_value << "\n";
            for(vector<FullGraph::Node>::const_iterator it = initial.begin(); it != initial.end(); ++it)
                cout << instance.graph.id(*it) << " ";

            cout << "\n- Opti. Cycle[" << cycle_id << "]: " << cycle_value << "\n"
                 <<  "(" << instance.graph.id(node) << ") ";
            for(deque< FullGraph::Node >::iterator it = ll.begin(); it != ll.end(); ++it)
                cout << instance.graph.id(*it) << " ";
            cout << endl;
            #endif
        }
        #ifdef FULLDEBUG
        else if(chromosome[cycle_id] == 2) {
            cout << "\n- Degenerated Cycle[" << cycle_id << "] (Size 2): " << cycle_value << endl;
        }
        else if(chromosome[cycle_id] == 1) {
            cout << "\n- Degenerated Cycle[" << cycle_id << "] (Size 1): 0" << endl;
        }
        #endif
        #endif
    }

    // Here, we deal with the internal cycle.
    initial.clear();
    opt.clear();
    for(unsigned i = 0; i < K; ++i)
        initial.push_back(instance.graph.nodeFromId(sorted_nodes[i].index - K));

    // Try improve the internal cycle using the Lin Kernighan heuristic.
    double internal_cycle_value = instance.alpha * concorde.getLinKernTour(initial, &opt);
    solution_value += internal_cycle_value;

    #ifdef DEBUG
    double previous_value = 0.0;
    vector<FullGraph::Node>::const_iterator it1, it2;
    for(it2 = initial.begin(); it2 != initial.end(); ) {
        it1 = it2++;
        if(it2 != initial.end())
            previous_value += instance.dist[instance.graph.edge(*it1,*it2)];
    }
    previous_value += instance.dist[instance.graph.edge(*it1,*(initial.begin()))];
    previous_value *= instance.alpha;
    original_solution_value += previous_value;

    #ifdef FULLDEBUG
    cout << "\n- Init. Internal: " << previous_value << "\n";
    for(vector<FullGraph::Node>::const_iterator it = initial.begin(); it != initial.end(); ++it)
        cout << instance.graph.id(*it) << " ";
    cout << "\n- Opti. Internal: " << internal_cycle_value << "\n";
    for(vector<FullGraph::Node>::const_iterator it = opt.begin(); it != opt.end(); ++it)
        cout << instance.graph.id(*it) << " ";
    #endif
    #endif

    //////////////////////////////////////////////////////////////////////////
    // Rotation Neighborhood
    //////////////////////////////////////////////////////////////////////////
    #ifdef FULLDEBUG
    cout << "\n\n|------- Rotating cycles ----------";
    cout.flush();
    #endif

    // Rotate each cycle looking for a better depot
    bool new_cycle = false;
    for(int i = 0; i < (int)opt.size(); ++i) {
        FullGraph::Node node = opt[i];
        FullGraph::Node intern_pred = opt[(i - 1 < 0)? (int)opt.size() - 1 : i - 1];
        FullGraph::Node intern_succ = opt[(i + 1 < (int)opt.size())? i + 1 : 0];

        deque< FullGraph::Node > &ll = (*cycles)[node];

        int best_deal = instance.dist[instance.graph.edge(intern_pred, node)] +
                        instance.dist[instance.graph.edge(node, intern_succ)];

        #ifdef FULLDEBUG
        cout << "\n\n> Rotating cycle " << i
             << "\n| Init. Depot: " << instance.graph.id(node)
             << "\t Pred: " << instance.graph.id(intern_pred)
             << "\t Succ: " << instance.graph.id(intern_succ)
             << "\t Value: " << best_deal;
        cout.flush();
        #endif

        deque< FullGraph::Node >::iterator itBest_node = ll.end();

        // Looking for a node in external cycle to be a better depot.
        for(deque< FullGraph::Node >::iterator itNode = ll.begin(); itNode != ll.end(); ++itNode) {
            int deal = instance.dist[instance.graph.edge(intern_pred, *itNode)] +
                       instance.dist[instance.graph.edge(*itNode, intern_succ)];

            #ifdef FULLDEBUG
            cout << "\n| Next node: " << instance.graph.id(*itNode)
                 << "\t Deal: " << deal
                 << "\t Best deal: " << best_deal;
            cout.flush();
            #endif

            if(best_deal > deal) {
                best_deal = deal;
                itBest_node = itNode;
            }
        }

        // If we found a better node, make it the new depot.
        if(itBest_node != ll.end()) {
            opt[i] = *itBest_node;
            (*cycle_sizes)[*itBest_node] = (*cycle_sizes)[node];
            (*cycle_sizes)[node] = 0;

            deque< FullGraph::Node > &ll2 = (*cycles)[opt[i]];
            ll2.swap(ll);

            ll2.push_back(node);
            rotate(ll2.begin(), itBest_node, ll2.end());
            ll2.pop_front();

            new_cycle = true;
//            solution_value += instance.alpha * (best_deal -
//                    (instance.dist[instance.graph.edge(intern_pred, node)] +
//                     instance.dist[instance.graph.edge(node, intern_succ)]));
            #ifdef FULLDEBUG
            cout << "\n|> New Depot: " << instance.graph.id(opt[i])
                 << "\n";

            for(deque< FullGraph::Node >::iterator itNode = ll.begin(); itNode != ll.end(); ++itNode)
                cout << instance.graph.id(*itNode) << " ";

            cout.flush();
            #endif
        }
    }

    // Do one more linkern to get better internal cycle.
    if(new_cycle) {
        initial.swap(opt);
        solution_value -= internal_cycle_value;
        internal_cycle_value = instance.alpha * concorde.getLinKernTour(initial, &opt);
        solution_value += internal_cycle_value;
    }

    //////////////////////////////////////////////////////////////////////////
    // Adjustment of the chromosome keys.
    //////////////////////////////////////////////////////////////////////////

    // Adjust the cycle size keys.
    for(unsigned i = 0; i < K; ++i)
        chromosome[i] = (*cycle_sizes)[opt[i]];

    // Adjust the remain keys.
    offset = K;
    for(unsigned i = 0; i < opt.size(); ++i) {
        FullGraph::Node node = opt[i];
        deque< FullGraph::Node > &ll = (*cycles)[node];

        // Depot/gateway key.
        chromosome[instance.graph.id(node) + K] = sorted_nodes[i].key;

        // External node keys
        for(deque< FullGraph::Node >::iterator it = ll.begin(); it != ll.end(); ++it)
            chromosome[instance.graph.id(*it) + K] = sorted_nodes[offset++].key;
    }

    // Clean up
    #ifdef _OPENMP
    #pragma omp critical
    {
    #endif
    delete cycles;
    delete cycle_sizes;
    #ifdef _OPENMP
    }
    #endif

    #ifdef DEBUG
    cout << "\n\n> Original Chromosome: " << original_solution_value << "\n";
    for(Population::Chromosome::iterator it = original_chr.begin(); it != original_chr.end(); ++it)
        cout << *it << " ";
        //cout << setiosflags(ios_base::right) << setw (4) << *it;

    cout << "\n\n> Improved Chromosome: " << solution_value << "\n";
    for(Population::Chromosome::iterator it = chromosome.begin(); it != chromosome.end(); ++it)
        cout << *it << " ";
        //cout << setiosflags(ios_base::right) << setw (4) << *it;
    cout << "\n------------------------------------------------------" << endl;
    #endif

    return solution_value;
}

//--------------------------------[ Decode ]----------------------------------//

void K_IMDMTSP_Decoder::getSolutionFromChromosome(const Population::Chromosome &chromosome,
                                                  IMDMTSP_Solution *solution) const {
    #ifdef DEBUG
    cout << "\n------------------------------------------------------\n"
         << "> Get Solution from Chromosome\n";
    copy(chromosome.begin(), chromosome.end(), ostream_iterator< Population::Allele >(cout, " "));
    cout << "\n" << endl;
    #endif

    const IMDMTSP_Instance &sol_inst = solution->instance;

    // Clean the solution
    solution->K = K;
    solution->cycles.clear();
    solution->cycles.resize(K);
    solution->value = 0.0;

    // Take the keys and sort then.
    PairVector sorted_nodes;
    sorted_nodes.reserve(chromosome.size() - K);

    for(unsigned i = K; i < chromosome.size(); ++i)
        sorted_nodes.push_back(KeyIndexPair(chromosome[i],i));

    sort(sorted_nodes.begin(), sorted_nodes.end(), greater<KeyIndexPair>());

    // Now, for each depot/gateway, build the cycle.
    FullGraph::Node current_node, last_node;
    unsigned offset = K;
    for(unsigned i = 0; i < K; ++i) {
        const unsigned cycle_size = chromosome[i] - 1;
        deque< FullGraph::Node > &ll = solution->cycles[i];

        // Begin with the depot.
        last_node = sol_inst.graph.nodeFromId(sorted_nodes[i].index - K);
        ll.push_back(last_node);

        #ifdef DEBUG
        cout << "> Cycle " << i << ": (" << sol_inst.graph.id(last_node) << ") ";
        double cycle_value = 0.0;
        #endif

        // Put all other nodes and sum the value of cycle edges.
        for(unsigned j = 0; j < cycle_size; ++j) {
            current_node = sol_inst.graph.nodeFromId(sorted_nodes[j + offset].index - K);
            ll.push_back(current_node);

            int edge_value = sol_inst.dist[sol_inst.graph.edge(last_node, current_node)];
            solution->value += ((cycle_size != 1)? edge_value : (2 * edge_value));

            last_node = current_node;

            #ifdef DEBUG
            cout << sol_inst.graph.id(last_node) << " ";
            cycle_value += ((cycle_size != 1)? edge_value : (2 * edge_value));
            #endif
        }

        if(cycle_size > 1)
            solution->value += sol_inst.dist[sol_inst.graph.edge(last_node, *(ll.begin()))];

        offset += cycle_size;

        #ifdef DEBUG
        if(cycle_size > 1)
            cycle_value += sol_inst.dist[sol_inst.graph.edge(last_node, *(ll.begin()))];

        cout << " | value: " << cycle_value << endl;
        #endif
    }

    // Sum the internal cycle.
    #ifdef DEBUG
    double cycle_value = 0.0;
    #endif
    vector< deque< FullGraph::Node > >::iterator it1, it2;
    for(it2 = solution->cycles.begin(); it2 != solution->cycles.end(); ) {
        it1 = it2++;
        if(it2 != solution->cycles.end()) {
            solution->value += instance.alpha * sol_inst.dist[sol_inst.graph.edge(*(it1->begin()),*(it2->begin()))];

            #ifdef DEBUG
            cycle_value += instance.alpha * sol_inst.dist[sol_inst.graph.edge(*(it1->begin()),*(it2->begin()))];
            #endif
        }
    }
    solution->value += instance.alpha * sol_inst.dist[sol_inst.graph.edge(*(it1->begin()),*(solution->cycles.begin()->begin()))];

    #ifdef DEBUG
    cycle_value += instance.alpha * sol_inst.dist[sol_inst.graph.edge(*(it1->begin()),*(solution->cycles.begin()->begin()))];
    cout << "> Internal cycle: " << cycle_value
         << "\n\n> Value: " << solution->value
         << "\n------------------------------------------------------" << endl;
    #endif
}
