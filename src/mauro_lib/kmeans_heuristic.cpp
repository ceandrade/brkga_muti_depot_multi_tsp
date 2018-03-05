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
#include "../include/concorde_wrapper.hpp"

#include <list>

Solution* kmeans_heuristic(Instance *instance)
{
	int K = instance->K;
	FullGraph &graph = instance->graph;
	//TimeStamp time_stamp;
	//time_stamp.stamp();

	vector<PointType> centroids(K);
	FullGraph::NodeMap<int> centroid_of(graph);
	vector<int> size(K);

	set<int> initial_centroids;

	while (int(initial_centroids.size()) < K)
		initial_centroids.insert(rnd[instance->n]);

	//initialization
	{
		int k = 0;
        for (auto it = initial_centroids.begin();
		//for (auto it = initial_centroids.begin();
				it != initial_centroids.end(); ++it, k++)
		{
			centroids[k] = instance->coords[graph(*it)];
		}
	}

	for (int i = 0; i < 20; i++)
	{
		//assignment step
		for (FullGraph::NodeIt v(graph); v != INVALID; ++v)
		{
			DistanceType min_norm_square = numeric_limits<DistanceType>::max();
			int closest_centroid = 0;
			PointType point = instance->coords[v];

			for (int k = 0; k < K; k++)
			{
				DistanceType d = (point - centroids[k]).normSquare();
				if (d < min_norm_square)
				{
					min_norm_square = d;
					closest_centroid = k;
				}
			}
			centroid_of[v] = closest_centroid;
		}

		//update step
		for (int k = 0; k < K; k++)
		{
			centroids[k] = dim2::makePoint(0, 0);
			size[k] = 0;
		}

		for (FullGraph::NodeIt v(graph); v != INVALID; ++v)
		{
			int c = centroid_of[v];
			size[c]++;
			centroids[c] += instance->coords[v];
		}
		for (int k = 0; k < K; k++)
			if (size[k] != 0)
				centroids[k] /= size[k];
	}

	//////////////////////////////////////////////////////////////////////////
	// Fix empty centroids (Carlos Eduardo).
	// Random choose non-empty centroid and get the closest node.
	//////////////////////////////////////////////////////////////////////////

	for (int k = 0; k < K; ++k) {
	    if (size[k] == 0) {
	        // Random choose a non-empty centroid.
	        int chosen_centroid;

	        do chosen_centroid = rnd.integer(K);
	        while(size[chosen_centroid] < 2);

	        // Take the nodes from this centroid.
	        list<FullGraph::Node> nodes;
            for (FullGraph::NodeIt v(graph); v != INVALID; ++v)
                if(centroid_of[v] == chosen_centroid)
                    nodes.push_back(v);

	        // Take the farther node.
	        DistanceType min_norm_square = numeric_limits<DistanceType>::max();
	        FullGraph::Node closest_node = *nodes.begin();

	        for(list<FullGraph::Node>::iterator v = nodes.begin(); v != nodes.end(); ++v) {
	            PointType point = instance->coords[*v];
                DistanceType d = (point - centroids[chosen_centroid]).normSquare();
                if (d < min_norm_square) {
                    min_norm_square = d;
                    closest_node = *v;
                }
	        }

            centroid_of[closest_node] = k;
            size[k] = 1;
            --size[chosen_centroid];
	    }
	}

	//////////////////////////////////////////////////////////////////////////

	vector<vector<FullGraph::Node> > subcycles(K + 1);

	for (int k = 0; k <= K; k++)
		subcycles[k].clear();

	for (FullGraph::NodeIt v(graph); v != INVALID; ++v)
	{
		int c = centroid_of[v];
		subcycles[c].push_back(v);
	}

//	for (int k = 0; k < K; k++)
//		for (FullGraph::NodeIt v(graph); v != INVALID; ++v)
//			if (centroid_of[v] == k)
//			{
//			    cout << "\n" << k << "> " << instance->graph.id(v);
//
//				//add node to the internal cycle
//				subcycles[K].push_back(v);
//				break;
//			}
	for(int k = 0; k < K; k++) {
	    subcycles[K].push_back(*subcycles[k].begin());
	}

//	for (int k = 0; k <= K; k++)
//	{
//		int cycle_size = int(subcycles[k].size());
//		if (cycle_size < 3)
//			return NULL;
//	}
	return new Solution(instance, subcycles);
}

