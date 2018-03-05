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

Instance::Instance(FullGraph &_graph, FullGraph::NodeMap<PointType> &_coords,
             FullGraph::EdgeMap<DistanceType> &_dist):

     graph(_graph), coords(_coords), dist(_dist) {}

void Instance::load(string filename)
{
	ConcordeWrapper concorde;
	concorde.loadTSP(const_cast<char *>(filename.c_str()), &graph, &dist,
			&coords, &has_coords);

	this->n = graph.nodeNum();
}

