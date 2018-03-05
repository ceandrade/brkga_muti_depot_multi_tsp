/******************************************************************************
 * imdmtsp_instance.hpp: Interface for IMDMTSP Instance class.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2012-2018, Carlos Eduardo de Andrade.
 *     All Rights Reserved.
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

#ifndef IMDMTSP_INSTANCE_HPP
#define IMDMTSP_INSTANCE_HPP

#include <lemon/full_graph.h>
#include <lemon/dim2.h>
using namespace lemon;

/**
 * \brief Interconnected Multi-Depot Multi Traveling Salesman Problem Instance
 *        description.
 * \author Carlos Eduardo de Andrade <andrade@ic.unicamp.br>
 * \date 2012
 *
 * This class represents a generic IMDMTSP instance as base for variations of
 * the problem.
 */
class IMDMTSP_Instance {
    public:
        /** \name Constructor and Destructor */
        //@{
        /** \brief Default Constructor. */
        IMDMTSP_Instance();

        /** \brief Destructor. */
        virtual ~IMDMTSP_Instance();
        //@}

        /** \brief Make an image of the instance if it has coordinates. */
        void draw();

    public:
        /** \name Data Members. */
        //@{
        /// Instance name
        std::string name;

        /// The full graph.
        FullGraph graph;

        /// The edges weight.
        FullGraph::EdgeMap<int> dist;

        /// Indicates if we have the coordinates of the nodes.
        bool has_coords;

        /// The 2D coordinates of each node.
        FullGraph::NodeMap< dim2::Point<double> > coords;

        /// Cost factor of internal/depot cycle.
        double alpha;
        //@}
};

#endif // IMDMTSP_INSTANCE_HPP

