/******************************************************************************
 * imdmtsp_solution.hpp: Interface for IMDMTSP Solution.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2012-2018, Carlos Eduardo de Andrade.
 *     All Rights Reserved.
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

#ifndef IMDMTSP_SOLUTION_HPP
#define IMDMTSP_SOLUTION_HPP

#include <iostream>
#include <vector>
#include <deque>
using std::vector;
using std::deque;

#include "imdmtsp_instance.hpp"

/**
 * \brief Interconnected Multi-Depot Multi Traveling Salesman Solution.
 *        description.
 * \author Carlos Eduardo de Andrade <andrade@ic.unicamp.br>
 * \date 2012
 *
 * This class represents a solution to IMDMTSP given an instance.
 */
class IMDMTSP_Solution {
    public:
        /** \name Constructor and Destructor */
        //@{
        /** \brief Default Constructor.
         * \param instance K-IMDMTSP instance.
         */
        IMDMTSP_Solution(const IMDMTSP_Instance &instance);

        /** \brief Destructor. */
        virtual ~IMDMTSP_Solution();
        //@}

        /** \name Informational methods */
        //@{
        /** \brief Make an image of the instance if it has coordinates. */
        void draw(const string &file_name);

        /** \brief Print in text format. */
        friend std::ostream& operator<<(std::ostream &os, const IMDMTSP_Solution &solution);
        //@}

    public:
        /** \name Data Members. */
        //@{
        /** The number of external cycles.
         * \warning MUST BE 3 <= K <= (num_nodes / 3)
         */
        unsigned K;

        /// The cycles. The first node of each cycle is the depot/gateway.
        vector< deque< FullGraph::Node > > cycles;

        /// The solution value.
        double value;

        /// The K-IMDMTSP instance which this solution is.
        const IMDMTSP_Instance &instance;
        //@}
};

#endif // IMDMTSP_SOLUTION_HPP
