/******************************************************************************
 * k_imdmtsp_instance.hpp: Interface for K-IMDMTSP Instance class.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2012-2018, Carlos Eduardo de Andrade.
 *     All Rights Reserved.
 *
 * Created on : Jun 18, 2012 by andrade
 * Last update: Jul 17, 2012 by andrade
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

#ifndef K_IMDMTSP_INSTANCE_HPP
#define K_IMDMTSP_INSTANCE_HPP

#include "imdmtsp_instance.hpp"

/**
 * \brief K-Interconnected Multi-Depot Multi Traveling Salesman Problem Instance
 *        description.
 * \author Carlos Eduardo de Andrade <andrade@ic.unicamp.br>
 * \date 2012
 *
 * This class represents a K-IMDMTSP instance.
 */
class K_IMDMTSP_Instance: public IMDMTSP_Instance {
    public:
        /** \name Constructor and Destructor */
        //@{
        /** \brief Default Constructor. */
        K_IMDMTSP_Instance();

        /** \brief Destructor. */
        virtual ~K_IMDMTSP_Instance();

    public:
        /** \name Data Members. */
        //@{

        /** The number of vertices in the interior cycle.
         * \warning MUST BE K <= (num_nodes / 3)
         */
        unsigned K;

        /** Maximum external cycle size.
         * \warning MUST BE \f$\lceil num\_nodes/K \rceil\f$.
         */
        unsigned max_cycle_size;
        //@}
};

#endif // K_IMDMTSP_INSTANCE_HPP

