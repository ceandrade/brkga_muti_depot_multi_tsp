/******************************************************************************
 * imdmtsp_decoder.hpp: Interface for K-IMDMTSP Decoder class.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2012-2018, Carlos Eduardo de Andrade.
 *     All Rights Reserved.
 *
 * Created on : Jun 06, 2012 by andrade
 * Last update: Jul 23, 2012 by andrade
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

#ifndef IMDMTSP_DECODER_HPP
#define IMDMTSP_DECODER_HPP

#include <deque>

#include "brkga_decoder.hpp"
#include "mtrand.hpp"
#include "k_imdmtsp_instance.hpp"
#include "imdmtsp_solution.hpp"
#include "concorde_wrapper.hpp"

/**
 * \brief K-Interconnected Multi-Depot Multi Traveling Salesman Problem Decoder
 *        description.
 * \author Carlos Eduardo de Andrade <andrade@ic.unicamp.br>
 * \date 2011
 *
 * This class contains the decoder method that takes a chromosome from
 * BRKGA algorithm and transforms to a solution of K-IMDMTSP.
 */
class K_IMDMTSP_Decoder: public BRKGA_Decoder {
    public:
        /** \name Constructor and Destructor */
        //@{
        /** \brief Default Constructor.
         * \param instance K-IMDMTSP instance.
         * \param rng the random number generator.
         * \param concorde the Concorde Wrapper used to local search routines.
         */
        K_IMDMTSP_Decoder(const K_IMDMTSP_Instance &instance, MTRand &rng,
                          ConcordeWrapper &concorde);

        /** \brief Destructor. */
        virtual ~K_IMDMTSP_Decoder();
        //@}

        /** \name Mandatory Method for BRKGA */
        //@{
        /** \brief Extracts and improve a solution given a chromosome.
         *
         * This method decodes a IMDMTSP solution from a chromosome, doing
         * local search procedures to improve the solution. Such improvements
         * are reflected on chromosome in the ending of the process.
         *
         * \param chromosome A vector of doubles represent a problem solution.
         * \return a double with fitness value.
         */
        virtual double decode(Population::Chromosome &chromosome);
        //@}

        /** \name Member Methods */
        //@{
        /** \brief This method extracts a solution from given chromosome.
         *  \param[in] chromosome to be decoded.
         *  \param[out] solution a reference to a new solution.
         */
        void getSolutionFromChromosome(const Population::Chromosome &chromosome,
                                       IMDMTSP_Solution *solution) const;
        //@}

    protected:
        /** \name Local helper classes and structures. */
        //@{
        /// Used to sort the key/alleles
        class KeyIndexPair {
            public:
                unsigned key;
                unsigned index;

                KeyIndexPair(unsigned _key = 0, unsigned _index = 0):
                    key(_key), index(_index) {}

                inline bool operator>(const KeyIndexPair &x) const {
                    return this->key > x.key;
                }

            #ifdef FULLDEBUG
                friend ostream& operator<<(ostream& output, const KeyIndexPair& kip) {
                    output << "(" <<  kip.key << ", " << kip.index <<")";
                    return output;
                }
            #endif
        };

        typedef vector< KeyIndexPair > PairVector;
        //@}

        /** \name Local attributes */
        //@{
        /// The K-IMDMTSP instance used to decode solutions.
        const K_IMDMTSP_Instance &instance;

        /// The Random Number Generator.
        MTRand &rng;

        /// The Concorde Wrapper.
        ConcordeWrapper &concorde;

        /// Some local instance data
        const unsigned K;
        const unsigned max_cycle_size;

        /// Keep the sorted nodes in memory.
        vector< PairVector > sorted_nodes_per_thread;
        //@}
};

#endif //IMDMTSP_DECODER_HPP
