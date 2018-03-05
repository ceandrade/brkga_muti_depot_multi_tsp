/******************************************************************************
 * dummy_decoder.hpp: Interface for a dummy decoder for BRKGA.
 *
 * Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
 *
 * (c) Copyright 2012-2018, Carlos Eduardo de Andrade. All Rights Reserved.
 *
 * Created on : Jun 20, 2012 by andrade
 * Last update: Jun 20, 2012 by andrade
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

#ifndef DUMMY_DECODER_HPP
#define DUMMY_DECODER_HPP

#include "brkga_decoder.hpp"

/**
 * \brief Dummy Decoder.
 *
 * \author Carlos Eduardo de Andrade <andrade@ic.unicamp.br>
 * \date 2012
 *
 * This class is a dummy decoder for BRKGA. The decode method only
 * print the chromosome.
 */
class DummyDecoder: public BRKGA_Decoder {
    public:
        /** \name Constructor and Destructor */
        //@{
        /** \brief Default Constructor. */
        DummyDecoder();

        /** \brief Destructor. */
        virtual ~DummyDecoder();
        //@}

        /** \name Mandatory Method for BRKGA */
        //@{
        /** \brief Only print the chromosome.
         *
         * \param chromosome A vector of doubles represent a problem solution.
         * \return 0.0.
         */
        virtual double decode(Population::Chromosome& chromosome) const;
        //@}
};

#endif //DUMMY_DECODER_HPP

