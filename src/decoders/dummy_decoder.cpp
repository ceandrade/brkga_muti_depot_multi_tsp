/******************************************************************************
 * dummy_decoder.cpp: Implements a dummy decoder for BRKGA.
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

#include <iostream>
#include <algorithm>
#include <iterator>
using namespace std;

#include "dummy_decoder.hpp"

//------------------[ Default Constructor and Destructor ]--------------------//

DummyDecoder::DummyDecoder() {}

DummyDecoder::~DummyDecoder() {}

//--------------------------------[ Decode ]----------------------------------//

double DummyDecoder::decode(Population::Chromosome& chromosome) const {
    cout << "\n";
    copy(chromosome.begin(), chromosome.end(), ostream_iterator<Population::Allele>(cout, " "));
    cout << endl;
    return 0.0;
}

