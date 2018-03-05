################################################################################
# Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
#
# (c) Copyright 2012-2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
#  This code is released under LICENSE.md.
#
#  Created on : Jul 02, 2012 by andrade
#  Last update: Jul 02, 2012 by andrade
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
###############################################################################

###############################################################################
# Lib and include definitions
###############################################################################

include ../Makefile.inc

LIBNAME = mauro_lib.a

OBJS = \
	./instance.o \
	./solution.o \
	./kmeans_heuristic.o
	
CONCORDE_PATH = ../concorde

###############################################################################
# Build Rules
###############################################################################

all: $(LIBNAME)

.PHONY: all clean
.SUFFIXES: .cpp .o

$(LIBNAME): $(OBJS)
	@echo "--> Building $<..."
	$(AR) $(ARFLAGS) $(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBNAME)

.cpp.o:
	@echo "--> Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -I$(CONCORDE_PATH)/INCLUDE -c $< -o $@
	@echo

clean:
	@echo "--> Cleaning compiled..."
	rm *.o *.a
	
