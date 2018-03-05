################################################################################
# Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
#
# (c) Copyright 2012-2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on : Jul 09, 2011 by andrade
# Last update: Mar 06, 2016 by andrade
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
# User Defines
###############################################################################

# Optimization switch
#OPT = opt
#USER_DEFINES += -DNDEBUG

# Set debug mode
USER_DEFINES += -DDEBUG
USER_DEFINES += -DFULLDEBUG 

# Set tuning mode
#USER_DEFINES += -DTUNING

###############################################################################

ifeq ($(shell uname -s),Darwin)
    SYSTEM = macosx
else
    SYSTEM = linux
endif

include Makefile.inc

###############################################################################
# Build options
###############################################################################

EXE = brkga-kimdmtsp
#EXE = teste

###############################################################################
# The user source files
###############################################################################

# Include dir
USER_INCDIRS = \
	./include \
	./brkga \
	./lib \
	./concorde/INCLUDE

# Source directories (to clean up object files)
SRC_DIRS = . \
	./problem_defs/ \
	./decoders/ \
	./concorde_wrapper/

# Object files
OBJS = \
	./brkga/population.o \
	./problem_defs/imdmtsp_instance.o \
	./problem_defs/k_imdmtsp_instance.o \
	./problem_defs/imdmtsp_solution.o \
	./decoders/k_imdmtsp_decoder.o \
	./decoders/dummy_decoder.o \
	./concorde_wrapper/concorde_wrapper.o

# Main function
BRKGA_OBJ += ./main_brkga.o

###############################################################################
# Lib and include definitions
###############################################################################

##############################
# Concorde libs
##############################

CONCORDE_PATH = ./concorde
CONCORDE_LIBS = \
	$(CONCORDE_PATH)/LINKERN/linkern.a \
	$(CONCORDE_PATH)/EDGEGEN/edgegen.a \
	$(CONCORDE_PATH)/KDTREE/kdtree.a \
	$(CONCORDE_PATH)/FMATCH/fmatch.a \
	$(CONCORDE_PATH)/UTIL/util.a

##############################
# Mauro's Lib
##############################

MAURO_LIB = ./mauro_lib/mauro_lib.a

###############################################################################
# Lib and include definitions
###############################################################################

ifeq ($(SYSTEM), macosx)
	INCLUDE_PATH = /opt/local
	USER_LIBDIRS = /opt/local/lib 	
#	USER_LIBDIRS = /opt/local/my_boost/lib 
#	USER_INCDIRS += /usr/local/my_boost/include \
#					/usr/local/include
else
	INCLUDE_PATH = ~/my_usr
	USER_LIBDIRS = ~/my_usr/lib
endif

USER_INCDIRS += \
	$(INCLUDE_PATH)/include \
	$(INCLUDE_PATH)/include/boost

LIBS = -lm -lemon

#Boost libraries
USER_DEFINES+= -DBOOST_FILESYSTEM_NO_DEPRECATED

ifeq ($(SYSTEM), macosx)
	LIBS += -lboost_system-mt -lboost_chrono-mt -lboost_timer-mt -lboost_filesystem-mt
else
	LIBS += -lboost_system -lboost_chrono -lboost_timer -lboost_filesystem
endif

##############################
# Consolidate paths
##############################

# Consolidate include paths
USER_INCDIRS += #\
	#$(CPLEXINCDIRS)
	
# Consolidate include paths
USER_LIBDIRS += #\
	#$(CPLEXLIBDIRS)

# Libraries necessary to link.
LIBS += -lemon
	
# Compiler flags
USER_FLAGS += 

###############################################################################
# Compiler defs
###############################################################################

# Necessary Include dirs
# Put -I in front of dirs
INCLUDES = `for i in $(USER_INCDIRS); do echo $$i | sed -e 's/^/-I/'; done`

# Necessary library dirs
# Put -L in front of dirs
LIBDIRS = `for i in $(USER_LIBDIRS); do echo $$i | sed -e 's/^/-L/'; done`

###############################################################################
# Build Rules
###############################################################################

all: $(EXE)

.PHONY: all concorde mauro_lib strip clean doc docclean depclean
.SUFFIXES: .cpp .o

$(EXE): concorde mauro_lib $(OBJS) $(BRKGA_OBJ)
	@echo "--> Linking objects... "
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJS) $(BRKGA_OBJ) $(CONCORDE_LIBS) $(MAURO_LIB) $(LIBDIRS) $(LIBS) -o $@
	@echo

.cpp.o:
	@echo "--> Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(USER_DEFINES) -c $< -o $@
	@echo
	
concorde:
	@echo "|---------------------[ Compiling Concorde ]---------------------|"
	make -C concorde
	@echo "|-----------------[ End of compiling Concorde ]------------------|"

mauro_lib:
	@echo "|--------------------[ Compiling Mauro Lib ]---------------------|"
	make -C mauro_lib
	@echo "|----------------[ End of compiling Mauro Lib ]------------------|"

strip: $(EXE)
	@echo "-->Stripping $<..."
	strip $<
	@echo

doc:
	doxygen Doxyfile

clean:
	@echo "--> Cleaning compiled..."
	rm -rf `for i in $(SRC_DIRS); do echo $$i*~; done` 
	rm -rf Debug
	make -C concorde clean
	rm -f concorde/EDGEGEN/edgegen
	rm -f concorde/FMATCH/fmatch
	rm -f concorde/LINKERN/linkern
	make -C mauro_lib clean
	rm *o
	rm -rf doc
	rm -rf $(EXE) $(OBJS) $(BRKGA_OBJ) *.dSYM
	
