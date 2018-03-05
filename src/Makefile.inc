################################################################################
# Author: Carlos Eduardo de Andrade <ce.andrade@gmail.com>
#
# (c) Copyright 2012-2018, Carlos Eduardo de Andrade. All Rights Reserved.
#
# This code is released under LICENSE.md.
#
# Created on : Jun 07, 2012 by andrade
# Last update: Apr 26, 2014 by andrade
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
# User Flags
###############################################################################

SYSTEM = macosx

# Optimization switch
#OPT = opt
#USER_DEFINES += -DNDEBUG

# Set debug mode
USER_DEFINES += -DDEBUG
USER_DEFINES += -DFULLDEBUG 

# Set tuning mode
#USER_DEFINES += -DTUNING

# Define the usage of C++11. This code compiles only in modern compilers,
# for instance, gcc >= 4.8
USER_FLAGS += -std=c++11 

# Machine type
MACHINE = native

# Compiler flags for debugging
ifneq ($(OPT), opt)
	USER_FLAGS += -ggdb3 -fexceptions
else
	USER_FLAGS += -O3 -march=$(MACHINE) -mtune=$(MACHINE) \
		-fomit-frame-pointer -ftracer \
		-funroll-loops -fpeel-loops -fprefetch-loop-arrays
endif

# Flags
#USER_FLAGS += -Wall -Wextra -Wcast-align -Wcast-qual -Wdisabled-optimization \
#				-Wformat=2 -Winit-self -Wmissing-format-attribute -Wshadow \
#				-Wpointer-arith -Wredundant-decls -Wstrict-aliasing=2 \
#				-Wfloat-equal -Weffc++

# Lemon has serious problems with this flag and GCC ignores the 
# "ignore" pragma.
#USER_FLAGS += -Winline
				
# Used to analyze format of structures 
#USER_FLAGS += -Wpacked -Wpadded

# Paranoid!!!
#USER_FLAGS += -Wconversion
#USER_FLAGS += -Wunsafe-loop-optimizations # doesn't work very well with C++11 foreach

USER_FLAGS += -pthread -fopenmp 

###############################################################################
# Lib and include definitions
###############################################################################

ifeq ($(SYSTEM), macosx)
	INCLUDE_PATH = /opt/local
	USER_LIBDIRS = /opt/local/my_boost/lib 
	USER_INCDIRS += /usr/local/my_boost/include \
					/usr/local/include
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
LIBS += -lboost_system -lboost_chrono -lboost_timer -lboost_filesystem

###############################################################################
# Compiler and linker defs
###############################################################################

# C Compiler command and flags
CC = gcc
CXX = g++
#ifneq ($(SYSTEM), macosx)
#    CC = gcc-4.8
#	CXX = g++-4.8
#endif

CFLAGS = $(USER_FLAGS)
CXXFLAGS = $(USER_FLAGS)

# Linker
LD = ld

# Lib maker commands
AR = ar
ARFLAGS	= rv
RANLIB = ranlib

# Other includes
RM = rm
SHELL = /bin/bash

###############################################################################
# Compiler defs
###############################################################################

# Necessary Include dirs
# Put -I in front of dirs
INCLUDES = `for i in $(USER_INCDIRS); do echo $$i | sed -e 's/^/-I/'; done`

# Necessary library dirs
# Put -L in front of dirs
LIBDIRS = `for i in $(USER_LIBDIRS); do echo $$i | sed -e 's/^/-L/'; done`
