#========================================
# USER RESERVED -- The following are reserved for users to set on the
# command line.  Makefiles should not set these.  These variables are
# for C/C++ compilation, and linking.
## -pg is for use with gprof.  Use on CFLAGS and on LDFLAGS
#CFLAGS		= -Wall -pg
#CFLAGS 	= -Winline  
#LDFLAGS	= -pg

# OPTIMIZE with the -O option.  Override from the command line for
# building debug versions.
#
OPTFLAGS	= -O3 -funroll-loops -DNDEBUG=1

#========================================
## To compile against the HMMoC-BFloat-Algebra library, use:
ALGEBRA_CFLAGS 	= -I.
ALGEBRA_LDFLAGS = -L.
ALGEBRA_LIBS	= -l$(ALGEBRA_LIB_NAME)


#========================================
### For Boost support, eg "g++ -I/sw/include/boost-1_44_0/ main.cpp Algebra.o -L/sw/lib/ -lboost_serialization": (but here we're assuming you've created links "boost-include" and "boost-lib" in the current dir)..
#BOOST_CFLAGS	= -DNBOOST_SERIALIZATION
BOOST_CFLAGS 	= -I./boost-include
BOOST_LDFLAGS 	= -L./boost-lib
BOOST_LIBS	= -lboost_serialization

###==============================================

INCS = Algebra.hpp
SOURCES = Algebra.cpp
OBJS = $(SOURCES:.cpp=.o)

ALGEBRA_LIB_NAME = HMMoC-BFloat-Algebra
ALGEBRA_LIB = $(LIB_PFX)$(ALGEBRA_LIB_NAME)$(LIB_SFX)


TESTS_INCS = $(INCS)
TESTS_SOURCES = Tests.cpp
TESTS_OBJS = $(TESTS_SOURCES:.cpp=.o)

TESTS = tests


default: $(ALGEBRA_LIB)

$(ALGEBRA_LIB): $(OBJS)
	     $(AR) $(ALGEBRA_LIB) $(OBJS)

$(TESTS): $(TESTS_SOURCES) $(TESTS_INCS) $(TESTS_OBJS)
	   $(CXX_LINK) $(ALGEBRA_LDFLAGS) $(ALGEBRA_LIBS) -o $(TESTS) $(TESTS_OBJS)

all: $(ALGEBRA_LIB) $(TESTS)

check: $(TESTS)
test: $(TESTS)

## Recompile if the includes are modified ...
$(OBJS): $(SOURCES) $(INCS)

.PHONY: clean
clean:
	rm -f $(OBJS) $(TESTS_OBJS) $(ALGEBRA_LIB) $(TESTS)

#========================================
# FILE EXTENSIONS.  Extensions and prefixes for different types of
# files change from platform to platform.  Hide these in macros so
# that we can more easily cut and paste between makefiles.
o		= .o
EXE_SFX		= 
SCRIPT_SFX 	= 
LIB_PFX		= lib
LIB_SFX		= .a
LIB_SHARED_SFX	= .so
TMPLIB		= libtemp.a

# FILE TOOLS
AR 	= ar rcs
CHMOD 	= chmod
CP	= cp
GREP	= grep
MKDIR 	= mkdir
MUNCH 	= stepmunch
MV	= mv
NM 	= nm
RANLIB	= ranlib
RM 	= rm -f
RMDIR 	= rm -rf
STRIP	= strip
UNZIP 	= unzip
ZIP 	= zip


#========================================
# ANSI C Compile and Link
#
CC		= gcc
CC_COMPILE	= $(CC) -c $(OPTFLAGS) $(CFLAGS) $(CC_CFLAGS) $(CC_SYSCFLAGS)
CC_LINK		= $(CC) $(LDFLAGS) $(CC_LDFLAGS) $(CC_SYSLDFLAGS) $(CC_LIBS)
CC_CFLAGS 	= $(BOOST_CFLAGS)
CC_LDFLAGS	= $(BOOST_LDFLAGS)
CC_LIBS		= $(BOOST_LIBS)

# Global system things used for compilation, static linking, etc.
CC_SYSCFLAGS 	= -I.
CC_SYSLDFLAGS 	=
CC_SYSLIBS	=

#========================================
# C++ Compile and Link
#
CXX		= g++
CXX_COMPILE	= $(CXX) -c  $(OPTFLAGS) $(CFLAGS) $(CXX_CFLAGS) $(CXX_SYSCFLAGS)
CXX_LINK	= $(CXX) $(LDFLAGS) $(CXX_LDFLAGS) $(CXX_SYSLDFLAGS) $(CXX_LIBS)
CXX_CFLAGS 	= $(BOOST_CFLAGS)
CXX_LDFLAGS	= $(BOOST_LDFLAGS)
CXX_LIBS	= $(BOOST_LIBS)

# The force flags are used for C/C++ compilers that select the
# language based on the file naming conventions.  Some C++ source
# may be in files with C naming conventions.
CXX_FORCE	= 

# System Flags -- Things for static linking or making sure that the
# compiler understands that a file is a C++ file or whatever.  These
# usually change from platform to platform.
CXX_SYSCFLAGS 	= -I.
CXX_SYSLDFLAGS 	= 
CXX_SYSLIBS	= 

# Compilation Rules -- Repeat the rules for all of the different
# naming conventions.
#
.cxx.o:	; $(CXX_COMPILE) $<
.cpp.o:	; $(CXX_COMPILE) $<
.cc.o:	; $(CXX_COMPILE) $<
.C.o:	; $(CXX_COMPILE) $<

.cxx:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)
.cpp:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)
.cc:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)
.C:	
	$(CXX_COMPILE) $<
	$(CXX_LINK) -o $@ $*.o $(LIBRARIES)

# for legacy reasons also compile .c as c++
.c.o:	; $(CXX_COMPILE) $(CXX_FORCE) $<
.c:	
	$(CXX_COMPILE) $(CXX_FORCE) $<

