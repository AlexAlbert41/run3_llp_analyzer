#common definitions go in here
SHELL = /bin/sh
CXX           = g++
CXXFLAGS      = -O -Wall -fPIC -g -ansi -Wextra
CXX11FLAGS    = -std=c++0x
CXX14FLAGS    = -std=c++17
LD            = g++
LDFLAGS       = -O
LIBDIR       = $(WORKDIR)

OS := $(shell uname -s)
ifeq ($(OS),Darwin)
        SOFLAGS       = -dynamiclib
else
        SOFLAGS       = -shared
endif

ifndef ROOTSYS
$(error *** Please set up Root)
endif

#find out about the installed ROOT
ROOTCONFIG   := $(ROOTSYS)/bin/root-config
ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags) $(shell $(ROOTCONFIG) --libs --nonew --glibs)

#append to options
CXXFLAGS     += $(ROOTCFLAGS)
LDFLAGS      += $(ROOTLDFLAGS)

# fast jet appde
CXXFLAGS += $(shell ../fastjet-install/bin/fastjet-config --cxxflags) 



LIBS = -lRooFitCore -lRooFit -lMinuit -lHtml -lMathMore -lPyROOT -lFoam -lRooStats -lTreePlayer -lTMVA
LIBS += $(shell ../fastjet-install/bin/fastjet-config --libs) -I${PWD}/fastjet/lib
ifneq ($(CMSSW_BASE),)
        ROOTFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')
        ROOTFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

        CINTINCLUDES := -I$(ROOTFIT_INCLUDE)
        CXXFLAGS     += -I$(ROOTFIT_INCLUDE)
        LDFLAGS      += -L$(ROOTFIT_LIBDIR)
endif
