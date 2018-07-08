ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
EXTRALIBS += $(shell $(ROOTSYS)/bin/root-config --libs)
LDLIBS += $(shell $(ROOTSYS)/bin/root-config --libs)

#ROOTCFLAGS    = $(shell /usr/bin/root-config --cflags)
#ROOTLIBS      = $(shell /usr/bin/root-config --libs)
#ROOTGLIBS     = $(shell /usr/bin/root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -Wall -fPIC -Wno-deprecated

NGLIBS         = -Wl,--no-as-needed $(ROOTGLIBS) 
NGLIBS        += -lMinuit

CXXFLAGS      += $(ROOTCFLAGS)
CXX           += -I./
LIBS           = $(ROOTLIBS) 

GLIBS          = $(filter-out -lNew, $(NGLIBS))

CXX	      += -I./obj/
OUTLIB	      = ./obj/
.SUFFIXES: .cc,.C
.PREFIXES: ./obj/

#----------------------------------------------------#

all: ascii2root csv2root analysis1 analysis2 ascii2root_v2 plot plot_quad plot2 convert csv2root_blm ascii2root_v3 convert_v3 ascii2root_v4 ascii2root_v5 convert_v5 analysis_v5 ascii2root_v6 analysis_v6 analysis_v7 csv2root_bpm analysis_v8 convert_v5_2 analysis_v5_2 GetUT_GMT analysis_v9 ascii2root_trk convert_trk analysis_trk trackreco_trk ascii2root_common convert_common analysis_common

clean:
	rm -f *~
	rm -f ascii2root
	rm -f convert
	rm -f convert_v3
	rm -f convert_v5
	rm -f ascii2root_v2
	rm -f ascii2root_v3
	rm -f ascii2root_v4
	rm -f ascii2root_v5
	rm -f ascii2root_v6
	rm -f analysis1
	rm -f analysis2
	rm -f analysis3
	rm -f analysis4
	rm -f analysis_v5
	rm -f analysis_v6
	rm -f analysis_v7
	rm -f analysis_v8
	rm -f analysis_v9
	rm -f csv2root
	rm -f plot
	rm -f plot2
	rm -f plot_quad
	rm -f csv2root_blm
	rm -f csv2root_bpm
	rm -f convert_v5_2
	rm -f analysis_v5_2
	rm -f GetUT_GMT
	rm -f ascii2root_trk
	rm -f convert_trk
	rm -f analysis_trk
	rm -f trackreco_trk
