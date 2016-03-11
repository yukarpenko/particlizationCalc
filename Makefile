ODIR           = obj

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXX           = g++
F             = gfortran
CXXFLAGS      = -fPIC -O3 -std=c++11
LD            = g++
LDFLAGS       = -O3
FFLAGS        = -fPIC $(ROOTCFLAGS) -O3

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

_HYDROO        = DecayChannel.o ParticlePDG2.o DatabasePDG2.o UKUtility.o gen.o \
                particle.o main.o
 
# VPATH = src:../UKW
HYDROO = $(patsubst %,$(ODIR)/%,$(_HYDROO))

TARGET = calc
#------------------------------------------------------------------------------

$(TARGET): $(HYDROO)
		$(LD) $(LDFLAGS) $^ -o $@ $(LIBS) -lgfortran
		@echo "$@ done"

clean:
		@rm -f $(ODIR)/*.o $(TARGET)

$(ODIR)/%.o: src/%.cpp src/const.h
		$(CXX) $(CXXFLAGS) -c $< -o $@

