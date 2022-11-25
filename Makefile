PROGRAM = stv
CPLEX=/Applications/CPLEX_Studio221
CPLEXLIB=$(CPLEX)/cplex/lib/x86-64_osx/static_pic/
CONCERTLIB=$(CPLEX)/concert/lib/x86-64_osx/static_pic/

RM = rm -rf
OBJDIR = obj

BASEDIRS = \
	-I. \
	-I$(CPLEX)/cplex/include \
	-I$(CPLEX)/concert/include \
    -I/opt/homebrew/include 


INCLUDEDIRS = $(BASEDIRS)

CXX = /opt/homebrew/opt/llvm/bin/clang++
LD =
SUFFIX = o

CXXFLAGS = -Wall -pedantic -g $(INCLUDEDIRS) -m64 -fPIC \
	-fexceptions -DNEBUG -DIL_STD -Wno-long-long \
	-Wno-attributes -Wno-ignored-attributes -fpermissive -Wno-sign-compare \
    -Wno-unused-private-field


LDFLAGS = -L/opt/homebrew/lib -lboost_system  -lboost_filesystem -lboost_thread-mt \
	-L$(CPLEXLIB) -lilocplex -lcplex \
	-L$(CONCERTLIB) -lconcert -lm -fopenmp


RENAME = -o

CXXSOURCES = \
	STV.cpp \
	sim_stv.cpp \
	model.cpp \
	tree_stv.cpp \
	stv_distance.cpp 
	
CXXOBJECTS = $(patsubst %.cpp, $(OBJDIR)/%.$(SUFFIX), $(CXXSOURCES))

all : $(PROGRAM)

$(PROGRAM) : $(CXXOBJECTS)
	$(CXX) -o ${@} $(CXXOBJECTS) $(LD) $(LDFLAGS) 

$(OBJDIR)/%.$(SUFFIX) : %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(RENAME) $(@D)/$(@F) -c $(<)

clean:
	$(RM) $(CXXOBJECTS) $(PROGRAM) $(OBJDIR)


