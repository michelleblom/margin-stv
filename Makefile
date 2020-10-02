PROGRAM = stv
CPLEX129=/opt/ibm/ILOG/CPLEX_Studio129
CPLEX129LIB=$(CPLEX129)/cplex/lib/x86-64_linux/static_pic/
CONCERT129LIB=$(CPLEX129)/concert/lib/x86-64_linux/static_pic/

RM = rm -rf
OBJDIR = obj

BASEDIRS = \
	-I. \
	-I$(CPLEX129)/cplex/include \
	-I$(CPLEX129)/concert/include 


INCLUDEDIRS = $(BASEDIRS)

CXX = g++
LD =
SUFFIX = o

CXXFLAGS = -Wall -pedantic -g $(INCLUDEDIRS) -m64 -fPIC \
	-fexceptions -DNEBUG -DIL_STD -Wno-long-long \
	-Wno-attributes -Wno-ignored-attributes -fpermissive -Wno-sign-compare


LDFLAGS =   -lboost_system  -lboost_filesystem -lboost_thread \
	-L$(CPLEX129LIB) -lilocplex -lcplex \
	-L$(CONCERT129LIB) -lconcert -m64 -lm -lpthread -ldl -lrt


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


