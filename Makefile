#This should point to the location of the ITensor library
#ITENSOR_LIBRARY_DIR=$(HOME)/software/itensor
ITENSOR_LIBRARY_DIR=$(ITENSOR)
#These are sample MPI flags and should be adjusted
#for the specific version and location of MPI on your computer
#With some MPI distributions, you can use the command
#  mpicxx --showme:compile
#  mpicxx --showme:link
#to see the compiler flags (CFLAGS and LFLAGS respectively) for your MPI installation 
LIBRARY_DIR=$(ITENSOR)
include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk
MPICOM=mpicxx -m64 -std=c++17 -fconcepts -fPIC
MPI_CFLAGS=
#-I/share/opt/pgi/18.7/linux86-64/2018/mpi/openmpi/include/
MPI_LFLAGS=
#-I/share/opt/pgi/18.7/linux86-64/2018/mpi/openmpi/include/ -L/usr/local/opt/libevent/lib -L//share/opt/pgi/18.7/linux86-64/2018/mpi/openmpi/lib/ -lmpi_cxx -lmpi

ifdef app
APP=$(app)
else
APP=pdmrg
endif

HEADERS=parallel_tdvp.h partition.h

OBJECTS=$(APP).o

include $(ITENSOR_LIBRARY_DIR)/this_dir.mk
include $(ITENSOR_LIBRARY_DIR)/options.mk

#Mappings --------------
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

#Rules ------------------


%.o: %.cc $(HEADERS)
	$(MPICOM) $(MPI_CFLAGS) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(HEADERS)
	$(MPICOM) $(MPI_CFLAGS) -c $(CCGFLAGS) -o $@ $<

%.o: %.cpp $(HEADERS)
	$(MPICOM) $(MPI_CFLAGS) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cpp $(HEADERS)
	$(MPICOM) $(MPI_CFLAGS) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

default: build

build: $(APP)
debug: $(APP)-g

$(APP): $(OBJECTS) 
	$(MPICOM) $(MPI_LFLAGS) $(CCFLAGS) $(OBJECTS) -o $(APP) $(LIBFLAGS)

$(APP)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(MPICOM) $(MPI_LFLAGS) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS)

new_f: new_f.o
	$(MPICOM) $(MPI_LFLAGS) $(CCFLAGS) new_f.o -o new_f $(LIBFLAGS)


mkdebugdir:
	mkdir -p .debug_objs

clean:
	rm -fr *.o $(APP) $(APP)-g .debug_objs
