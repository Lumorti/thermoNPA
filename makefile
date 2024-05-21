CXX=g++
CXXFLAGS=-fmax-errors=3 -O3 -march=native
#CXXFLAGS=-g -fmax-errors=3 -Og -march=native -Wall
LIBSGUROBI= -I${GUROBIHOME}/include -L${GUROBIHOME}/lib -lgurobi_c++ -lgurobi110
LIBSEIGEN= -I${EIGENHOME}
LIBSMOSEK= -I${MSKHOME}/h -L${MSKHOME}/bin -Wl,-rpath-link,${MSKHOME}/bin -Wl,-rpath=${MSKHOME}/bin -lmosek64 -lfusion64
LIBS=$(LIBSEIGEN) $(LIBSMOSEK) $(LIBSGUROBI)
EXTRAFLAGS= -lgurobi_c++ -lgurobi110 -lmosek64 -lfusion64
MAIN=src/main.cpp
FILES=src/mon.cpp src/poly.cpp src/printing.cpp src/utils.cpp src/mosek.cpp src/gurobi.cpp
ASOBJ=$(FILES:.cpp=.o)

all: run

run: $(MAIN) $(ASOBJ)
	$(CXX) $(CXXFLAGS) -o run $(MAIN) $(ASOBJ) $(LIBS) $(EXTRAFLAGS)

src/mosek.o: src/mosek.cpp src/mosek.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBSEIGEN) $(LIBSMOSEK)

src/gurobi.o: src/gurobi.cpp src/gurobi.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBSEIGEN) $(LIBSGUROBI)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBSEIGEN)

clean:
	rm -f run src/*.o

