CXX=g++
CXXFLAGS=-fmax-errors=3 -O3 -march=native -fopenmp
#CXXFLAGS=-g -fmax-errors=3 -Wall
LIBSGUROBI= -I${GUROBIHOME}/include -L${GUROBIHOME}/lib
LIBSEIGEN= -I${EIGENHOME}
LIBSMOSEK= -I${MSKHOME}/h -L${MSKHOME}/bin -Wl,-rpath-link,${MSKHOME}/bin -Wl,-rpath=${MSKHOME}/bin
LIBS=$(LIBSEIGEN) $(LIBSMOSEK) $(LIBSGUROBI)
EXTRAFLAGS= -lgurobi_c++ -lgurobi110 -lmosek64 -lfusion64
FILES=$(wildcard src/*.cpp)
ASOBJ=$(FILES:.cpp=.o)

all: run

run: $(ASOBJ)
	$(CXX) $(CXXFLAGS) -o run $(ASOBJ) $(LIBS) $(EXTRAFLAGS)

src/optMOSEK.o: src/optMOSEK.cpp src/optMOSEK.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBSEIGEN) $(LIBSMOSEK)

src/optGurobi.o: src/optGurobi.cpp src/optGurobi.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBSEIGEN) $(LIBSGUROBI)

src/mon.o: src/mon.cpp src/mon.h src/settings.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

src/main.o: src/main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBSEIGEN)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBSEIGEN)

clean:
	rm -f run src/*.o

