CXX=g++
CXXFLAGS=-fmax-errors=3 -O3 -march=native
#CXXFLAGS=-g -fmax-errors=3 -Og -march=native -Wall
LIBSEIGEN= -I${EIGENHOME}
LIBSMOSEK= -I${MSKHOME}/h -L${MSKHOME}/bin -Wl,-rpath-link,${MSKHOME}/bin -Wl,-rpath=${MSKHOME}/bin -lmosek64 -lfusion64
LIBS=$(LIBSEIGEN) $(LIBSMOSEK)
MAIN=src/main.cpp
FILES=src/mon.cpp src/poly.cpp src/printing.cpp src/utils.cpp src/mosek.cpp
ASOBJ=$(FILES:.cpp=.o)

all: run

run: $(MAIN) $(ASOBJ)
	$(CXX) $(CXXFLAGS) -o run $(MAIN) $(ASOBJ) $(LIBS)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LIBS)

clean:
	rm -f run src/*.o

