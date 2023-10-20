include Makevars

all: smt hmap khmap kdive hsib smt ksearch dsearch main

main: main.cpp smt.cpp smt.h smt_utils.cpp smt_utils.h
	$(CXX) $(CXXFLAGS) -o main main.cpp smt.cpp smt_utils.cpp $(LDFLAGS) $(LIBS)

dsearch: dsearch.cpp smt_operations.cpp smt_operations.h smt_utils.h smt_utils.cpp
	$(CXX) $(CXXFLAGS) -o dsearch dsearch.cpp smt_operations.cpp smt_utils.cpp $(LDFLAGS) $(LIBS)

ksearch: ksearch.cpp smt_operations.cpp smt_operations.h smt_utils.h smt_utils.cpp
	$(CXX) $(CXXFLAGS) -o ksearch ksearch.cpp smt_operations.cpp smt_utils.cpp $(LDFLAGS) $(LIBS)

smt: run_smt.cpp smt.cpp smt_utils.cpp smt.h smt_utils.h
	$(CXX) $(CXXFLAGS) -o smt run_smt.cpp smt.cpp smt_utils.cpp $(LDFLAGS) $(LIBS)
	
hmap: run_hmap.cpp hmap.cpp smt_utils.cpp hmap.h smt_utils.h
	$(CXX) $(CXXFLAGS) -o hmap run_hmap.cpp hmap.cpp smt_utils.cpp $(LDFLAGS) $(LIBS)

khmap: khmap.cpp smt_operations.cpp smt_utils.cpp smt_operations.h smt_utils.h
	$(CXX) $(CXXFLAGS) -o khmap khmap.cpp smt_operations.cpp smt_utils.cpp $(LDFLAGS) $(LIBS)
	
kdive: kdive.cpp smt_operations.cpp smt_utils.cpp smt_operations.h smt_utils.h
	$(CXX) $(CXXFLAGS) -o kdive kdive.cpp smt_operations.cpp smt_utils.cpp $(LDFLAGS) $(LIBS)

hsib: hsib.cpp smt_operations.cpp smt_utils.cpp smt_operations.h smt_utils.h
	$(CXX) $(CXXFLAGS) -o hsib hsib.cpp smt_operations.cpp smt_utils.cpp $(LDFLAGS) $(LIBS)

clean:
	rm -f smt hmap khmap kdive hsib ksearch dsearch main *.o
