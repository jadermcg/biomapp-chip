include Makevars

CXXFLAGS += -I ../utils
UTILS = ../utils

all: em oops zoops

em: em.cpp oops.cpp zoops.cpp oops.h zoops.h em_utils.cpp em_utils.h $(UTILS)/utils.cpp $(UTILS)/prob_utils.cpp $(UTILS)/utils.h $(UTILS)/prob_utils.h
	$(CXX) $(CXXFLAGS) -o em em.cpp oops.cpp zoops.cpp em_utils.cpp $(UTILS)/utils.cpp $(UTILS)/prob_utils.cpp $(LIBS)

oops: run_oops.cpp oops.h oops.cpp $(UTILS)/utils.cpp $(UTILS)/prob_utils.cpp $(UTILS)/utils.h $(UTILS)/prob_utils.h 
	$(CXX) $(CXXFLAGS) -o oops run_oops.cpp oops.cpp $(UTILS)/utils.cpp $(UTILS)/prob_utils.cpp $(LIBS)

zoops: run_zoops.cpp zoops.h zoops.cpp $(UTILS)/utils.cpp $(UTILS)/prob_utils.cpp $(UTILS)/utils.h $(UTILS)/prob_utils.h
	$(CXX) $(CXXFLAGS) -o zoops run_zoops.cpp zoops.cpp $(UTILS)/utils.cpp $(UTILS)/prob_utils.cpp $(LIBS)

clean:
	rm -f em oops zoops *.o