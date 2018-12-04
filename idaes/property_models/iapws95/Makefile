ASL=$(ASL_BUILD)
BOOST=$(BOOST_HEADER)

CC = gcc
CXX = g++
RM = rm
CP = cp

ifeq ($(ASL),)
    CFLAGS = -Wall -Wextra -O3
else
    CFLAGS = -Wall -Wextra -O3 -I$(ASL)
endif

ifeq ($(BOOST),)
    CXXFLAGS = -c $(CFLAGS) -std=c++11 -fPIC
else
    CXXFLAGS = -c $(CFLAGS) -std=c++11 -fPIC -I$(BOOST)
endif

LDFLAGS = -shared -lm
LDFLAGS_EXE = -lm

OBJECTS = iapws95.o iapws95_phi.o iapws95_asl_funcs.o iapws95_memo.o

ALL: iapws95.so iapws95_tests

iapws95.o: iapws95.cpp
		$(CXX) $(CXXFLAGS) iapws95.cpp -o iapws95.o

iapws95_phi.o: iapws95_phi.cpp
	  $(CXX) $(CXXFLAGS) iapws95_phi.cpp -o iapws95_phi.o

iapws95_memo.o: iapws95_memo.cpp
	  $(CXX) $(CXXFLAGS) iapws95_memo.cpp -o iapws95_memo.o

iapws95_asl_funcs.o: iapws95_asl_funcs.cpp
	  $(CXX) $(CXXFLAGS) iapws95_asl_funcs.cpp -o iapws95_asl_funcs.o

iapws95_tests.o: iapws95_tests.cpp
		$(CXX) $(CXXFLAGS) iapws95_tests.cpp -o iapws95_tests.o

iapws95.so: $(OBJECTS)
	  $(CXX) $(LDFLAGS) $(OBJECTS) -o iapws95.so

iapws95_tests: iapws95_tests.o $(OBJECTS)
		$(CXX) $(LDFLAGS_EXE) iapws95_tests.o $(OBJECTS) -o iapws95_tests

clean:
	rm -f *.o
	rm -f *.so
	rm -f iapws95_tests
