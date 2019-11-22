include ../Makefile.in

ifeq ($(BOOST),)
    CXXFLAGS = -c $(CFLAGS) -std=c++11 -fPIC
else
    CXXFLAGS = -c $(CFLAGS) -std=c++11 -fPIC -I$(BOOST)
endif

OBJECTS = iapws95_external.o iapws95_phi.o iapws95_asl_funcs.o iapws95_memo.o

ALL: iapws95_external.so

iapws95_external.o: iapws95_external.cpp
		$(CXX) $(CXXFLAGS) iapws95_external.cpp -o iapws95_external.o

iapws95_phi.o: iapws95_phi.cpp
	  $(CXX) $(CXXFLAGS) iapws95_phi.cpp -o iapws95_phi.o

iapws95_memo.o: iapws95_memo.cpp
	  $(CXX) $(CXXFLAGS) iapws95_memo.cpp -o iapws95_memo.o

iapws95_asl_funcs.o: iapws95_asl_funcs.cpp
	  $(CXX) $(CXXFLAGS) iapws95_asl_funcs.cpp -o iapws95_asl_funcs.o

iapws95_tests.o: iapws95_tests.cpp
		$(CXX) $(CXXFLAGS) iapws95_tests.cpp -o iapws95_tests.o

iapws95_external.so: $(OBJECTS)
	  $(LINK) $(LDFLAGS) $(OBJECTS) -o iapws95_external.so

clean:
	rm -f *.o
	rm -f *.so
	rm -f iapws95_tests
