GCC = g++ "-std=c++11" -O3
uclique: uc_main.o uncertain_clique.o
	$(GCC) -o uclique uncertain_clique.o uc_main.o
uc_main.o: uc_main.cpp
	$(GCC) -c uc_main.cpp
uncertain_clique.o: uncertain_clique.cpp uncertain_clique.h
	$(GCC) -c uncertain_clique.cpp

clean:
	rm -f *.o