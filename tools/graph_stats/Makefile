
CXXFLAGS=-Wall -Wextra -g -O3

graph_stats : graph_stats.o
	clang++ $(CXXFLAGS) -o $@ $^

%.o : %.cc
	clang++ -c $(CXXFLAGS) -o $@ $^

clean :
	rm -f graph_stats *.o

