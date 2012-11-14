
all : graph_stats comp-exp-rate


CFLAGS=-Wall -Wextra -g -O3
CXXFLAGS=-Wall -Wextra -g -O3

graph_stats : graph_stats.o
	clang++ $(CXXFLAGS) -o $@ $^

%.o : %.cc
	clang++ -c $(CXXFLAGS) -o $@ $^

%.o : %.c
	clang -c $(CFLAGS) -o $@ $^

comp-exp-rate : comp-exp-rate.o
	clang $(CFLAGS) -o $@ $^

clean :
	rm -f graph_stats comp-exp-rate *.o

