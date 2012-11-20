
all : graph_stats comp-exp-rate


CFLAGS=-Wall -Wextra -g -O3
CXXFLAGS=-Wall -Wextra -g -O3

graph_stats : graph_stats.o
	g++ $(CXXFLAGS) -o $@ $^

%.o : %.cc
	g++ -c $(CXXFLAGS) -o $@ $^

%.o : %.c
	gcc -c $(CFLAGS) -o $@ $^

comp-exp-rate : comp-exp-rate.o rng.o
	gcc $(CFLAGS) -o $@ $^

clean :
	rm -f graph_stats comp-exp-rate *.o

