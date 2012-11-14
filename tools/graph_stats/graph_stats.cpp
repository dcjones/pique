
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <vector>
#include <map>
#include <cstdio>

typedef std::pair<unsigned int, unsigned int> edge_t;

typedef boost::compressed_sparse_row_graph<
                    boost::bidirectionalS,
                    //boost::directedS,
                    void,
                    unsigned int,
                    boost::no_property,
                    unsigned int,
                    unsigned int> graph_t;

typedef std::map<unsigned int, unsigned int> component_map_t;
typedef boost::associative_property_map<component_map_t> component_prop_map_t;

void mm_parse_error()
{
    fprintf(stderr, "Error: incorrectly formatted matrix market file.\n");
    exit(1);
}

int main(int argc, char* argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: graph_stats adjmat.mm\n");
        return 1;
    }

    FILE* f = fopen(argv[1], "r");
    if (f == NULL) {
        fprintf(stderr, "Error: can't open %s.\n", argv[1]);
        return 1;
    }

    /* Note: this is not a parser for matrix market files, but
     * specifically the matrices produced by pique. */
    std::vector<edge_t> edges;
    char buffer[512];

    /* Read header */
    fgets(buffer, sizeof(buffer), f);
    if (strcmp(buffer, "%%MatrixMarket matrix coordinate integer general\n") != 0) {
        mm_parse_error();
    }

    unsigned int n, n2, m;
    fscanf(f, "%u %u %u\n", &n, &n2, &m);

    unsigned int u, v, w;
    while ((fgets(buffer, sizeof(buffer), f)) > 0) {
        sscanf(buffer, "%u %u %u\n", &u, &v, &w);
        edges.push_back(std::make_pair(u - 1, v - 1));
        edges.push_back(std::make_pair(v - 1, u - 1));
    }
    fclose(f);

    graph_t G(boost::edges_are_unsorted_multi_pass, edges.begin(), edges.end(), n);

    component_map_t components;
    boost::connected_components(G, component_prop_map_t(components));
    std::map<unsigned int, unsigned int> sizes;

    component_map_t::iterator i;
    std::map<unsigned int, unsigned int>::iterator j;
    for (i = components.begin(); i != components.end(); ++i) {
        j = sizes.find(i->second);
        if (j == sizes.end()) {
            sizes[i->second] = 1;
        }
        else j->second++;
    }

    for (j = sizes.begin(); j != sizes.end(); ++j) {
        printf("%u\t%u\n", j->first, j->second);
    }

    return 0;
}


