//
// Created by 王润安 on 2023/11/16.
//

#ifndef MLGWORK_GRAPH_H
#define MLGWORK_GRAPH_H

#include "../constant.h"
#include "vector"

using namespace std;

struct edge {
    NODE_TYPE s{0};
    NODE_TYPE t{0};

    edge() = default;

    edge(NODE_TYPE s_, NODE_TYPE t_) : s(s_), t(t_) {};

    bool operator<(const edge &e) const {
        return s < e.s || (s == e.s && t < e.t);
    }
};


class Graph {
public:
    Graph() = default;

    ~Graph();

    NODE_TYPE *edges;
    EDGE_TYPE *index;

    NODE_TYPE nodesNum{0};
    EDGE_TYPE edgesNum{0};
    NODE_TYPE maxDegree{0};

    void BuildFromEdgeLst(vector<edge> &edge_buf, NODE_TYPE num_of_vtx, EDGE_TYPE num_of_edge);

    NODE_TYPE getDegreeOfVertex(NODE_TYPE v) const;

    EDGE_TYPE start(NODE_TYPE v) const;

    EDGE_TYPE end(NODE_TYPE v) const;

    EDGE_TYPE posLargerThanV(NODE_TYPE v) const;

    bool hasEdge(NODE_TYPE u, NODE_TYPE v) const;

    double getMemUsage() const;

};


#endif //MLGWORK_GRAPH_H
