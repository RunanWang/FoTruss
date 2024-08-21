//
// Created by 王润安 on 2023/11/16.
//

#include "Graph.h"
#include "algorithm"
#include "../Utils/logging.h"

Graph::~Graph() {
    delete[] index;
    delete[] edges;
}

// Using CSR to store graph.
void Graph::BuildFromEdgeLst(vector<edge> &edge_buf, NODE_TYPE num_of_vtx, EDGE_TYPE num_of_edge) {
    EDGE_TYPE i, j, pNodeFirst;
    NODE_TYPE v;

    nodesNum = num_of_vtx;
    index = new EDGE_TYPE[nodesNum + 1];
    edges = new NODE_TYPE[num_of_edge + nodesNum + 1];
    // edges[index[v]:index[v+1]] = [deg_of_v(index ptr to here), nbr1, nbr2, ...]

    std::sort(edge_buf.begin(), edge_buf.end());

    i = 0; // index of edge_buf
    j = 0; // index of edges in CSR
    for (v = 0; v < nodesNum; v++) {
        index[v] = j;
        if (i >= num_of_edge || edge_buf[i].s > v) {
            // 这个点没有邻居，存一个度数=0
            edges[j++] = 0;
        } else {
            pNodeFirst = j; // index to store vtx degree
            j++;
            edges[j] = edge_buf[i].t;
            j++;
            i++;
            while (i < num_of_edge && edge_buf[i].s == v) {
                if (edge_buf[i].t == edge_buf[i - 1].t || edge_buf[i].s == edge_buf[i].t) {
                    i++;  // duplicated pr self edge
                } else {
                    edges[j] = edge_buf[i].t;
                    j++;
                    i++;
                }
            }
            edges[pNodeFirst] = j - pNodeFirst - 1;
            edgesNum += edges[pNodeFirst];
            if (edges[pNodeFirst] > maxDegree) maxDegree = edges[pNodeFirst];
        }
    }
    // 最后edges和index分别+1，方便访问时避免数组越界
    edges[j] = 0;
    index[nodesNum] = j;
}

NODE_TYPE Graph::getDegreeOfVertex(NODE_TYPE v) const {
    return edges[index[v]];
}

EDGE_TYPE Graph::start(NODE_TYPE v) const {
    return index[v] + 1;
}

EDGE_TYPE Graph::end(NODE_TYPE v) const {
    return index[v + 1];
}

EDGE_TYPE Graph::posLargerThanV(NODE_TYPE v) const {
    auto start_id = start(v);
    if (edges[start_id] > v) return start_id;
    auto end_id = end(v) - 1;
    if (start_id == end_id) return start_id;
    while (edges[start_id] < v) {
        auto mid_id = (start_id + end_id) / 2;
        if (mid_id == start_id or mid_id == end_id) break;
        if (edges[mid_id] > v) end_id = mid_id;
        else start_id = mid_id;
    }
    if (edges[start_id] > v) return start_id;
    else if (edges[end_id] > v) return end_id;
    else return end_id + 1;
}

bool Graph::hasEdge(NODE_TYPE u, NODE_TYPE v) const {
    auto fromNode = edges[index[u]] > edges[index[v]] ? v : u;
    auto toNode = edges[index[u]] > edges[index[v]] ? u : v;
    auto start_id = start(fromNode);
    auto end_id = end(fromNode) - 1;
    if (start_id == end_id) return false;
    while (start_id < end_id) {
        auto mid_id = (start_id + end_id) / 2;
        if (mid_id == start_id or mid_id == end_id) break;
        if (edges[mid_id] == toNode) return true;
        else if (edges[mid_id] > toNode) {
            end_id = mid_id;
        } else {
            // edges[mid_id] < toNode
            start_id = mid_id;
        }
    }
    if (edges[start_id] == toNode or edges[end_id] == toNode) return true;
    else return false;
}

double Graph::getMemUsage() const {
    double memInMB = 0;
    memInMB += double((nodesNum + 1) * sizeof(EDGE_TYPE)) / 1000.0 / 1000.0;                // index
    memInMB += double((edgesNum + nodesNum + 1) * sizeof(EDGE_TYPE)) / 1000.0 / 1000.0;     // nbs
    return memInMB;
}