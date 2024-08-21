//
// Created by 王润安 on 2024/1/2.
//

#include "MLGWithSchema.h"
#include "../Algorithm/Algo.h"
#include "../Utils/logging.h"
#include "algorithm"

MLGWithSchema::~MLGWithSchema() {
    delete[] graph_layers;
    delete[] layerOrder;
    delete[] layerIdToLayerName;
    delete sumGraph;
    delete schema;
    layerNameToLayerId.clear();
    map<LAYER_TYPE, LAYER_TYPE>().swap(layerNameToLayerId);
}

void MLGWithSchema::LoadFromFile(const filesystem::path &input_file) {
    vector<edge> *edge_buf;
    vector<edge> sum_buf;
    vector<edgeWithLayer> total_buf;
    EDGE_TYPE *edge_size;
    EDGE_TYPE total_size = 0;
    EDGE_TYPE sum_size = 0;

    LAYER_TYPE layer, layerId, li;
    NODE_TYPE fromNode, toNode;

    ifstream fin;
    fin.open(input_file);
    fin >> layersNum >> nodesNum >> maxNodesNum;
    nodesNum += 1;
    edge_buf = new vector<edge>[layersNum];
    edge_size = new EDGE_TYPE[layersNum];
    graph_layers = new Graph[layersNum];
    sumGraph = new Graph();
    schema = new Schema();
    layerIdToLayerName = new LAYER_TYPE[layersNum];
    for (layer = 0; layer < layersNum; layer++) {
//        edge_buf[layer] = *new vector<edge>();
        edge_size[layer] = 0;
    }
    li = 0;
    while (not fin.eof()) {
        fin >> layer >> fromNode >> toNode;
        if (layerNameToLayerId.find(layer) == layerNameToLayerId.end()) {
            layerNameToLayerId[layer] = li;
            layerIdToLayerName[li] = layer;
            li++;
            layerId = layerNameToLayerId[layer];
        } else {
            layerId = layerNameToLayerId[layer];
        }
        edge_size[layerId] += 2;
        edge_buf[layerId].emplace_back(fromNode, toNode);
        edge_buf[layerId].emplace_back(toNode, fromNode);
        sum_buf.emplace_back(fromNode, toNode);
        sum_buf.emplace_back(toNode, fromNode);
        sum_size += 2;
        if (fromNode < toNode) {
            total_buf.emplace_back(fromNode, toNode, layerId);
        } else {
            total_buf.emplace_back(toNode, fromNode, layerId);
        }
        total_size += 1;
    }

    LOG(INFO) << fmt::format("Converting EdgeLists Into CSR...");
    for (layer = 0; layer < layersNum; layer++) {
        graph_layers[layer].BuildFromEdgeLst(edge_buf[layer], nodesNum, edge_size[layer]);
    }
    schema->BuildFromEdgeLst(total_buf, nodesNum, total_size, layersNum);
    sumGraph->BuildFromEdgeLst(sum_buf, nodesNum, sum_size);

    LOG(INFO) << fmt::format("Cleaning Up...");
    delete[] edge_size;
    for (layer = 0; layer < layersNum; layer++) {
        edge_buf[layer].clear();
        edge_buf[layer].shrink_to_fit();
    }
    total_buf.clear();
    total_buf.shrink_to_fit();
    delete[] edge_buf;
}

void MLGWithSchema::PrintStatistics() {
    LAYER_TYPE layer;
    LOG(INFO) << fmt::format("Dataset Info >>>>>>");
    LOG(INFO) << fmt::format("Layers = {}.", layersNum);
    LOG(INFO) << fmt::format("Nodes = {}.", nodesNum);
    for (layer = 0; layer < layersNum; layer++) {
        auto edgeNum = graph_layers[layer].edgesNum / 2;
        double den = double(edgeNum) / nodesNum;
//        LOG(INFO) << fmt::format("On layer {}, edges = {}, density = {:.3f}.", layerIdToLayerName[layer], edgeNum, den);
    }

}

bool MLGWithSchema::hasEdge(LAYER_TYPE layer, NODE_TYPE v1, NODE_TYPE v2) const {
//    auto g = &graph_layers[layer];
//    auto fromNode = g->edges[g->index[v1]] > g->edges[g->index[v2]] ? v2 : v1;
//    auto toNode = g->edges[g->index[v1]] > g->edges[g->index[v2]] ? v1 : v2;
//    for (EDGE_TYPE i = g->start(fromNode); i < g->end(fromNode); i++) {
//        if (g->edges[i] == toNode) {
//            return true;
//        }
//    }
//    return false;
    return graph_layers[layer].hasEdge(v1, v2);
}

inline EDGE_TYPE calInd(EDGE_TYPE edge, LAYER_TYPE layer, LAYER_TYPE layerNum) {
    return edge * layerNum + layer;
}

void MLGWithSchema::genTriangle(EDGE_TYPE *triangle) {
    EDGE_TYPE edge = 0;
    LAYER_TYPE layer = 0;
    NODE_TYPE u, v;
    for (edge = 0; edge < schema->edgesNum; edge++) {
        for (layer = 0; layer < layersNum; layer++) {
            triangle[calInd(edge, layer, layersNum)] = 0;
        }
    }
    for (edge = 0; edge < schema->edgesNum - 1; edge++) {
        for (layer = 0; layer < layersNum; layer++) {
            if (not schema->edgeInSchema[edge][layer]) continue;
            schema->idToEdge(edge, &u, &v);
            auto u_start = graph_layers[layer].posLargerThanV(u);
            auto u_end = graph_layers[layer].end(u);
            auto uList = &graph_layers[layer].edges[u_start];
            auto uLen = u_end - u_start;
            auto v_start = graph_layers[layer].posLargerThanV(v);
            auto v_end = graph_layers[layer].end(v);
            auto vList = &graph_layers[layer].edges[v_start];
            auto vLen = v_end - v_start;
            auto commonNeighbor = new NODE_TYPE[max(uLen, vLen)];
            NODE_TYPE commonNeighborLen = 0;
            commonNeighborLen = intersect(uList, vList, uLen, vLen, commonNeighbor);
            triangle[calInd(edge, layer, layersNum)] += commonNeighborLen;
            for (NODE_TYPE neighbor_j = 0; neighbor_j < commonNeighborLen; neighbor_j++) {
                auto edge1 = schema->edgeToId(u, commonNeighbor[neighbor_j]);
                auto edge2 = schema->edgeToId(v, commonNeighbor[neighbor_j]);
                triangle[calInd(edge1, layer, layersNum)] += 1;
                triangle[calInd(edge2, layer, layersNum)] += 1;
            }
            delete[] commonNeighbor;
        }
    }
}

void MLGWithSchema::genTriangle(EDGE_TYPE *triangle, const EDGE_TYPE *edgeInSubgraph) {
    EDGE_TYPE edge = 0;
    LAYER_TYPE layer = 0;
    NODE_TYPE u, v;
    for (edge = 0; edge < schema->edgesNum; edge++) {
        for (layer = 0; layer < layersNum; layer++) {
            triangle[calInd(edge, layer, layersNum)] = 0;
        }
    }
    for (edge = 0; edge < schema->edgesNum - 1; edge++) {
        if (edgeInSubgraph[edge] == 0) continue;
        for (layer = 0; layer < layersNum; layer++) {
            if (not schema->edgeInSchema[edge][layer]) continue;
            schema->idToEdge(edge, &u, &v);
            auto u_start = graph_layers[layer].posLargerThanV(u);
            auto u_end = graph_layers[layer].end(u);
            auto uList = &graph_layers[layer].edges[u_start];
            auto uLen = u_end - u_start;
            auto v_start = graph_layers[layer].posLargerThanV(v);
            auto v_end = graph_layers[layer].end(v);
            auto vList = &graph_layers[layer].edges[v_start];
            auto vLen = v_end - v_start;
            auto commonNeighbor = new NODE_TYPE[max(uLen, vLen)];
            NODE_TYPE commonNeighborLen = 0;
            commonNeighborLen = intersect(uList, vList, uLen, vLen, commonNeighbor);
            for (NODE_TYPE neighbor_j = 0; neighbor_j < commonNeighborLen; neighbor_j++) {
                auto edge1 = schema->edgeToId(u, commonNeighbor[neighbor_j]);
                auto edge2 = schema->edgeToId(v, commonNeighbor[neighbor_j]);
                if (edgeInSubgraph[edge1] == 0 or edgeInSubgraph[edge2] == 0) {
                    continue;
                }
                triangle[calInd(edge, layer, layersNum)]++;
                triangle[calInd(edge1, layer, layersNum)] += 1;
                triangle[calInd(edge2, layer, layersNum)] += 1;
            }
            delete[] commonNeighbor;
        }
    }
}


double MLGWithSchema::getMemUsage() const {
    LAYER_TYPE layer;
    double memInMB = 0;
    // Graph of each layer
    for (layer = 0; layer < layersNum; layer++) {
        memInMB += graph_layers[layer].getMemUsage();
    }
    // sumgraph
    memInMB += sumGraph->getMemUsage();
    return memInMB;
}


void Schema::BuildFromEdgeLst(vector<edgeWithLayer> &edge_buf, NODE_TYPE num_of_vtx, EDGE_TYPE num_of_edge,
                              LAYER_TYPE num_of_layer) {
    EDGE_TYPE i, j, edgeNum;
    NODE_TYPE v;
    LAYER_TYPE l;

    nodesNum = num_of_vtx;
    index = new EDGE_TYPE[nodesNum + 1];
    edges = new NODE_TYPE[num_of_edge + 1];
    fromNodes = new NODE_TYPE[num_of_edge + 1];
    // edges[index[v]:index[v+1]] = [deg_of_v(index ptr to here), nbr1, nbr2, ...]

    std::sort(edge_buf.begin(), edge_buf.end());

    i = 0; // index of edge_buf
    j = 0; // index of edges in CSR
    for (v = 0; v < nodesNum; v++) {
        index[v] = j;
        if (i >= num_of_edge || edge_buf[i].s > v) {
            // 这个点没有邻居，存一个度数=0
            continue;
        } else {
//            j++;
            fromNodes[j] = edge_buf[i].s;
            edges[j] = edge_buf[i].t;
            j++;
            i++;
            while (i < num_of_edge && edge_buf[i].s == v) {
                if (edge_buf[i].t == edge_buf[i - 1].t || edge_buf[i].s == edge_buf[i].t) {
                    i++;  // duplicated pr self edge
                } else {
                    fromNodes[j] = edge_buf[i].s;
                    edges[j] = edge_buf[i].t;
                    j++;
                    i++;
                }
            }
            edgeNum = j - index[v];
            edgesNum += edgeNum;
            if (edgeNum > maxDegree) maxDegree = edgeNum;
        }
    }
    // 最后edges和index分别+1，方便访问时避免数组越界
    edges[j] = 0;
    index[nodesNum] = j;
    fromNodes[j] = 0;

    edgeInSchema = new bool *[edgesNum + 1];
    for (EDGE_TYPE edge = 0; edge < edgesNum; edge++) {
        edgeInSchema[edge] = new bool[num_of_layer];
        for (LAYER_TYPE tl = 0; tl < num_of_layer; tl++) {
            edgeInSchema[edge][tl] = false;
        }
    }
    for (auto edge: edge_buf) {
        EDGE_TYPE id = edgeToId(edge.s, edge.t);
        edgeInSchema[id][edge.l] = true;
    }
}

NODE_TYPE Schema::getDegreeOfVertex(NODE_TYPE v) const {
    return index[v + 1] - index[v];
}

EDGE_TYPE Schema::start(NODE_TYPE v) const {
    return index[v];
}

EDGE_TYPE Schema::edgeToId(NODE_TYPE u, NODE_TYPE v) {
    NODE_TYPE fromNode = min(u, v);
    NODE_TYPE toNode = max(u, v);
    auto start_id = start(fromNode);
    auto end_id = end(fromNode) - 1;
    while (start_id < end_id) {
        auto mid_id = (start_id + end_id) / 2;
        if (mid_id == start_id or mid_id == end_id) break;
        if (edges[mid_id] == toNode) return mid_id;
        else if (edges[mid_id] > toNode) {
            end_id = mid_id;
        } else {
            // edges[mid_id] < toNode
            start_id = mid_id;
        }
    }
    if (edges[start_id] == toNode) return start_id;
    else if (edges[end_id] == toNode) return end_id;
    else return 0;
}

void Schema::idToEdge(EDGE_TYPE edgeId, NODE_TYPE *u, NODE_TYPE *v) {
    *v = edges[edgeId];
    *u = fromNodes[edgeId];
}

Schema::~Schema() {
//    delete[] index;
//    delete[] edges;
    for (EDGE_TYPE edge = 0; edge < edgesNum; edge++) {
        delete[] edgeInSchema[edge];
    }
    delete[] edgeInSchema;
    delete[] fromNodes;
}

void Schema::buildFromCompact(Schema *sparseSchema, const EDGE_TYPE *trussness, LAYER_TYPE num_of_layer) {
    EDGE_TYPE i, j;
    NODE_TYPE v;
    LAYER_TYPE l;

    nodesNum = sparseSchema->nodesNum;
    index = new EDGE_TYPE[nodesNum + 1];
    edgesNum = 0;
    for (i = 0; i < sparseSchema->edgesNum; i++) {
        if (trussness[i] != 0) edgesNum++;
    }

    edges = new NODE_TYPE[edgesNum + 1];
    fromNodes = new NODE_TYPE[edgesNum + 1];
    edgeInSchema = new bool *[edgesNum + 1];
    for (EDGE_TYPE edge = 0; edge < edgesNum; edge++) {
        edgeInSchema[edge] = new bool[num_of_layer];
        for (LAYER_TYPE tl = 0; tl < num_of_layer; tl++) {
            edgeInSchema[edge][tl] = false;
        }
    }

    // 需要compact的：
    // NODE_TYPE *edges;
    // TODO EDGE_TYPE *index;
    // TODO NODE_TYPE *fromNodes{nullptr};
    // bool **edgeInSchema{nullptr};

    for (i = 0; i < nodesNum; i++) {
        index[i] = 0;
    }

    j = 0; // 压缩后的边id
    NODE_TYPE lastNodeId = nodesNum + 10;
    for (i = 0; i < sparseSchema->edgesNum; i++) {
        if (lastNodeId != sparseSchema->fromNodes[i]) {
            lastNodeId = sparseSchema->fromNodes[i];
            index[sparseSchema->fromNodes[i]] = j;
        }
        if (trussness[i] != 0) {
            edges[j] = sparseSchema->edges[i];
            fromNodes[j] = sparseSchema->fromNodes[i];
            for (LAYER_TYPE tl = 0; tl < num_of_layer; tl++) {
                edgeInSchema[j][tl] = sparseSchema->edgeInSchema[i][tl];
            }
            j++;
        }
    }

    // 可能有的点没有attach边，那么它的index没有更新，我们更新一下
    for (i = 1; i < nodesNum; i++) {
        if (index[i] == 0) {
            index[i] = index[i - 1];
        }
    }


}