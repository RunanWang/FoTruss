//
// Created by 王润安 on 2024/1/2.
//

#ifndef MLG_CS_MLGWITHSCHEMA_H
#define MLG_CS_MLGWITHSCHEMA_H


#include <filesystem>
#include "Graph.h"
#include "map"
#include "unordered_map"
#include "vector"
#include "fstream"
#include "../Utils/logging.h"

struct edgeWithLayer {
    NODE_TYPE s{0};
    NODE_TYPE t{0};
    LAYER_TYPE l{0};

    edgeWithLayer() = default;

    edgeWithLayer(NODE_TYPE s_, NODE_TYPE t_, LAYER_TYPE l_) : s(s_), t(t_), l(l_) {};

    bool operator<(const edgeWithLayer &e) const {
        return s < e.s || (s == e.s && t < e.t);
    }
};

class Schema : public Graph {
public:
    NODE_TYPE *fromNodes{nullptr};
    bool **edgeInSchema{nullptr};

    ~Schema();

    EDGE_TYPE edgeToId(NODE_TYPE u, NODE_TYPE v);

    void idToEdge(EDGE_TYPE edgeId, NODE_TYPE *u, NODE_TYPE *v);

    void BuildFromEdgeLst(vector<edgeWithLayer> &edge_buf, NODE_TYPE num_of_vtx, EDGE_TYPE num_of_edge,
                          LAYER_TYPE num_of_layer);

    void buildFromCompact(Schema *sparseSchema, const EDGE_TYPE *trussness, LAYER_TYPE num_of_layer);

    NODE_TYPE getDegreeOfVertex(NODE_TYPE v) const;

    EDGE_TYPE start(NODE_TYPE v) const;
};

class MLGWithSchema {
public:
    MLGWithSchema() = default;

    ~MLGWithSchema();

    void LoadFromFile(const filesystem::path &input_path);

    void PrintStatistics();

    bool hasEdge(LAYER_TYPE layer, NODE_TYPE v1, NODE_TYPE v2) const;

    void getEdgeSchema();

    void genTriangle(EDGE_TYPE *triangle);

    void genTriangle(EDGE_TYPE *triangle, const EDGE_TYPE *edgeInSubgraph);

    double getMemUsage() const;

    Graph *graph_layers{nullptr};
    Graph *sumGraph{nullptr};
    map<LAYER_TYPE, LAYER_TYPE> layerNameToLayerId;
    LAYER_TYPE *layerIdToLayerName{nullptr};
    LAYER_TYPE *layerOrder{nullptr};

    LAYER_TYPE layersNum{0};
    NODE_TYPE nodesNum{0};
    NODE_TYPE maxNodesNum{0};

    Schema *schema{nullptr};
};


#endif //MLG_CS_MLGWITHSCHEMA_H
