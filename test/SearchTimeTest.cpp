//
// Created by 王润安 on 2024/4/4.
//

#include "testHeader.h"
#include "../Algorithm/TrussOpt.h"
#include "../Algorithm/Algo.h"
#include "../Utils/logging.h"
#include "../Utils/Timer.h"
#include "../constant.h"
#include "algorithm"
#include "queue"

bool triangleConnected(NODE_TYPE node, const set<NODE_TYPE> &subgraph, MLGWithSchema *mlg) {
    NODE_TYPE u, v, w;
    vector<NODE_TYPE> subgraphList;
    subgraphList.assign(subgraph.begin(), subgraph.end());
    for (auto i = 0; i < subgraphList.size(); i++) {
        u = subgraphList[i];
        auto u_start = mlg->sumGraph->start(u);
        auto u_end = mlg->sumGraph->end(u);
        auto uList = &mlg->sumGraph->edges[u_start];
        auto uLen = u_end - u_start;
        for (auto j = i + 1; j < subgraphList.size(); j++) {
            v = subgraphList[j];
            auto v_start = mlg->sumGraph->start(v);
            auto v_end = mlg->sumGraph->end(v);
            auto vList = &mlg->sumGraph->edges[v_start];
            auto vLen = v_end - v_start;
            auto commonNeighbor = new NODE_TYPE[max(uLen, vLen)];
            NODE_TYPE commonNeighborLen = intersect(uList, vList, uLen, vLen, commonNeighbor);
            for (NODE_TYPE neighbor_i = 0; neighbor_i < commonNeighborLen; neighbor_i++) {
                w = commonNeighbor[neighbor_i];
                if (w == node) {
                    return true;
                }
            }
            delete[]commonNeighbor;
        }
    }
    return false;
}


void nodesTriangleConnectedByGivenEdgeBasic(MLGWithSchema *mlg, NODE_TYPE u, NODE_TYPE v, vector<NODE_TYPE> &nbs) {
    // 最基础的连通
    NODE_TYPE w;
    auto u_start = mlg->sumGraph->start(u);
    auto u_end = mlg->sumGraph->end(u);
    auto uList = &mlg->sumGraph->edges[u_start];
    auto uLen = u_end - u_start;
    auto v_start = mlg->sumGraph->start(v);
    auto v_end = mlg->sumGraph->end(v);
    auto vList = &mlg->sumGraph->edges[v_start];
    auto vLen = v_end - v_start;
    auto commonNeighbor = new NODE_TYPE[max(uLen, vLen)];
    NODE_TYPE commonNeighborLen = intersect(uList, vList, uLen, vLen, commonNeighbor);
    for (NODE_TYPE neighbor_i = 0; neighbor_i < commonNeighborLen; neighbor_i++) {
        w = commonNeighbor[neighbor_i];
        nbs.emplace_back(w);
    }
    delete[] commonNeighbor;
}


void getCommunityForGivenCombo(NODE_TYPE queryNode, MLGWithSchema *mlg, EDGE_TYPE *trussness,
                               vector<vector<NODE_TYPE>> &ans) {
    NODE_TYPE u, v, w;
    EDGE_TYPE e, e1, e2;
    auto neighbors = new NODE_TYPE[mlg->sumGraph->maxDegree];
    NODE_TYPE neighborSize = 0;
    set<NODE_TYPE> trussNodes;
    vector<EDGE_TYPE> trussEdges;
    auto queryNodeNeighborStart = mlg->sumGraph->start(queryNode);
    auto queryNodeNeighborEnd = mlg->sumGraph->end(queryNode);
    ans.clear();
    vector<NODE_TYPE> tempCommunity;
    // 先找到对应的最大的k
    EDGE_TYPE trussnessBound = 0;
    for (EDGE_TYPE neighbor_i = queryNodeNeighborStart; neighbor_i < queryNodeNeighborEnd; neighbor_i++) {
        auto neighbor = mlg->sumGraph->edges[neighbor_i];
        auto edge = mlg->schema->edgeToId(queryNode, neighbor);
        trussnessBound = max(trussness[edge], trussnessBound);
    }
    if (trussnessBound == 0) return;
    set<EDGE_TYPE> edgesAttachToQueryNode;

    // 然后找到所有满足k的社区
    for (EDGE_TYPE neighbor_i = queryNodeNeighborStart; neighbor_i < queryNodeNeighborEnd; neighbor_i++) {
        // 从一个连接查询点的边作为出发
        auto neighbor = mlg->sumGraph->edges[neighbor_i];
        auto edge = mlg->schema->edgeToId(queryNode, neighbor);
        if (trussness[edge] < trussnessBound or edgesAttachToQueryNode.find(edge) != edgesAttachToQueryNode.end()) {
            continue;
        }
        edgesAttachToQueryNode.insert(edge);
        set<EDGE_TYPE> currentEdgesInC;
        queue<EDGE_TYPE> toInsertEdges;
        toInsertEdges.push(edge);
        while (not toInsertEdges.empty()) {
            e = toInsertEdges.front();
            toInsertEdges.pop();
            if (currentEdgesInC.find(e) != currentEdgesInC.end()) {
                continue;
            }
            edgesTriangleConnectedByGivenEdgeBasic(mlg, e, neighbors, neighborSize);
            mlg->schema->idToEdge(e, &u, &v);
            for (NODE_TYPE j = 0; j < neighborSize; j++) {
                w = neighbors[j];
                e1 = mlg->schema->edgeToId(u, w);
                e2 = mlg->schema->edgeToId(v, w);
                if (trussness[e1] >= trussnessBound and trussness[e2] >= trussnessBound) {
                    toInsertEdges.push(e1);
                    toInsertEdges.push(e2);
                }
                if (w == queryNode) {
                    edgesAttachToQueryNode.insert(e1);
                    edgesAttachToQueryNode.insert(e2);
                }
            }
            currentEdgesInC.insert(e);
        }
        ans.resize(ans.size() + 1);
        // edges to nodes
        set<NODE_TYPE> nodeInCommunity;
        nodeInCommunity.insert(queryNode);
        for (auto ee: currentEdgesInC) {
            mlg->schema->idToEdge(ee, &u, &v);
            nodeInCommunity.insert(u);
            nodeInCommunity.insert(v);
        }
        ans[ans.size() - 1].assign(nodeInCommunity.begin(), nodeInCommunity.end());
    }
}

void getCommunity(NODE_TYPE queryNode, MLGWithSchema *mlg, EDGE_TYPE **trussness, ComboIndex *ci, double &bestden,
                  Timer *timer) {

    bestden = 0;
    vector<NODE_TYPE> community;
    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        // 循环全部可能的combo，找到其中density最大的一个
        if (ci->tomb[i]) {
            continue;
        }
        vector<vector<NODE_TYPE>> qus;
        auto queryNodeNeighborStart = mlg->sumGraph->start(queryNode);
        auto queryNodeNeighborEnd = mlg->sumGraph->end(queryNode);
        EDGE_TYPE trussnessBound = 0;
        for (EDGE_TYPE neighbor_i = queryNodeNeighborStart; neighbor_i < queryNodeNeighborEnd; neighbor_i++) {
            auto neighbor = mlg->sumGraph->edges[neighbor_i];
            auto edge = mlg->schema->edgeToId(queryNode, neighbor);
            trussnessBound = max(trussness[i][edge], trussnessBound);
        }
        timer->startTimer();
        getCommunityForGivenCombo(queryNode, mlg, trussness[i], qus);
        timer->endTimer();
        if (not qus.empty()) {
            LOG(INFO) << fmt::format("Found all trusses for {} in combo {}.", queryNode, ci->idToComboName[i]);
        }
        for (const auto &qu: qus) {
            double den = density(qu, mlg, DEN_BETA);
//            LOG(INFO) << fmt::format("Combo = {}-{}, den = {}, community = [{}]", ci->idToComboName[i], trussnessBound,
//                                     den, fmt::join(qu.begin(), qu.end(), ", "));
            if (den > bestden) {
                bestden = den;
            }
        }
    }
}


void searchCommByTrussnessForGivenQueries(const filesystem::path &input_dir_path, const string &dataset_name,
                                          const vector<NODE_TYPE> &queries) {
    auto path = filesystem::path(input_dir_path);
    auto graph_path = path.append(fmt::format("{}.txt", dataset_name));
    auto mlg = new MLGWithSchema();
    mlg->LoadFromFile(graph_path);
    auto maxComboSize = LAYER_TYPE(pow(2, mlg->layersNum)) + mlg->layersNum * LAYER_TYPE(pow(2, mlg->layersNum - 1));
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *ci = new ComboIndex();
    FoTrussDecompositionAll(mlg, trussness, ci);
    double sumDen = 0.0;
    double sumDiam = 0.0;
    double den, diam;
    auto timer = new Timer();
    LOG(INFO) << fmt::format("Finding Communities for q in Queries = [{}]", fmt::join(queries, ", "));
    for (auto q: queries) {
        LOG(INFO) << fmt::format("Finding community for query vertex {}...", q);
        getCommunity(q, mlg, trussness, ci, den, timer);
        sumDen += den;
    }
    double avgDen = sumDen / double(queries.size());
    double avgDiam = sumDiam / double(queries.size());
    double avgQueryTime = timer->getTimerSecond() / double(queries.size());
    LOG(INFO) << fmt::format("Avg Density of {} is {:.4f}, Diameter is {}, cost time = {:.4f}s.", dataset_name, avgDen,
                             avgDiam, avgQueryTime);
    delete ci;
}