//
// Created by 王润安 on 2024/1/2.
//

#include "TrussOpt.h"
#include "Algo.h"
#include "../DataStructure/Bucket.h"
#include "set"
#include "../Utils/logging.h"
#include "../Utils/Timer.h"
#include "vector"
#include "queue"
#include "../DataStructure/UnionFind.h"
#include "../constant.h"
#include "iostream"
#include "math.h"

using namespace std;

inline EDGE_TYPE calInd(EDGE_TYPE edge, LAYER_TYPE layer, LAYER_TYPE layerNum) {
    return edge * layerNum + layer;
}


EDGE_TYPE FoTrussDecompositionForGivenCombo(MLGWithSchema *mlg, EDGE_TYPE *trussness, EDGE_TYPE *initTriangles,
                                            EDGE_TYPE *fatherRes, const Combo &combo) {
    // Init
    LAYER_TYPE lamb = combo.lamb;
    LAYER_TYPE layer = 0;
    NODE_TYPE u, v, w;
    EDGE_TYPE edge, edge1, edge2, newK;
    EDGE_TYPE maxK = 0;
    EDGE_TYPE maxTrussness = 0;

    LOG(INFO) << fmt::format("Decomposition for combo {}", combo.toString());
    auto triangles = new EDGE_TYPE[mlg->schema->edgesNum * mlg->layersNum];    // 存每条边目前的三角形数目
    memcpy(triangles, initTriangles, mlg->schema->edgesNum * mlg->layersNum * sizeof(EDGE_TYPE));
//    mlg->genTriangle(triangles);
    auto subGraph = new bool[mlg->schema->edgesNum];                           // 存每条边目前是否在子图中
    auto toCheck = set<EDGE_TYPE>();                                        // 需要top-k检查的临时存储边set
    bool *layerIsFocus = new bool[mlg->layersNum];                          // 一个focus-layer的map,接下来对它初始化
    auto commonNeighbor = new NODE_TYPE[mlg->sumGraph->maxDegree * 2];
    for (layer = 0; layer < mlg->layersNum; layer++) {
        layerIsFocus[layer] = false;
    }
    for (layer = 0; layer < combo.focusNum; layer++) {
        layerIsFocus[combo.focusLayers[layer]] = true;
    }

    // triangles和trussness初始化
    for (edge = 0; edge < mlg->schema->edgesNum; edge++) {
        trussness[edge] = findTopK(&triangles[calInd(edge, 0, mlg->layersNum)], mlg->layersNum, lamb);
        trussness[edge] = min(trussness[edge], fatherRes[edge]);
        for (layer = 0; layer < combo.focusNum; layer++) {
            LAYER_TYPE focusLayer = combo.focusLayers[layer];
            trussness[edge] = min(trussness[edge], triangles[calInd(edge, focusLayer, mlg->layersNum)]);
        }
        maxK = max(maxK, trussness[edge]);
        subGraph[edge] = true;
    }
    auto bucket = new Bucket(mlg->schema->edgesNum, maxK);
    bucket->construct(trussness);

    for (NODE_TYPE k = 0; k < maxK + 1; k++) {
        auto toRemoveEdges = bucket->getBucket(k);
        if (bucket->length[k] != 0) {
            maxTrussness = k;
        }
        for (NODE_TYPE i = 0; i < bucket->length[k]; i++) {
            edge = toRemoveEdges[i];
            subGraph[edge] = false;
            trussness[edge] = k;
            for (layer = 0; layer < mlg->layersNum; layer++) {
                if (not mlg->schema->edgeInSchema[edge][layer]) continue;
                mlg->schema->idToEdge(edge, &u, &v);
                auto u_start = mlg->graph_layers[layer].start(u);
                auto u_end = mlg->graph_layers[layer].end(u);
                auto uList = &mlg->graph_layers[layer].edges[u_start];
                auto uLen = u_end - u_start;
                auto v_start = mlg->graph_layers[layer].start(v);
                auto v_end = mlg->graph_layers[layer].end(v);
                auto vList = &mlg->graph_layers[layer].edges[v_start];
                auto vLen = v_end - v_start;
                NODE_TYPE commonNeighborLen = intersect(uList, vList, uLen, vLen, commonNeighbor);
                for (NODE_TYPE neighbor_i = 0; neighbor_i < commonNeighborLen; neighbor_i++) {
                    w = commonNeighbor[neighbor_i];
                    edge1 = mlg->schema->edgeToId(u, w);
                    edge2 = mlg->schema->edgeToId(v, w);
                    if (subGraph[edge1] and subGraph[edge2]) {
                        // edge1
                        triangles[calInd(edge1, layer, mlg->layersNum)] -= 1;
                        newK = triangles[calInd(edge1, layer, mlg->layersNum)];
                        newK = min(newK, fatherRes[edge1]);
                        if (layerIsFocus[layer] and newK < trussness[edge1] and newK > k) {
                            bucket->changeItemBucket(edge1, trussness[edge1], newK);
                            trussness[edge1] = newK;
                        }
                        if (newK + 1 == trussness[edge1])
                            toCheck.insert(edge1);

                        // edge2
                        triangles[calInd(edge2, layer, mlg->layersNum)] -= 1;
                        newK = triangles[calInd(edge2, layer, mlg->layersNum)];
                        newK = min(newK, fatherRes[edge2]);
                        if (layerIsFocus[layer] and newK < trussness[edge2] and newK > k) {
                            bucket->changeItemBucket(edge2, trussness[edge2], newK);
                            trussness[edge2] = newK;
                        }
                        if (newK + 1 == trussness[edge2])
                            toCheck.insert(edge2);
                    }
                }
            }
            for (auto checkEdge: toCheck) {
//                if (not subGraph[checkEdge] or trussness[checkEdge] <= k) continue;
                auto newTopK = findTopK(&triangles[calInd(checkEdge, 0, mlg->layersNum)], mlg->layersNum, lamb);
                newTopK = min(newTopK, fatherRes[checkEdge]);
                if (newTopK != trussness[checkEdge] and newTopK <= k) {
                    bucket->changeItemBucket(checkEdge, trussness[checkEdge], k);
                    trussness[checkEdge] = k;
                } else if (newTopK < trussness[checkEdge]) {
                    bucket->changeItemBucket(checkEdge, trussness[checkEdge], newTopK);
                    trussness[checkEdge] = newTopK;
                } else {
                    continue;
                }
            }
            toCheck.clear();
        }
    }
    LOG(INFO) << fmt::format("max trussness = {}", maxTrussness + 2);
//    for (edge = 0; edge < mlg->edgeSchemaNum; edge++) {
//        if (subGraph[edge]) LOG(ERROR)
//                    << fmt::format("Error: after peeling, edge {} still in graph with trussness = {}.", edge,
//                                   trussness[edge]);
//    }
    delete[] triangles;
    delete[] subGraph;
    delete[] layerIsFocus;
    delete[] commonNeighbor;
    delete bucket;
    return maxTrussness;
}


EDGE_TYPE FoTrussDecompositionForGivenComboUsingFather(MLGWithSchema *mlg, EDGE_TYPE *trussness, EDGE_TYPE *fatherRes,
                                                       const Combo &combo) {
    // Init
    LAYER_TYPE lamb = combo.lamb;
    LAYER_TYPE layer = 0;
    NODE_TYPE u, v, w;
    EDGE_TYPE edge, edge1, edge2, newK;
    EDGE_TYPE maxK = 0;
    EDGE_TYPE maxTrussness = 0;

    LOG(INFO) << fmt::format("Decomposition for combo {}.", combo.toString());
    auto triangles = new EDGE_TYPE[mlg->schema->edgesNum * mlg->layersNum];    // 存每条边目前的三角形数目
    auto subGraph = new bool[mlg->schema->edgesNum];                           // 存每条边目前是否在子图中
    auto toCheck = set<EDGE_TYPE>();                                        // 需要top-k检查的临时存储边set
    auto commonNeighbor = new NODE_TYPE[mlg->sumGraph->maxDegree * 2];
    bool *layerIsFocus = new bool[mlg->layersNum];                          // 一个focus-layer的map,接下来对它初始化
    for (layer = 0; layer < mlg->layersNum; layer++) {
        layerIsFocus[layer] = false;
    }
    for (layer = 0; layer < combo.focusNum; layer++) {
        layerIsFocus[combo.focusLayers[layer]] = true;
    }

    // triangles和trussness初始化
    mlg->genTriangle(triangles, fatherRes);
    for (edge = 0; edge < mlg->schema->edgesNum; edge++) {
        trussness[edge] = findTopK(&triangles[calInd(edge, 0, mlg->layersNum)], mlg->layersNum, lamb);
        trussness[edge] = min(trussness[edge], fatherRes[edge]);
        for (layer = 0; layer < combo.focusNum; layer++) {
            LAYER_TYPE focusLayer = combo.focusLayers[layer];
            trussness[edge] = min(trussness[edge], triangles[calInd(edge, focusLayer, mlg->layersNum)]);
        }
        maxK = max(maxK, trussness[edge]);
        if (fatherRes[edge] == 0) { subGraph[edge] = false; }
        else { subGraph[edge] = true; }
    }
    auto bucket = new Bucket(mlg->schema->edgesNum, maxK);
    bucket->construct(trussness);

    for (NODE_TYPE k = 0; k < maxK + 1; k++) {
        auto toRemoveEdges = bucket->getBucket(k);
        if (bucket->length[k] != 0) {
            maxTrussness = k;
        }
        for (NODE_TYPE i = 0; i < bucket->length[k]; i++) {
            edge = toRemoveEdges[i];
            if (not subGraph[edge]) continue; // 这个edge被father已经移除了，它的trussness=0，并且triangle里没有计数，不用再peel
            subGraph[edge] = false;
//            trussness[edge] = k;
            for (layer = 0; layer < mlg->layersNum; layer++) {
                if (not mlg->schema->edgeInSchema[edge][layer]) continue;
                mlg->schema->idToEdge(edge, &u, &v);
                auto u_start = mlg->graph_layers[layer].start(u);
                auto u_end = mlg->graph_layers[layer].end(u);
                auto uList = &mlg->graph_layers[layer].edges[u_start];
                auto uLen = u_end - u_start;
                auto v_start = mlg->graph_layers[layer].start(v);
                auto v_end = mlg->graph_layers[layer].end(v);
                auto vList = &mlg->graph_layers[layer].edges[v_start];
                auto vLen = v_end - v_start;
                NODE_TYPE commonNeighborLen = intersect(uList, vList, uLen, vLen, commonNeighbor);
                for (NODE_TYPE neighbor_i = 0; neighbor_i < commonNeighborLen; neighbor_i++) {
                    w = commonNeighbor[neighbor_i];
                    edge1 = mlg->schema->edgeToId(u, w);
                    edge2 = mlg->schema->edgeToId(v, w);
                    if (subGraph[edge1] and subGraph[edge2]) {
                        // edge1
                        triangles[calInd(edge1, layer, mlg->layersNum)] -= 1;
                        newK = triangles[calInd(edge1, layer, mlg->layersNum)];
                        newK = min(newK, fatherRes[edge1]);
                        if (layerIsFocus[layer] and newK < trussness[edge1] and newK > k) {
                            bucket->changeItemBucket(edge1, trussness[edge1], newK);
                            trussness[edge1] = newK;
                        }
                        if (newK + 1 == trussness[edge1])
                            toCheck.insert(edge1);

                        // edge2
                        triangles[calInd(edge2, layer, mlg->layersNum)] -= 1;
                        newK = triangles[calInd(edge2, layer, mlg->layersNum)];
                        newK = min(newK, fatherRes[edge2]);
                        if (layerIsFocus[layer] and newK < trussness[edge2] and newK > k) {
                            bucket->changeItemBucket(edge2, trussness[edge2], newK);
                            trussness[edge2] = newK;
                        }
                        if (newK + 1 == trussness[edge2])
                            toCheck.insert(edge2);
                    }
                }
            }
            for (auto checkEdge: toCheck) {
//                if (not subGraph[checkEdge] or trussness[checkEdge] <= k) continue;
                auto newTopK = findTopK(&triangles[calInd(checkEdge, 0, mlg->layersNum)], mlg->layersNum, lamb);
                newTopK = min(newTopK, fatherRes[checkEdge]);
                if (newTopK != trussness[checkEdge] and newTopK <= k) {
                    bucket->changeItemBucket(checkEdge, trussness[checkEdge], k);
                    trussness[checkEdge] = k;
                } else if (newTopK < trussness[checkEdge]) {
                    bucket->changeItemBucket(checkEdge, trussness[checkEdge], newTopK);
                    trussness[checkEdge] = newTopK;
                } else {
                    continue;
                }
            }
            toCheck.clear();
        }
    }
    LOG(INFO) << fmt::format("max trussness = {}", maxTrussness + 2);
    delete[] triangles;
    delete[] subGraph;
    delete[] layerIsFocus;
    delete bucket;
    delete[] commonNeighbor;
    return maxTrussness;
}


void FoTrussDecompositionOpt(MLGWithSchema *mlg, const filesystem::path &outPath) {
    auto maxComboSize = LAYER_TYPE(pow(2, mlg->layersNum)) + mlg->layersNum * LAYER_TYPE(pow(2, mlg->layersNum - 1));
    auto result = new EDGE_TYPE *[maxComboSize];
    auto ci = new ComboIndex();
    FoTrussDecompositionAll(mlg, result, ci);
    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        delete[] result[i];
    }
    delete[] result;
    delete ci;
}


void FoTrussDecompositionBasic(MLGWithSchema *mlg, EDGE_TYPE **result, ComboIndex *comboIndex) {
    auto timer = new Timer();
    timer->startTimer();
    Combo *c;
    EDGE_TYPE maxTrussness = 0;                                     // 用于判断分解的truss是否为空
    EDGE_TYPE schemaNum = mlg->schema->edgesNum;                    // Edge-schema的总数
    LAYER_TYPE fatherSize;
    LAYER_TYPE comboSize = 0;
    bool checkWithMaxCombo = false;

    // Init
    EDGE_TYPE *trussness;                                           // 用于暂存分解结果
    auto fatherRes = new EDGE_TYPE[schemaNum];                      // 父亲combo的trussness min值
    auto fatherIds = new LAYER_TYPE[mlg->layersNum + 2];            // 每个combo最多也就layersNum个父亲
    auto triangles = new EDGE_TYPE[schemaNum * mlg->layersNum];     // 存每条边目前的三角形数目，防止多次计算
    mlg->genTriangle(triangles);

    c = new Combo(mlg->layersNum);
    c->lamb = mlg->layersNum;
    trussness = new EDGE_TYPE[schemaNum];
    for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
        fatherRes[edge] = mlg->maxNodesNum;
    }
    maxTrussness = FoTrussDecompositionForGivenCombo(mlg, trussness, triangles, fatherRes, c);
    if (maxTrussness != 0) {
        result[comboSize] = trussness;
        comboSize++;
        comboIndex->insert(c);
        comboIndex->rootComboId = 1;
        checkWithMaxCombo = true;
    } else {
        delete[] trussness;
    }
    delete c;

    auto possibleCombos = new queue<Combo *>();
    c = new Combo(mlg->layersNum);
    possibleCombos->push(c);

    while (not possibleCombos->empty()) {
        auto pCombo = possibleCombos->front();
        possibleCombos->pop();
        bool allFatherInIndex = comboIndex->allFatherInIndex(pCombo);
        if (not allFatherInIndex) {
            continue;
        }
        comboIndex->getOneHopFatherIds(pCombo, fatherIds, fatherSize);
        for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
            fatherRes[edge] = mlg->maxNodesNum;
        }

        trussness = new EDGE_TYPE[schemaNum];
        maxTrussness = FoTrussDecompositionForGivenCombo(mlg, trussness, triangles, fatherRes, pCombo);

        bool sameWithMaxCombo = true;
        if (checkWithMaxCombo) {
            for (EDGE_TYPE e = 0; e < mlg->schema->edgesNum; e++) {
                if (trussness[e] != result[0][e]) {
                    sameWithMaxCombo = false;
                    break;
                }
            }
        }

        if (maxTrussness <= MAX_PRUNE_TRUSSNESS or (checkWithMaxCombo and sameWithMaxCombo)) {
            delete[] trussness;
            delete pCombo;
        } else {
            result[comboSize] = trussness;
            comboSize++;
            comboIndex->insert(pCombo);
            for (LAYER_TYPE i = 0; i < fatherSize; i++) {
                LAYER_TYPE fatherId = fatherIds[i];
                bool sameWithFather = true;
                for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
                    if (result[fatherId][edge] != trussness[edge]) {
                        sameWithFather = false;
                        break;
                    }
                }
                if (sameWithFather) {
                    comboIndex->setTomb(fatherId);
                }
            }
            if (pCombo->focusNum < pCombo->lamb and pCombo->lamb < pCombo->layerNum) {
                LAYER_TYPE countStart = 0;
                if (pCombo->focusNum != 0) {
                    countStart = pCombo->focusLayers[pCombo->focusNum - 1] + 1;
                }
                for (LAYER_TYPE i = countStart; i < mlg->layersNum; i++) {
                    c = new Combo(pCombo);
                    c->focusNum++;
                    c->focusLayers[c->focusNum - 1] = i;
                    possibleCombos->push(c);
                }
            }
            if (pCombo->focusNum == 0 and pCombo->lamb <= mlg->layersNum) {
                c = new Combo(pCombo);
                c->lamb++;
                possibleCombos->push(c);
            }
            delete pCombo;
        }
    }
    timer->endTimer();
    LOG(INFO) << fmt::format("Decomposition cost (basic): {:.4f}s.", timer->getTimerSecond());

    auto result_size = double(comboSize * mlg->schema->edgesNum * sizeof(EDGE_TYPE)) / 1000.0 / 1000.0;

    LOG(INFO) << fmt::format("Non-empty meaningful Trussness size: {:.4f}MB.", result_size);
    LOG(INFO) << fmt::format("Non-empty meaningful combo num: {}.", comboSize);
    LOG(INFO) << fmt::format("Pruned sons = {}, Pruned fathers = {}.", comboIndex->prunedSonComboNum(mlg->layersNum),
                             comboIndex->prunedFatherComboNum());


    // De-allocate
    delete[] triangles;
    delete[] fatherRes;
    delete[] fatherIds;
    delete timer;
}


void FoTrussDecompositionPlus(MLGWithSchema *mlg, EDGE_TYPE **result, ComboIndex *comboIndex) {
    auto timer = new Timer();
    timer->startTimer();
    Combo *c;
    EDGE_TYPE maxTrussness = 0;                                     // 用于判断分解的truss是否为空
    EDGE_TYPE schemaNum = mlg->schema->edgesNum;                    // Edge-schema的总数
    LAYER_TYPE fatherSize;
    LAYER_TYPE comboSize = 0;
    bool checkWithMaxCombo = false;

    // Init
    EDGE_TYPE *trussness;                                           // 用于暂存分解结果
    auto fatherRes = new EDGE_TYPE[schemaNum];                      // 父亲combo的trussness min值
    auto fatherIds = new LAYER_TYPE[mlg->layersNum + 2];            // 每个combo最多也就layersNum个父亲
    auto triangles = new EDGE_TYPE[schemaNum * mlg->layersNum];     // 存每条边目前的三角形数目，防止多次计算
    mlg->genTriangle(triangles);

    c = new Combo(mlg->layersNum);
    c->lamb = mlg->layersNum;
    trussness = new EDGE_TYPE[schemaNum];
    for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
        fatherRes[edge] = mlg->maxNodesNum;
    }
    maxTrussness = FoTrussDecompositionForGivenCombo(mlg, trussness, triangles, fatherRes, c);
    if (maxTrussness != 0) {
        result[comboSize] = trussness;
        comboSize++;
        comboIndex->insert(c);
        comboIndex->rootComboId = 1;
        checkWithMaxCombo = true;
    } else {
        delete[] trussness;
    }
    delete c;

    auto possibleCombos = new queue<Combo *>();
    c = new Combo(mlg->layersNum);
    possibleCombos->push(c);

    while (not possibleCombos->empty()) {
        auto pCombo = possibleCombos->front();
        possibleCombos->pop();
        bool allFatherInIndex = comboIndex->allFatherInIndex(pCombo);
        if (not allFatherInIndex) {
            continue;
        }
        comboIndex->getOneHopFatherIds(pCombo, fatherIds, fatherSize);
        for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
            fatherRes[edge] = mlg->maxNodesNum;
        }
        if (fatherSize == 0) {
            trussness = new EDGE_TYPE[schemaNum];
            maxTrussness = FoTrussDecompositionForGivenCombo(mlg, trussness, triangles, fatherRes, pCombo);
        } else {
            EDGE_TYPE count = 0;
            for (LAYER_TYPE i = 0; i < fatherSize; i++) {
                LAYER_TYPE fatherId = fatherIds[i];
                for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
                    fatherRes[edge] = min(fatherRes[edge], result[fatherId][edge]);
                    if (fatherRes[edge] == 0) count++;
                }
            }
            trussness = new EDGE_TYPE[schemaNum];
            maxTrussness = FoTrussDecompositionForGivenComboUsingFather(mlg, trussness, fatherRes, pCombo);
        }

        bool sameWithMaxCombo = true;
        if (checkWithMaxCombo) {
            for (EDGE_TYPE e = 0; e < mlg->schema->edgesNum; e++) {
                if (trussness[e] != result[0][e]) {
                    sameWithMaxCombo = false;
                    break;
                }
            }
        }

        if (maxTrussness <= MAX_PRUNE_TRUSSNESS or (checkWithMaxCombo and sameWithMaxCombo)) {
            delete[] trussness;
            delete pCombo;
        } else {
            result[comboSize] = trussness;
            comboSize++;
            comboIndex->insert(pCombo);
            for (LAYER_TYPE i = 0; i < fatherSize; i++) {
                LAYER_TYPE fatherId = fatherIds[i];
                bool sameWithFather = true;
                for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
                    if (result[fatherId][edge] != trussness[edge]) {
                        sameWithFather = false;
                        break;
                    }
                }
                if (sameWithFather) {
                    comboIndex->setTomb(fatherId);
                }
            }
            if (pCombo->focusNum < pCombo->lamb and pCombo->lamb < pCombo->layerNum) {
                LAYER_TYPE countStart = 0;
                if (pCombo->focusNum != 0) {
                    countStart = pCombo->focusLayers[pCombo->focusNum - 1] + 1;
                }
                for (LAYER_TYPE i = countStart; i < mlg->layersNum; i++) {
                    c = new Combo(pCombo);
                    c->focusNum++;
                    c->focusLayers[c->focusNum - 1] = i;
                    possibleCombos->push(c);
                }
            }
            if (pCombo->focusNum == 0 and pCombo->lamb <= mlg->layersNum) {
                c = new Combo(pCombo);
                c->lamb++;
                possibleCombos->push(c);
            }
            delete pCombo;
        }
    }
    timer->endTimer();
    LOG(INFO) << fmt::format("Decomposition cost (Plus): {:.4f}s.", timer->getTimerSecond());

    auto result_size = double(comboSize * mlg->schema->edgesNum * sizeof(EDGE_TYPE)) / 1000.0 / 1000.0;

    LOG(INFO) << fmt::format("Non-empty meaningful Trussness size: {:.4f}MB.", result_size);
    LOG(INFO) << fmt::format("Non-empty meaningful combo num: {}.", comboSize);
    LOG(INFO) << fmt::format("Pruned sons = {}, Pruned fathers = {}.", comboIndex->prunedSonComboNum(mlg->layersNum),
                             comboIndex->prunedFatherComboNum());


    // De-allocate
    delete[] triangles;
    delete[] fatherRes;
    delete[] fatherIds;
    delete timer;
}


void FoTrussDecompositionAll(MLGWithSchema *mlg, EDGE_TYPE **result, ComboIndex *comboIndex) {
    auto timer = new Timer();
    timer->startTimer();
    Combo *c;
    EDGE_TYPE maxTrussness = 0;                                     // 用于判断分解的truss是否为空
    EDGE_TYPE schemaNum = mlg->schema->edgesNum;                    // Edge-schema的总数
    LAYER_TYPE fatherSize;
    LAYER_TYPE comboSize = 0;
    bool checkWithMaxCombo = false;

    // Init
    EDGE_TYPE *trussness;                                           // 用于暂存分解结果
    auto fatherRes = new EDGE_TYPE[schemaNum];                      // 父亲combo的trussness min值
    auto fatherIds = new LAYER_TYPE[mlg->layersNum + 2];            // 每个combo最多也就layersNum个父亲
    auto triangles = new EDGE_TYPE[schemaNum * mlg->layersNum];     // 存每条边目前的三角形数目，防止多次计算
    mlg->genTriangle(triangles);

    c = new Combo(mlg->layersNum);
    c->lamb = mlg->layersNum;
    trussness = new EDGE_TYPE[schemaNum];
    for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
        fatherRes[edge] = mlg->maxNodesNum;
    }
    maxTrussness = FoTrussDecompositionForGivenCombo(mlg, trussness, triangles, fatherRes, c);
    if (maxTrussness != 0) {
        result[comboSize] = trussness;
        comboSize++;
        comboIndex->insert(c);
        comboIndex->rootComboId = 1;
        checkWithMaxCombo = true;
    } else {
        delete[] trussness;
    }
    delete c;

    auto possibleCombos = new queue<Combo *>();
    c = new Combo(mlg->layersNum);
    possibleCombos->push(c);

    while (not possibleCombos->empty()) {
        auto pCombo = possibleCombos->front();
        possibleCombos->pop();
        bool allFatherInIndex = comboIndex->allFatherInIndex(pCombo);
        if (not allFatherInIndex) {
            continue;
        }
        comboIndex->getOneHopFatherIds(pCombo, fatherIds, fatherSize);
        for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
            fatherRes[edge] = mlg->maxNodesNum;
        }
        if (fatherSize == 0) {
            trussness = new EDGE_TYPE[schemaNum];
            maxTrussness = FoTrussDecompositionForGivenCombo(mlg, trussness, triangles, fatherRes, pCombo);
        } else {
            EDGE_TYPE count = 0;
            for (LAYER_TYPE i = 0; i < fatherSize; i++) {
                LAYER_TYPE fatherId = fatherIds[i];
                for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
                    fatherRes[edge] = min(fatherRes[edge], result[fatherId][edge]);
                    if (fatherRes[edge] == 0) count++;
                }
            }
            // 这里是可以继续优化的，using father需要重新计算三角形数目，但可以省掉一些peeling，或许可以计算一下两者的cost
            if (count > schemaNum / 2) {
                trussness = new EDGE_TYPE[schemaNum];
                maxTrussness = FoTrussDecompositionForGivenComboUsingFather(mlg, trussness, fatherRes, pCombo);
            } else {
                trussness = new EDGE_TYPE[schemaNum];
                maxTrussness = FoTrussDecompositionForGivenCombo(mlg, trussness, triangles, fatherRes, pCombo);
            }
        }

        bool sameWithMaxCombo = true;
        if (checkWithMaxCombo) {
            for (EDGE_TYPE e = 0; e < mlg->schema->edgesNum; e++) {
                if (trussness[e] != result[0][e]) {
                    sameWithMaxCombo = false;
                    break;
                }
            }
        }

        if (maxTrussness <= MAX_PRUNE_TRUSSNESS or (checkWithMaxCombo and sameWithMaxCombo)) {
            delete[] trussness;
            delete pCombo;
        } else {
            result[comboSize] = trussness;
            comboSize++;
            comboIndex->insert(pCombo);
            for (LAYER_TYPE i = 0; i < fatherSize; i++) {
                LAYER_TYPE fatherId = fatherIds[i];
                bool sameWithFather = true;
                for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
                    if (result[fatherId][edge] != trussness[edge]) {
                        sameWithFather = false;
                        break;
                    }
                }
                if (sameWithFather) {
                    comboIndex->setTomb(fatherId);
                }
            }
            if (pCombo->focusNum < pCombo->lamb and pCombo->lamb < pCombo->layerNum) {
                LAYER_TYPE countStart = 0;
                if (pCombo->focusNum != 0) {
                    countStart = pCombo->focusLayers[pCombo->focusNum - 1] + 1;
                }
                for (LAYER_TYPE i = countStart; i < mlg->layersNum; i++) {
                    c = new Combo(pCombo);
                    c->focusNum++;
                    c->focusLayers[c->focusNum - 1] = i;
                    possibleCombos->push(c);
                }
            }
            if (pCombo->focusNum == 0 and pCombo->lamb <= mlg->layersNum) {
                c = new Combo(pCombo);
                c->lamb++;
                possibleCombos->push(c);
            }
            delete pCombo;
        }
    }
    timer->endTimer();
    LOG(INFO) << fmt::format("Decomposition cost: {:.4f}s.", timer->getTimerSecond());

    auto result_size = double(comboSize * mlg->schema->edgesNum * sizeof(EDGE_TYPE)) / 1000.0 / 1000.0;

    LOG(INFO) << fmt::format("Non-empty meaningful Trussness size: {:.4f}MB.", result_size);
    LOG(INFO) << fmt::format("Calculated combo num: {}.", comboSize);
    LOG(INFO) << fmt::format("Pruned sons = {}, Pruned fathers = {}, Useful Combos = {}.",
                             comboIndex->prunedSonComboNum(mlg->layersNum),
                             comboIndex->prunedFatherComboNum(),
                             comboSize - comboIndex->prunedFatherComboNum());


    // De-allocate
    delete[] triangles;
    delete[] fatherRes;
    delete[] fatherIds;
    delete timer;
}


void edgesTriangleConnectedByGivenEdgeBasic(MLGWithSchema *mlg, EDGE_TYPE edge, NODE_TYPE *ret, NODE_TYPE &size) {
    // 最基础的连通
    NODE_TYPE u, v, w;
    size = 0;
    mlg->schema->idToEdge(edge, &u, &v);
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
        ret[size] = w;
        size++;
    }
    delete[] commonNeighbor;
}


void edgesTriangleConnectedByGivenEdgeShare(MLGWithSchema *mlg, EDGE_TYPE edge, NODE_TYPE *ret, NODE_TYPE &size) {
    NODE_TYPE u, v, w;
    size = 0;
    set<NODE_TYPE> commonNeighborSet;
    mlg->schema->idToEdge(edge, &u, &v);
    for (LAYER_TYPE layer = 0; layer < mlg->layersNum; layer++) {
        auto u_start = mlg->graph_layers[layer].start(u);
        auto u_end = mlg->graph_layers[layer].end(u);
        auto uList = &mlg->graph_layers[layer].edges[u_start];
        auto uLen = u_end - u_start;
        auto v_start = mlg->graph_layers[layer].start(v);
        auto v_end = mlg->graph_layers[layer].end(v);
        auto vList = &mlg->graph_layers[layer].edges[v_start];
        auto vLen = v_end - v_start;
        auto commonNeighbor = new NODE_TYPE[max(uLen, vLen)];
        NODE_TYPE commonNeighborLen = intersect(uList, vList, uLen, vLen, commonNeighbor);
        for (NODE_TYPE neighbor_i = 0; neighbor_i < commonNeighborLen; neighbor_i++) {
            w = commonNeighbor[neighbor_i];
            commonNeighborSet.insert(w);
        }
        delete[] commonNeighbor;
    }
    for (auto node: commonNeighborSet) {
        ret[size] = node;
        size++;
    }
}


void TriangleConnectedTruss(MLGWithSchema *mlg, EDGE_TYPE *trussness, EqualTree *eqTree, Schema *usingSchema) {
//    auto eqTree = new EqualTree();
    NODE_TYPE u, v, w;
    EDGE_TYPE edge, e1, e2, pvt, compactId, superNodeId, mt;
    EDGE_TYPE maxTrussness = 0;
    EDGE_TYPE schemaNum = usingSchema->edgesNum;
    auto edgeIdToCompactId = new EDGE_TYPE[schemaNum + 1];
    auto compactIdToEdgeId = new EDGE_TYPE[schemaNum + 1];
    for (EDGE_TYPE i = 0; i < schemaNum + 1; i++) {
        edgeIdToCompactId[i] = 0;
    }
    for (EDGE_TYPE i = 0; i < schemaNum + 1; i++) {
        compactIdToEdgeId[i] = 0;
    }
    auto parentPivots = new set<EDGE_TYPE>();
    auto neighbors = new NODE_TYPE[mlg->sumGraph->maxDegree];
    NODE_TYPE neighborSize = 0;

    // 先生成有序的边集
    for (edge = 0; edge < schemaNum; edge++) {
        maxTrussness = max(maxTrussness, trussness[edge]);
    }
    auto trussnessToNum = *new vector<EDGE_TYPE>(maxTrussness + 1, 0);
    for (edge = 0; edge < schemaNum; edge++) {
        trussnessToNum[trussness[edge]] += 1;
    }
    trussnessToNum[0] = 0;
    for (auto i = 2; i < trussnessToNum.size(); i++) {
        trussnessToNum[i] += trussnessToNum[i - 1];
    }
    for (auto i = trussnessToNum.size() - 1; i > 0; i--) {
        trussnessToNum[i] = trussnessToNum[i - 1] + 1;
    }
    for (edge = 0; edge < schemaNum; edge++) {
        if (trussness[edge] == 0) continue;
        edgeIdToCompactId[edge] = trussnessToNum[trussness[edge]];
        compactIdToEdgeId[trussnessToNum[trussness[edge]]] = edge;
        trussnessToNum[trussness[edge]] += 1;
    }
    trussnessToNum[0] = 1;

    auto uf = new UnionFind(trussnessToNum[maxTrussness] + 1);
    auto compactIdToSuperNodeId = new EDGE_TYPE[trussnessToNum[maxTrussness] + 1];
    for (EDGE_TYPE i = 0; i < trussnessToNum[maxTrussness] + 1; i++) {
        compactIdToSuperNodeId[i] = 0;
    }

    for (EDGE_TYPE k = maxTrussness; k > 0; k--) {
        EDGE_TYPE compactIdStart = trussnessToNum[k - 1];
        EDGE_TYPE compactIdEnd = trussnessToNum[k];
//        LOG(INFO) << fmt::format("generating equal tree for k={}.", k);
        // 寻找parent
//        LOG(INFO) << "find parent ...";
        for (compactId = compactIdStart; compactId < compactIdEnd; compactId++) {
            edge = compactIdToEdgeId[compactId];
//            if (trussness[edge] != k) LOG(ERROR) << fmt::format("trussness={} in bucket={}.", trussness[edge], k);
            if (DEFAULT_TRIANGLE_CONNECT_LEVEL == triangleConnectLevelEnum(basic)) {
                edgesTriangleConnectedByGivenEdgeBasic(mlg, edge, neighbors, neighborSize);
            } else if (DEFAULT_TRIANGLE_CONNECT_LEVEL == triangleConnectLevelEnum(shareEdge)) {
                edgesTriangleConnectedByGivenEdgeShare(mlg, edge, neighbors, neighborSize);
            } else {
                edgesTriangleConnectedByGivenEdgeBasic(mlg, edge, neighbors, neighborSize);
            }
            usingSchema->idToEdge(edge, &u, &v);
            for (NODE_TYPE j = 0; j < neighborSize; j++) {
                w = neighbors[j];
                e1 = usingSchema->edgeToId(u, w);
                e2 = usingSchema->edgeToId(v, w);
                mt = min(trussness[e1], trussness[e2]);
                mt = min(mt, trussness[edge]);
                if (trussness[edge] == mt) {
                    pvt = uf->getPivot(edgeIdToCompactId[e1]);
                    parentPivots->insert(pvt);
                    pvt = uf->getPivot(edgeIdToCompactId[e2]);
                    parentPivots->insert(pvt);
                }
            }
        }

        // 连通
//        LOG(INFO) << "union find ...";
        for (compactId = compactIdStart; compactId < compactIdEnd; compactId++) {
            edge = compactIdToEdgeId[compactId];
            if (DEFAULT_TRIANGLE_CONNECT_LEVEL == triangleConnectLevelEnum(basic)) {
                edgesTriangleConnectedByGivenEdgeBasic(mlg, edge, neighbors, neighborSize);
            } else if (DEFAULT_TRIANGLE_CONNECT_LEVEL == triangleConnectLevelEnum(shareEdge)) {
                edgesTriangleConnectedByGivenEdgeShare(mlg, edge, neighbors, neighborSize);
            } else {
                edgesTriangleConnectedByGivenEdgeBasic(mlg, edge, neighbors, neighborSize);
            }
            usingSchema->idToEdge(edge, &u, &v);
            for (NODE_TYPE j = 0; j < neighborSize; j++) {
                w = neighbors[j];
                e1 = usingSchema->edgeToId(u, w);
                e2 = usingSchema->edgeToId(v, w);
                if (trussness[e1] >= k and trussness[e2] >= k) {
                    uf->toUnion(compactId, edgeIdToCompactId[e1]);
                    uf->toUnion(compactId, edgeIdToCompactId[e2]);
                }
            }
        }

        // 建树
//        LOG(INFO) << "construct tree...";
        for (compactId = compactIdStart; compactId < compactIdEnd; compactId++) {
            pvt = uf->getPivot(compactId);
            if (compactIdToSuperNodeId[pvt] == 0) {
                superNodeId = eqTree->superNodeNum;
                eqTree->addSuperNode(k);
                compactIdToSuperNodeId[pvt] = superNodeId;
            }
            superNodeId = compactIdToSuperNodeId[pvt];
            compactIdToSuperNodeId[compactId] = superNodeId;
            EDGE_TYPE edgeId = compactIdToEdgeId[compactId];
            eqTree->addEdgeToSuperNode(superNodeId, edgeId);
        }

        // 树的父子关系
//        LOG(INFO) << "setup son-parent ...";
        for (auto e: *parentPivots) {
            pvt = uf->getPivot(e);
            auto child = compactIdToSuperNodeId[e];
            auto parent = compactIdToSuperNodeId[pvt];
            if (child != parent) {
                eqTree->addParent(child, parent);
                eqTree->addChild(parent, child);
            }
        }

        parentPivots->clear();
    }

    delete uf;
    delete[] compactIdToEdgeId;
    delete[] edgeIdToCompactId;
    delete[] compactIdToSuperNodeId;
    delete[] neighbors;
    eqTree->compactTree();
    trussnessToNum.clear();
    parentPivots->clear();
}


void queryNodeCommunity(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode, EDGE_TYPE trussnessBound,
                        vector<vector<NODE_TYPE>> &result) {
    result.clear();
    if (trussnessBound == 0) return;
    NODE_TYPE u, v;
    auto queryNodeNeighborStart = mlg->sumGraph->start(queryNode);
    auto queryNodeNeighborEnd = mlg->sumGraph->end(queryNode);
    vector<EDGE_TYPE> edgesAttachedToQueryNode;
    for (EDGE_TYPE neighbor_i = queryNodeNeighborStart; neighbor_i < queryNodeNeighborEnd; neighbor_i++) {
        auto neighbor = mlg->sumGraph->edges[neighbor_i];
        auto edge = mlg->schema->edgeToId(queryNode, neighbor);
        if (eqTree->getEdgeTrussness(edge) >= trussnessBound) {
            edgesAttachedToQueryNode.emplace_back(edge);
        }
    }
    auto superNodes = eqTree->edgesToSuperNodes(edgesAttachedToQueryNode);
    result.resize(superNodes.size());
    EDGE_TYPE resultId = 0;
    for (auto superNode: superNodes) {
        auto nodeInCommunity = new set<NODE_TYPE>();
        nodeInCommunity->insert(queryNode);
        auto temp = eqTree->getSubgraphBySuperNodeID(superNode);
        for (auto edge: temp) {
            mlg->schema->idToEdge(edge, &u, &v);
            nodeInCommunity->insert(u);
            nodeInCommunity->insert(v);
        }
        result[resultId].assign(nodeInCommunity->begin(), nodeInCommunity->end());
        resultId++;
    }
    superNodes.clear();
    superNodes.shrink_to_fit();
}


EDGE_TYPE getNodeMaxTrussness(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode) {
    NODE_TYPE u, v;
    EDGE_TYPE trussness = 0;
    auto queryNodeNeighborStart = mlg->sumGraph->start(queryNode);
    auto queryNodeNeighborEnd = mlg->sumGraph->end(queryNode);
    vector<EDGE_TYPE> edgesAttachedToQueryNode;
    for (EDGE_TYPE neighbor_i = queryNodeNeighborStart; neighbor_i < queryNodeNeighborEnd; neighbor_i++) {
        auto neighbor = mlg->sumGraph->edges[neighbor_i];
        auto edge = mlg->schema->edgeToId(queryNode, neighbor);
        edgesAttachedToQueryNode.emplace_back(edge);
    }
    for (auto edge: edgesAttachedToQueryNode) {
        auto t = eqTree->getEdgeTrussness(edge);
        trussness = max(trussness, t);
    }
    return trussness;
}


EDGE_TYPE getMaxTrussness(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode) {
    NODE_TYPE u, v;
    EDGE_TYPE maxTrussness = 0;
    auto queryNodeNeighborStart = mlg->sumGraph->start(queryNode);
    auto queryNodeNeighborEnd = mlg->sumGraph->end(queryNode);
    auto edgesAttachedToQueryNode = new vector<EDGE_TYPE>();
    auto result = new vector<NODE_TYPE>();
    auto nodeInCommunity = new set<NODE_TYPE>();
    for (EDGE_TYPE neighbor_i = queryNodeNeighborStart; neighbor_i < queryNodeNeighborEnd; neighbor_i++) {
        auto neighbor = mlg->sumGraph->edges[neighbor_i];
        auto edge = mlg->schema->edgeToId(queryNode, neighbor);
        maxTrussness = max(eqTree->getEdgeTrussness(edge), maxTrussness);
    }
    return maxTrussness;
}


void FoTrussIndexAll(MLGWithSchema *mlg, EqualTree **index, ComboIndex *ci, EDGE_TYPE **trussness) {
    LOG(INFO) << fmt::format("Build FoTruss Index For all non-empty FoTruss...");
    double totalIndexSize = 0;
    auto timer = new Timer();
    timer->startTimer();
    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        if (ci->tomb[i]) {
            index[i] = nullptr;
            continue;
        }
        index[i] = new EqualTree(mlg->schema->edgesNum);
        TriangleConnectedTruss(mlg, trussness[i], index[i], mlg->schema);
        delete[] trussness[i];
        totalIndexSize += index[i]->getMemUsage();
    }
    timer->endTimer();

    LOG(INFO) << fmt::format("Index Construction cost: {:.4f}s.", timer->getTimerSecond());
    LOG(INFO) << fmt::format("Index size: {:.4f}MB.", totalIndexSize);
    delete[] trussness;
    delete timer;
}


void FoTrussIndex(MLGWithSchema *mlg, const filesystem::path &outPath) {
    auto maxComboSize = LAYER_TYPE(pow(2, mlg->layersNum)) + mlg->layersNum * LAYER_TYPE(pow(2, mlg->layersNum - 1));
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *ci = new ComboIndex();
    FoTrussDecompositionAll(mlg, trussness, ci);
    auto **index = new EqualTree *[ci->comboNum + 1];
    auto **schemas = new Schema *[ci->comboNum + 1];
    FoTrussIndexAll(mlg, index, ci, trussness);
//    FoTrussIndexCompact(mlg, index, ci, trussness, schemas);
    double totalIndexSize = 0;
    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        if (ci->tomb[i]) {
            continue;
        }
        totalIndexSize += index[i]->getMemUsage();
        delete index[i];
    }
    delete ci;
}


void FoTrussIndexGreedy(MLGWithSchema *mlg, LAYER_TYPE maxComboSize, const vector<NODE_TYPE> &queries,
                        const vector<NODE_TYPE> &testQueries, vector<string> &combos, vector<EqualTree *> &index) {
    Combo *c;
    EDGE_TYPE maxTrussness;
    EDGE_TYPE schemaNum = mlg->schema->edgesNum;                    // Edge-schema的总数
    EDGE_TYPE *trussness;                                           // 用于暂存分解结果
    trussness = new EDGE_TYPE[schemaNum];
    auto fatherRes = new EDGE_TYPE[schemaNum];                      // 父亲combo的trussness min值
    auto triangles = new EDGE_TYPE[schemaNum * mlg->layersNum];     // 存每条边目前的三角形数目，防止多次计算
    mlg->genTriangle(triangles);
    vector<double> metrics;
    LAYER_TYPE minMetricId = 0;
    double minMetric = 0.0;
    LAYER_TYPE exploredComboNum = 0;

    auto possibleCombos = new queue<Combo *>();
    c = new Combo(mlg->layersNum);
    possibleCombos->push(c);

    while (not possibleCombos->empty()) {
        auto pCombo = possibleCombos->front();
        possibleCombos->pop();
        // 先得到combo的结果，并计算评价
        for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
            fatherRes[edge] = mlg->maxNodesNum;
        }
        maxTrussness = FoTrussDecompositionForGivenCombo(mlg, trussness, triangles, fatherRes, pCombo);
        exploredComboNum++;
        if (maxTrussness == 0) {
            LOG(INFO) << fmt::format("Combo {} Excluded because of empty core", pCombo->toString());
            delete pCombo;
            continue;
        }
        auto t = new EqualTree(mlg->schema->edgesNum);
        TriangleConnectedTruss(mlg, trussness, t, mlg->schema);
        vector<vector<NODE_TYPE>> qus;
        double sumMetric = 0.0;
        for (auto query: queries) {
            queryNodeCommunity(mlg, t, query, getMaxTrussness(mlg, t, query), qus);
            double bestDen = 0;
            for (const auto &qu: qus) {
                double den = density(qu, mlg, DEN_BETA);
                bestDen = max(bestDen, den);
            }
            sumMetric += bestDen;
        }
        // 然后对比combo和已有的combo，看是否可以代替。
        if (metrics.size() >= maxComboSize and sumMetric <= metrics[minMetricId]) {
            LOG(INFO) << fmt::format("Combo {} excluded because small metric {:.4f}", pCombo->toString(), sumMetric);
            delete pCombo;
            delete t;
            continue;
        }
        LOG(INFO) << fmt::format("Combo {} into heap with metric {:.4f}", pCombo->toString(), sumMetric);

        // 把这个代替的combo存下来
        if (metrics.size() < maxComboSize) {
            metrics.emplace_back(sumMetric);
            combos.emplace_back(pCombo->toString());
            index.emplace_back(t);
            minMetric = sumMetric;
            minMetricId = metrics.size() - 1;
            for (LAYER_TYPE l = 0; l < minMetricId; l++) {
                if (metrics[l] < minMetric) {
                    minMetricId = l;
                    minMetric = metrics[l];
                }
            }
        } else {
            metrics[minMetricId] = sumMetric;
            LOG(INFO) << fmt::format("Combo {} replace combo {}", pCombo->toString(), combos[minMetricId]);
            combos[minMetricId] = pCombo->toString();
            delete index[minMetricId];
            index[minMetricId] = t;
            minMetric = sumMetric;
            minMetricId = metrics.size() - 1;
            for (LAYER_TYPE l = 0; l < maxComboSize; l++) {
                if (metrics[l] < minMetric) {
                    minMetricId = l;
                    minMetric = metrics[l];
                }
            }
        }

        // 如果可以代替的话，把这个combo往孩子方向推一步加入到队列中
        if (pCombo->focusNum < pCombo->lamb and pCombo->lamb < pCombo->layerNum) {
            LAYER_TYPE countStart = 0;
            if (pCombo->focusNum != 0) {
                countStart = pCombo->focusLayers[pCombo->focusNum - 1] + 1;
            }
            for (LAYER_TYPE i = countStart; i < mlg->layersNum; i++) {
                c = new Combo(pCombo);
                c->focusNum++;
                c->focusLayers[c->focusNum - 1] = i;
                possibleCombos->push(c);
            }
        }
        if (pCombo->focusNum == 0 and pCombo->lamb <= mlg->layersNum) {
            c = new Combo(pCombo);
            c->lamb++;
            possibleCombos->push(c);
        }
    }

    double sumMetric = 0.0;
    for (auto query: testQueries) {
        double bestDen = 0;
        for (auto t: index) {
            vector<vector<NODE_TYPE>> qus;
            queryNodeCommunity(mlg, t, query, getMaxTrussness(mlg, t, query), qus);
            for (const auto &qu: qus) {
                double den = density(qu, mlg, DEN_BETA);
                bestDen = max(bestDen, den);
            }
        }
        sumMetric += bestDen;
    }
    double avgMetric = sumMetric / double(testQueries.size());

    LOG(INFO) << fmt::format("Explored Combo Num = {}", exploredComboNum);
    LOG(INFO) << fmt::format("Metric = {:.4f}", avgMetric);
    LOG(INFO) << fmt::format("Stored Below Combos:");
    for (LAYER_TYPE l = 0; l < maxComboSize; l++) {
        metrics[l] = metrics[l] / double(queries.size());
        LOG(INFO) << fmt::format("Combo Name = {}, Metric = {:.4f}", combos[l], metrics[l]);
    }

    delete[] triangles;
}


void FoTrussIndexGreedyFirstK(MLGWithSchema *mlg, LAYER_TYPE maxComboSize, const vector<NODE_TYPE> &testQueries) {
    Combo *c;
    EDGE_TYPE maxTrussness;
    EDGE_TYPE schemaNum = mlg->schema->edgesNum;                    // Edge-schema的总数
    EDGE_TYPE *trussness;                                           // 用于暂存分解结果
    trussness = new EDGE_TYPE[schemaNum];
    auto fatherRes = new EDGE_TYPE[schemaNum];                      // 父亲combo的trussness min值
    auto triangles = new EDGE_TYPE[schemaNum * mlg->layersNum];     // 存每条边目前的三角形数目，防止多次计算
    mlg->genTriangle(triangles);
    LAYER_TYPE exploredComboNum = 0;

    vector<string> combos;
    vector<EqualTree *> index;

    auto possibleCombos = new queue<Combo *>();
    c = new Combo(mlg->layersNum);
    possibleCombos->push(c);

    while (not possibleCombos->empty()) {
        auto pCombo = possibleCombos->front();
        possibleCombos->pop();
        // 先得到combo的结果，并计算评价
        for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
            fatherRes[edge] = mlg->maxNodesNum;
        }
        maxTrussness = FoTrussDecompositionForGivenCombo(mlg, trussness, triangles, fatherRes, pCombo);
        exploredComboNum++;
        if (maxTrussness == 0) {
            LOG(INFO) << fmt::format("Combo {} Excluded because of empty core", pCombo->toString());
            delete pCombo;
            continue;
        }
        auto t = new EqualTree(mlg->schema->edgesNum);
        TriangleConnectedTruss(mlg, trussness, t, mlg->schema);

        // 把这个代替的combo存下来
        if (index.size() < maxComboSize) {
            combos.emplace_back(pCombo->toString());
            index.emplace_back(t);
        } else {
            break;
        }

        // 如果可以代替的话，把这个combo往孩子方向推一步加入到队列中
        if (pCombo->focusNum < pCombo->lamb and pCombo->lamb < pCombo->layerNum) {
            LAYER_TYPE countStart = 0;
            if (pCombo->focusNum != 0) {
                countStart = pCombo->focusLayers[pCombo->focusNum - 1] + 1;
            }
            for (LAYER_TYPE i = countStart; i < mlg->layersNum; i++) {
                c = new Combo(pCombo);
                c->focusNum++;
                c->focusLayers[c->focusNum - 1] = i;
                possibleCombos->push(c);
            }
        }
        if (pCombo->focusNum == 0 and pCombo->lamb <= mlg->layersNum) {
            c = new Combo(pCombo);
            c->lamb++;
            possibleCombos->push(c);
        }
    }

    double sumMetric = 0.0;
    for (auto query: testQueries) {
        double bestDen = 0;
        for (auto t: index) {
            vector<vector<NODE_TYPE>> qus;
            queryNodeCommunity(mlg, t, query, getMaxTrussness(mlg, t, query), qus);
            for (const auto &qu: qus) {
                double den = density(qu, mlg, DEN_BETA);
                bestDen = max(bestDen, den);
            }
        }
        sumMetric += bestDen;
    }
    double avgMetric = sumMetric / double(testQueries.size());

    LOG(INFO) << fmt::format("Explored Combo Num = {}", exploredComboNum);
    LOG(INFO) << fmt::format("Metric = {:.4f}", avgMetric);
    LOG(INFO) << fmt::format("Stored Below Combos:");
    for (LAYER_TYPE l = 0; l < maxComboSize; l++) {
        LOG(INFO) << fmt::format("Combo Name = {}", combos[l]);
    }

    delete[] triangles;
}


void FoTrussIndexGreedyFirm(MLGWithSchema *mlg, LAYER_TYPE maxComboSize, const vector<NODE_TYPE> &testQueries) {
    Combo *c;
    EDGE_TYPE maxTrussness;
    EDGE_TYPE schemaNum = mlg->schema->edgesNum;                    // Edge-schema的总数
    EDGE_TYPE *trussness;                                           // 用于暂存分解结果
    trussness = new EDGE_TYPE[schemaNum];
    auto fatherRes = new EDGE_TYPE[schemaNum];                      // 父亲combo的trussness min值
    auto triangles = new EDGE_TYPE[schemaNum * mlg->layersNum];     // 存每条边目前的三角形数目，防止多次计算
    mlg->genTriangle(triangles);
    LAYER_TYPE exploredComboNum = 0;

    vector<string> combos;
    vector<EqualTree *> index;

    auto possibleCombos = new queue<Combo *>();
    c = new Combo(mlg->layersNum);
    possibleCombos->push(c);

    while (not possibleCombos->empty()) {
        auto pCombo = possibleCombos->front();
        possibleCombos->pop();
        // 先得到combo的结果，并计算评价
        for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
            fatherRes[edge] = mlg->maxNodesNum;
        }
        maxTrussness = FoTrussDecompositionForGivenCombo(mlg, trussness, triangles, fatherRes, pCombo);
        exploredComboNum++;
        if (maxTrussness == 0) {
            LOG(INFO) << fmt::format("Combo {} Excluded because of empty core", pCombo->toString());
            delete pCombo;
            continue;
        }
        auto t = new EqualTree(mlg->schema->edgesNum);
        TriangleConnectedTruss(mlg, trussness, t, mlg->schema);

        // 把这个代替的combo存下来
        if (index.size() < maxComboSize) {
            combos.emplace_back(pCombo->toString());
            index.emplace_back(t);
        } else {
            break;
        }

        // 如果可以代替的话，把这个combo往孩子方向推一步加入到队列中
        if (pCombo->focusNum == 0 and pCombo->lamb <= mlg->layersNum) {
            c = new Combo(pCombo);
            c->lamb++;
            possibleCombos->push(c);
        }
    }

    double sumMetric = 0.0;
    for (auto query: testQueries) {
        double bestDen = 0;
        for (auto t: index) {
            vector<vector<NODE_TYPE>> qus;
            queryNodeCommunity(mlg, t, query, getMaxTrussness(mlg, t, query), qus);
            for (const auto &qu: qus) {
                double den = density(qu, mlg, DEN_BETA);
                bestDen = max(bestDen, den);
            }
        }
        sumMetric += bestDen;
    }
    double avgMetric = sumMetric / double(testQueries.size());

    LOG(INFO) << fmt::format("Explored Combo Num = {}", exploredComboNum);
    LOG(INFO) << fmt::format("Metric = {:.4f}", avgMetric);
    LOG(INFO) << fmt::format("Stored Below Combos:");
    for (LAYER_TYPE l = 0; l < maxComboSize; l++) {
        LOG(INFO) << fmt::format("Combo Name = {}", combos[l]);
    }

    delete[] triangles;
}


void saveTrussness(MLGWithSchema *mlg, EDGE_TYPE **trussness, const filesystem::path &outPath) {
    ofstream outfile;
    LOG(INFO) << fmt::format("Storing file into {}", outPath.string());
    outfile.open(outPath, ios::out | ios::trunc);
    for (LAYER_TYPE lamb = 1; lamb < mlg->layersNum + 1; lamb++) {
        outfile << fmt::format("lambda={}", lamb) << endl;
        outfile << "[";
        for (EDGE_TYPE j = 0; j < mlg->schema->edgesNum; j++) {
            outfile << trussness[lamb][j] << ", ";
        }
        outfile << "]" << endl;
    }
    outfile.close();
}
