//
// Created by 王润安 on 2023/12/5.
//

#include "Algo.h"
#include "algorithm"
#include "iterator"
#include <cmath>
#include "queue"
#include "../Utils/logging.h"
#include "vector"

inline void swap(NODE_TYPE a[], LAYER_TYPE i, LAYER_TYPE j) {
    NODE_TYPE temp = a[i];
    a[i] = a[j];
    a[j] = temp;
}

inline void refreshHeapAfterInsert(NODE_TYPE *heap, LAYER_TYPE heapSize) {
    LAYER_TYPE father = 0;
    LAYER_TYPE selectedChild, leftChild, rightChild;
    while (father < heapSize) {
        leftChild = 2 * father + 1;
        rightChild = 2 * father + 2;
        if (leftChild >= heapSize) {
            // father就是最下一层的节点，不用再换了
            break;
        } else {
            if (rightChild >= heapSize) {
                // father只有左节点，不用换就直接结束，否则换完结束
                if (heap[father] > heap[leftChild]) {
                    swap(heap, father, leftChild);
                }
                break;
            } else {
                // father左右节点都有，在左右节点中选较小的那个换上来
                selectedChild = heap[leftChild] > heap[rightChild] ? rightChild : leftChild;
                if (heap[father] > heap[selectedChild]) {
                    swap(heap, father, selectedChild);
                    father = selectedChild;
                } else {
                    // father更小，说明father处于合适的位置，结束
                    break;
                }
            }
        }
    }
}

inline NODE_TYPE findByHeap(const NODE_TYPE *degreeList, LAYER_TYPE layerNum, LAYER_TYPE k) {
    //用前k个元素建立最小堆
    NODE_TYPE i = 0;
    NODE_TYPE tempLargestHeap[k];
    for (i = 0; i < k; i++) {
        tempLargestHeap[i] = 0;
    }
    for (i = 0; i < layerNum; i++) {
        if (degreeList[i] > tempLargestHeap[0]) {
            tempLargestHeap[0] = degreeList[i];
            refreshHeapAfterInsert(tempLargestHeap, k);
        }
    }
    return tempLargestHeap[0];
}

NODE_TYPE findTopK(const NODE_TYPE *degreeList, LAYER_TYPE layerNum, LAYER_TYPE k) {
    return findByHeap(degreeList, layerNum, k);
}

NODE_TYPE intersect(const NODE_TYPE *list1, const NODE_TYPE *list2, NODE_TYPE len1, NODE_TYPE len2, NODE_TYPE *ans) {
    NODE_TYPE i = 0;
    NODE_TYPE j = 0;
    NODE_TYPE ansLen = 0;
    while (i < len1 and j < len2) {
        if (list1[i] == list2[j]) {
            ans[ansLen] = list1[i];
            ansLen++;
            i++;
            j++;
        } else if (list1[i] < list2[j]) {
            i++;
        } else {
            // list1[i] > list2[j]
            j++;
        }
//        if (list1[len1 - 1] < list2[j] or list2[len2 - 1] < list1[i]) break;
    }
    return ansLen;
}

double f1Score(vector<NODE_TYPE> community, vector<NODE_TYPE> gt_community) {
    std::set<NODE_TYPE> c(community.begin(), community.end());
    std::set<NODE_TYPE> gtc(gt_community.begin(), gt_community.end());
    auto inter = new std::vector<NODE_TYPE>();
    std::set_intersection(c.begin(), c.end(), gtc.begin(), gtc.end(),
                          std::back_inserter<std::vector<NODE_TYPE>>(*inter));
    double pre = double(inter->size()) / double(gt_community.size());
    double recall = double(inter->size()) / double(community.size());
    double f1 = 2 * pre * recall / (pre + recall);
    return f1;
}

double density(const vector<NODE_TYPE> &subgraph, MLGWithSchema *mlg, double beta) {
    double den = 0.0;
    vector<EDGE_TYPE> edgeSizeEachLayer(mlg->layersNum, 0);
    set<NODE_TYPE> graph(subgraph.begin(), subgraph.end());
    for (LAYER_TYPE layer = 0; layer < mlg->layersNum; layer++) {
        for (auto v: graph) {
            for (NODE_TYPE uid = mlg->graph_layers[layer].posLargerThanV(v);
                 uid < mlg->graph_layers[layer].end(v); uid++) {
                NODE_TYPE u = mlg->graph_layers[layer].edges[uid];
                if (graph.find(u) != graph.end()) {
                    edgeSizeEachLayer[layer]++;
                }
            }
        }
    }
    std::sort(edgeSizeEachLayer.begin(), edgeSizeEachLayer.end(), greater<>());
    for (LAYER_TYPE layer = 0; layer < mlg->layersNum; layer++) {
        double denThisLayer = double(edgeSizeEachLayer[layer]) / double(subgraph.size());
        double thisDen = denThisLayer * pow((layer + 1), beta);
        if (thisDen > den) den = thisDen;
    }
    return den;
}

double weightedDensity(const vector<NODE_TYPE> &subgraph, MLGWithSchema *mlg, const vector<double> &weight) {
    double den = 0.0;
    vector<EDGE_TYPE> edgeSizeEachLayer(mlg->layersNum, 0);
    set<NODE_TYPE> graph(subgraph.begin(), subgraph.end());
    for (LAYER_TYPE layer = 0; layer < mlg->layersNum; layer++) {
        for (auto v: graph) {
            for (NODE_TYPE uid = mlg->graph_layers[layer].posLargerThanV(v);
                 uid < mlg->graph_layers[layer].end(v); uid++) {
                NODE_TYPE u = mlg->graph_layers[layer].edges[uid];
                if (graph.find(u) != graph.end()) {
                    edgeSizeEachLayer[layer]++;
                }
            }
        }
    }
    for (LAYER_TYPE layer = 0; layer < mlg->layersNum; layer++) {
        double denThisLayer = double(edgeSizeEachLayer[layer]) / double(subgraph.size());
        den += denThisLayer * weight[layer];
    }
    return den;
}


void densityOfEachLayer(const vector<NODE_TYPE> &subgraph, MLGWithSchema *mlg, vector<double> &den) {
    den.clear();
    vector<EDGE_TYPE> edgeSizeEachLayer(mlg->layersNum, 0);
    set<NODE_TYPE> graph(subgraph.begin(), subgraph.end());
    for (LAYER_TYPE layer = 0; layer < mlg->layersNum; layer++) {
        for (auto v: graph) {
            for (NODE_TYPE uid = mlg->graph_layers[layer].posLargerThanV(v);
                 uid < mlg->graph_layers[layer].end(v); uid++) {
                NODE_TYPE u = mlg->graph_layers[layer].edges[uid];
                if (graph.find(u) != graph.end()) {
                    edgeSizeEachLayer[layer]++;
                }
            }
        }
    }
    for (LAYER_TYPE layer = 0; layer < mlg->layersNum; layer++) {
        double denThisLayer;
        if (subgraph.size() == 0) denThisLayer = 0;
        else denThisLayer = double(edgeSizeEachLayer[layer]) / double(subgraph.size());
        den.emplace_back(denThisLayer);
    }
}

NODE_TYPE diameter(const vector<NODE_TYPE> &subgraph, MLGWithSchema *mlg) {
    NODE_TYPE diam = 0;
    // 进行BFS，需要下面三个数据结构，前两个以subGID，后一个存mlgID
    vector<bool> visited;
    vector<NODE_TYPE> distance;
    queue<NODE_TYPE> toVisit;
    vector<bool> visitedInit;
    visitedInit.assign(mlg->sumGraph->nodesNum, true);
    for (auto srcNode: subgraph) {
        visitedInit[srcNode] = false;
    }

    for (auto srcNode: subgraph) {
        // srcNode是出发点，子图中每个点都要算一遍，初始化
        visited.assign(visitedInit.begin(), visitedInit.end());
        distance.assign(mlg->sumGraph->nodesNum, 0);
        visited[srcNode] = true;
        toVisit.push(srcNode);
        NODE_TYPE visitNum = 1;
        NODE_TYPE thisMaxDiam = 0;
        while (not toVisit.empty()) {
            if (subgraph.size() - visitNum + thisMaxDiam <= diam) {
                break;
            }
            NODE_TYPE visitNode = toVisit.front();
            toVisit.pop();
            for (NODE_TYPE neighborId = mlg->sumGraph->start(visitNode);
                 neighborId < mlg->sumGraph->end(visitNode); neighborId++) {
                // neighbor是mlgID
                NODE_TYPE neighbor = mlg->sumGraph->edges[neighborId];
                // subgraphNeighbor是subGID
                if (not visited[neighbor]) {
                    distance[neighbor] = distance[visitNode] + 1;
                    thisMaxDiam = max(thisMaxDiam, distance[neighbor]);
                    toVisit.push(neighbor);
                    visited[neighbor] = true;
                    visitNum++;
                }
            }
        }
        while (not toVisit.empty()) {
            toVisit.pop();
        }

        diam = max(thisMaxDiam, diam);
    }
    return diam + 1; // 包括src自己
}


struct visitPair {
    NODE_TYPE node;
    LAYER_TYPE layer;
    NODE_TYPE distance;
};


NODE_TYPE diameterMultilayer(const vector<NODE_TYPE> &subgraph, MLGWithSchema *mlg) {
    NODE_TYPE diam = 0;
    // 进行BFS，需要下面三个数据结构，前两个以subGID，后一个存mlgID
    vector<bool> visited;
    vector<NODE_TYPE> distance;
    queue<visitPair> toVisit;
    vector<bool> visitedInit;
    visitedInit.assign(mlg->sumGraph->nodesNum, true);
    for (auto srcNode: subgraph) {
        visitedInit[srcNode] = false;
    }

    for (auto srcNode: subgraph) {
        // srcNode是出发点，子图中每个点都要算一遍，初始化
        visited.assign(visitedInit.begin(), visitedInit.end());
        distance.assign(mlg->sumGraph->nodesNum, mlg->sumGraph->nodesNum + 10);
        visited[srcNode] = true;
        for (LAYER_TYPE l = 0; l < mlg->layersNum; l++) {
            visitPair vp{srcNode, l, 0};
            toVisit.push(vp);
        }

        NODE_TYPE thisMaxDiam = 0;
        while (not toVisit.empty()) {
            visitPair vp = toVisit.front();
            NODE_TYPE visitNode = vp.node;
            LAYER_TYPE visitLayer = vp.layer;
            NODE_TYPE visitDistance = vp.distance;
            toVisit.pop();
            for (NODE_TYPE neighborId = mlg->graph_layers[visitLayer].start(visitNode);
                 neighborId < mlg->graph_layers[visitLayer].end(visitNode); neighborId++) {
                // neighbor是mlgID
                NODE_TYPE neighbor = mlg->graph_layers[visitLayer].edges[neighborId];
                // subgraphNeighbor是subGID
                distance[neighbor] = min(visitDistance + 1, distance[neighbor]);
                thisMaxDiam = max(thisMaxDiam, distance[neighbor]);
                if (not visited[neighbor]) {
                    visitPair nextVp{neighbor, visitLayer, distance[neighbor]};
                    toVisit.push(nextVp);
                    for (LAYER_TYPE l = 0; l < mlg->layersNum; l++) {
                        if (l != visitLayer) {
                            visitPair nextNextVp{neighbor, l, distance[neighbor] + 1};
                            toVisit.push(nextNextVp);
                        }
                    }
                    visited[neighbor] = true;
                }
            }

        }
        diam = max(thisMaxDiam, diam);
    }
    return diam + 1; // 包括src自己
}
