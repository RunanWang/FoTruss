//
// Created by 王润安 on 2023/12/28.
//

#ifndef MLG_CS_EQUALTREE_H
#define MLG_CS_EQUALTREE_H

#include "../constant.h"
#include "vector"
#include "unordered_map"
#include "Combo.h"
#include "filesystem"

using namespace std;

class EqualTree {
public:
    vector<vector<EDGE_TYPE>> edgesInSuperNodes;
    vector<EDGE_TYPE> superNodeParent;
    vector<vector<EDGE_TYPE>> superNodeChild;
    EDGE_TYPE superNodeNum = 1;
    // superNodeID从1开始计数，0表示目前这个点没被加入树中

    vector<EDGE_TYPE> edgeIdToSuperNodeId;
    vector<EDGE_TYPE> superNodeIdToTrussness;

    explicit EqualTree(EDGE_TYPE edgeSize);

    ~EqualTree() = default;

    void addSuperNode(EDGE_TYPE trussness);

    void addEdgeToSuperNode(EDGE_TYPE superNodeId, EDGE_TYPE toAddEdgeId);

    void addChild(EDGE_TYPE thisSuperNodeId, EDGE_TYPE childSuperNodeId);

    void addParent(EDGE_TYPE thisSuperNodeId, EDGE_TYPE parentSuperNodeId);

    void compactTree();

    EDGE_TYPE getSuperNodeID(EDGE_TYPE edgeId);

    EDGE_TYPE getEdgeTrussness(EDGE_TYPE edgeId);

    vector<EDGE_TYPE> edgesToSuperNodes(const vector<EDGE_TYPE> &edges);

    vector<EDGE_TYPE> getSubgraphBySuperNodeID(EDGE_TYPE superNodeID);

    double getMemUsage();
};

#endif //MLG_CS_EQUALTREE_H
