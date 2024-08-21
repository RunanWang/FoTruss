//
// Created by 王润安 on 2024/1/2.
//

#ifndef MLG_CS_TRUSSOPT_H
#define MLG_CS_TRUSSOPT_H

#include <filesystem>
#include "../DataStructure/MLGWithSchema.h"
#include "../DataStructure/Combo.h"
#include "../DataStructure/EqualTree.h"

// Decomposition Algorithms
void FoTrussDecompositionOpt(MLGWithSchema *mlg, const filesystem::path &outPath);
void FoTrussDecompositionAll(MLGWithSchema *mlg, EDGE_TYPE **result, ComboIndex *comboIndex);
void FoTrussDecompositionPlus(MLGWithSchema *mlg, EDGE_TYPE **result, ComboIndex *comboIndex);
void FoTrussDecompositionBasic(MLGWithSchema *mlg, EDGE_TYPE **result, ComboIndex *comboIndex);

// Fotruss Index
void FoTrussIndex(MLGWithSchema *mlg, const filesystem::path &outPath);
void FoTrussIndexAll(MLGWithSchema *mlg, EqualTree **index, ComboIndex *ci, EDGE_TYPE **trussness);
void FoTrussIndexGreedy(MLGWithSchema *mlg, LAYER_TYPE maxComboSize, const vector<NODE_TYPE> &queries,
                        const vector<NODE_TYPE> &testQueries, vector<string> &combos, vector<EqualTree *> &index);
void FoTrussIndexGreedyFirstK(MLGWithSchema *mlg, LAYER_TYPE maxComboSize, const vector<NODE_TYPE> &testQueries);
void FoTrussIndexGreedyFirm(MLGWithSchema *mlg, LAYER_TYPE maxComboSize, const vector<NODE_TYPE> &testQueries);

// Triangle
void edgesTriangleConnectedByGivenEdgeBasic(MLGWithSchema *mlg, EDGE_TYPE edge, NODE_TYPE *ret, NODE_TYPE &size);
void TriangleConnectedTruss(MLGWithSchema *mlg, EDGE_TYPE *trussness, EqualTree *eqTree, Schema *usingSchema);

// Index Usage
EDGE_TYPE getNodeMaxTrussness(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode);
EDGE_TYPE getMaxTrussness(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode);
void queryNodeCommunity(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode, EDGE_TYPE trussnessBound,
                        vector<vector<NODE_TYPE>> &result);

// Storage
void saveTrussness(MLGWithSchema *mlg, EDGE_TYPE **trussness, const filesystem::path &outPath);

#endif //MLG_CS_TRUSSOPT_H
