//
// Created by 王润安 on 2023/12/5.
//

#ifndef MLGWORK_ALGO_H
#define MLGWORK_ALGO_H

#include "../constant.h"
#include "vector"
#include "set"
#include "../DataStructure/MLGWithSchema.h"

using namespace std;

NODE_TYPE findTopK(const NODE_TYPE *degreeList, LAYER_TYPE layerNum, LAYER_TYPE k);

/*
 * 求两个有序列表list1和list2的交集，长度分别为len1和len2，返回结果在ans中，ans的长度在return值里
 * */
NODE_TYPE intersect(const NODE_TYPE *list1, const NODE_TYPE *list2, NODE_TYPE len1, NODE_TYPE len2, NODE_TYPE *ans);


double f1Score(vector<NODE_TYPE> community, vector<NODE_TYPE> gt_community);


double density(const vector<NODE_TYPE> &subgraph, MLGWithSchema *mlg, double beta);


double weightedDensity(const vector<NODE_TYPE> &subgraph, MLGWithSchema *mlg, const vector<double> &weight);


void densityOfEachLayer(const vector<NODE_TYPE> &subgraph, MLGWithSchema *mlg, vector<double> &den);


NODE_TYPE diameter(const vector<NODE_TYPE> &subgraph, MLGWithSchema *mlg);


NODE_TYPE diameterMultilayer(const vector<NODE_TYPE> &subgraph, MLGWithSchema *mlg);

#endif //MLGWORK_ALGO_H
