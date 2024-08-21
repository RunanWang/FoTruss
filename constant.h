//
// Created by 王润安 on 2023/12/4.
//

#ifndef MLGWORK_CONSTANT_H
#define MLGWORK_CONSTANT_H

typedef unsigned int NODE_TYPE;
typedef unsigned int EDGE_TYPE;
typedef unsigned int LAYER_TYPE;

const float DEN_BETA = 1.0;

enum triangleConnectLevelEnum {
    basic = 1,
    shareEdge = 2
};

enum decompositionAlgoEnum {
    basicD = 1,
    plusD = 2,
    allD = 3
};

const int DEFAULT_DEC_ALGO_LEVEL = decompositionAlgoEnum(allD);

const int DEFAULT_TRIANGLE_CONNECT_LEVEL = triangleConnectLevelEnum(basic);

const LAYER_TYPE GREEDY_FOTRUSS_INDEX_MAX_COMBO_SIZE = 10;

const EDGE_TYPE MAX_PRUNE_TRUSSNESS = 0;

#endif //MLGWORK_CONSTANT_H
