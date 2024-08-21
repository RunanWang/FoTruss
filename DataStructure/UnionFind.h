//
// Created by 王润安 on 2023/12/28.
//

#ifndef MLG_CS_UNIONFIND_H
#define MLG_CS_UNIONFIND_H

#include "../constant.h"

class UnionFind {
// 这里我们是带有pivot的UF，需要输入的item的id带顺序
// id更小的item的trussness更小
// parent-自己的上一个节点的id，等于自己说明是根结点，初始时每个点都是根结点
// rank-用于维护树高
// pivot-当前这个连通部分序号最小的点ID，在初始化时每个点的pivot等于自己的ID，当union时，根结点会保留更小的ID作为自己的pivot
// 只有根结点的pivot有意义，发生union时并不会更新全部连通分量的pivot
public:
    EDGE_TYPE maxSize;
    EDGE_TYPE *parent;
    EDGE_TYPE *rank;
    EDGE_TYPE *pivot;

    UnionFind(EDGE_TYPE size);

    ~UnionFind();

    EDGE_TYPE find(EDGE_TYPE item);

    void toUnion(EDGE_TYPE x1, EDGE_TYPE x2);

    EDGE_TYPE getPivot(EDGE_TYPE item);
};


#endif //MLG_CS_UNIONFIND_H
