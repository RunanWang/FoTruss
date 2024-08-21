//
// Created by 王润安 on 2023/12/28.
//

#include "UnionFind.h"
#include "algorithm"

using namespace std;

UnionFind::UnionFind(EDGE_TYPE size) {
    maxSize = size;
    parent = new EDGE_TYPE[size + 1];
    rank = new EDGE_TYPE[size + 1];
    pivot = new EDGE_TYPE[size + 1];
    for (EDGE_TYPE i = 0; i < size; i++) {
        parent[i] = i;
        rank[i] = 0;
        pivot[i] = i;
    }
}

UnionFind::~UnionFind() {
    delete[] parent;
    delete[] rank;
    delete[] pivot;
}

EDGE_TYPE UnionFind::find(EDGE_TYPE x) {
    return x == parent[x] ? x : (parent[x] = find(parent[x]));
}

void UnionFind::toUnion(EDGE_TYPE x1, EDGE_TYPE x2) {
    EDGE_TYPE f1 = find(x1);
    EDGE_TYPE f2 = find(x2);
    EDGE_TYPE newPivot = min(pivot[f1], pivot[f2]);
    if (rank[f1] > rank[f2]) {
        parent[f2] = f1;
        pivot[f1] = newPivot;
    } else {
        parent[f1] = f2;
        pivot[f2] = newPivot;
        if (rank[f1] == rank[f2])
            ++rank[f2];
    }
}

EDGE_TYPE UnionFind::getPivot(EDGE_TYPE item) {
    EDGE_TYPE father = find(item);
    return pivot[father];
}