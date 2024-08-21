//
// Created by 王润安 on 2023/12/11.
//

#ifndef MLGWORK_BUCKET_H
#define MLGWORK_BUCKET_H

#include "vector"

using namespace std;

typedef unsigned int T;

class Bucket {
public:

    T **data;
    T *length;
    T *itemOffset;
    T bucketSize;
    T itemSize;

    Bucket(T maxItem, T maxBucket);

    ~Bucket();

    void construct(const T *itemList);

    T *getBucket(T bucketId) const;

    void changeItemBucket(T item, T oldBucket, T toBucket);
};

class DynamicBucket {
public:

    vector<T> *data;
    T *itemOffset;
    T bucketSize;
    T itemSize;

    DynamicBucket(T maxItem, T maxBucket);

    ~DynamicBucket();

    void construct(const T *itemList);

    void changeItemBucket(T item, T oldBucket, T toBucket);
};


#endif //MLGWORK_BUCKET_H
