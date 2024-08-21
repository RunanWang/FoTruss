//
// Created by 王润安 on 2023/12/11.
//

#include "Bucket.h"

Bucket::Bucket(T maxItem, T maxBucket) {
    bucketSize = maxBucket + 1;
    itemSize = maxItem;
    data = new T *[bucketSize];
    length = new T[bucketSize];
    itemOffset = new T[itemSize];
    for (T i = 0; i < itemSize; i++) {
        itemOffset[i] = 0;
    }
    for (T i = 0; i < bucketSize; i++) {
        length[i] = 0;
        data[i] = nullptr;
    }
}

Bucket::~Bucket() {
    for (T i = 0; i < bucketSize; i++) {
        delete[] data[i];
    }
    delete[] data;
    delete[] length;
    delete[] itemOffset;
}

void Bucket::construct(const T *itemList) {
    for (T i = 0; i < itemSize; i++) {
        length[itemList[i]]++;
    }
    T aboveSum = 0;
    for (T i = bucketSize - 1; i != 0; i--) {
        aboveSum += length[i];
        data[i] = new T[aboveSum];
        length[i] = 0;
    }
    aboveSum += length[0];
    data[0] = new T[aboveSum];
    length[0] = 0;
    for (T i = 0; i < itemSize; i++) {
        data[itemList[i]][length[itemList[i]]] = i;
        itemOffset[i] = length[itemList[i]];
        length[itemList[i]]++;
    }
}

T *Bucket::getBucket(T bucketId) const {
    return data[bucketId];
}

void Bucket::changeItemBucket(T item, T oldBucket, T toBucket) {
    if (toBucket == oldBucket) return;
    // 先从old中移除
    if (length[oldBucket] == 1) {
        length[oldBucket] = 0;
    } else {
        T swapNode = data[oldBucket][length[oldBucket] - 1];
        data[oldBucket][itemOffset[item]] = swapNode;
        itemOffset[swapNode] = itemOffset[item];
        length[oldBucket] -= 1;
    }
    // 加入新bucket
    data[toBucket][length[toBucket]] = item;
    itemOffset[item] = length[toBucket];
    length[toBucket]++;
}

DynamicBucket::DynamicBucket(T maxItem, T maxBucket) {
    bucketSize = maxBucket + 1;
    itemSize = maxItem;
    data = new vector<T>[bucketSize];
    itemOffset = new T[itemSize];
    for (T i = 0; i < itemSize; i++) {
        itemOffset[i] = 0;
    }
    for (T i = 0; i < bucketSize; i++) {
        data[i] = *new vector<T>();
    }
}

DynamicBucket::~DynamicBucket() {
    delete[] data;
    delete[] itemOffset;
}

void DynamicBucket::construct(const T *itemList) {
    for (T i = 0; i < itemSize; i++) {
        itemOffset[i] = data[itemList[i]].size();
        data[itemList[i]].emplace_back(i);
    }
}

void DynamicBucket::changeItemBucket(T item, T oldBucket, T toBucket) {
    if (toBucket == oldBucket) return;
    // 先从old中移除
    if (data[oldBucket].size() == 1) {
        data[oldBucket].pop_back();
    } else {
        T swapNode = data[oldBucket][data[oldBucket].size() - 1];
        data[oldBucket][itemOffset[item]] = swapNode;
        itemOffset[swapNode] = itemOffset[item];
        data[oldBucket].pop_back();
    }
    // 加入新bucket
    itemOffset[item] = data[toBucket].size();
    data[toBucket].emplace_back(item);
}