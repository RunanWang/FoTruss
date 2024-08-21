//
// Created by 王润安 on 2023/12/20.
//

#ifndef MLG_CS_COMBO_H
#define MLG_CS_COMBO_H

#include "vector"
#include "../constant.h"
#include "string"
#include "unordered_map"

using namespace std;

class Combo {
public:
    LAYER_TYPE layerNum{0};
    LAYER_TYPE lamb{0};
    LAYER_TYPE focusNum{0};
    LAYER_TYPE *focusLayers{nullptr};

    Combo(LAYER_TYPE layerNum);

    Combo(Combo *c);

    void copyCombo(Combo *c);

    ~Combo();

    string toString() const;
};

class ComboIndex {
public:
    // 这里我想的解决方案是id是唯一的，删除时只插入一个墓碑标记
    vector<string> idToComboName;
    unordered_map<string, LAYER_TYPE> comboNameToId;
    LAYER_TYPE comboNum = 0;
    LAYER_TYPE rootComboId = 0;
    vector<bool> tomb;
    vector<vector<LAYER_TYPE>> fathersId;

    ComboIndex();

    ~ComboIndex();

    bool allFatherInIndex(Combo *childCombo);

    void getOneHopFatherIds(Combo *childCombo, LAYER_TYPE *fatherIds, LAYER_TYPE &fatherSize);

    void insert(Combo *childCombo);

    void deleteLastInsert(Combo *childCombo);

    void setTomb(LAYER_TYPE comboId);

    LAYER_TYPE prunedSonComboNum(LAYER_TYPE layerNum) const;

    LAYER_TYPE prunedFatherComboNum();
};


#endif //MLG_CS_COMBO_H
