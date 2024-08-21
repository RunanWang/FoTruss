//
// Created by 王润安 on 2023/12/20.
//

#include "Combo.h"
#include <fmt/format.h>
#include "../Utils/logging.h"

Combo::Combo(LAYER_TYPE layerNum) {
    this->layerNum = layerNum;
    this->lamb = 1;
    this->focusNum = 0;
    focusLayers = new LAYER_TYPE[layerNum];
}

void Combo::copyCombo(Combo *c) {
    this->lamb = c->lamb;
    this->focusNum = c->focusNum;
    this->layerNum = c->layerNum;
    for (LAYER_TYPE i = 0; i < focusNum; i++) {
        this->focusLayers[i] = c->focusLayers[i];
    }
}

Combo::Combo(Combo *c) {
    this->lamb = c->lamb;
    this->focusNum = c->focusNum;
    this->layerNum = c->layerNum;
    focusLayers = new LAYER_TYPE[layerNum];
    for (LAYER_TYPE i = 0; i < focusNum; i++) {
        this->focusLayers[i] = c->focusLayers[i];
    }
}

Combo::~Combo() {
    delete[] focusLayers;
}

string Combo::toString() const {
    string focus;
    for (LAYER_TYPE i = 0; i < focusNum; i++) {
        if (focusLayers[i] < layerNum) {
            focus += fmt::format("{},", focusLayers[i]);
        }
    }
    return fmt::format("Background{}-Focus[{}]", lamb, focus);
}

ComboIndex::ComboIndex() {
//    idToComboName = *new vector<string>();
//    comboNameToId = *new unordered_map<string, LAYER_TYPE>();
}

ComboIndex::~ComboIndex() {
//    idToComboName.clear();
//    comboNameToId.clear();
//    idToComboName = *new vector<string>();
//    comboNameToId = *new unordered_map<string, LAYER_TYPE>();
}

void ComboIndex::getOneHopFatherIds(Combo *childCombo, LAYER_TYPE *fatherIds, LAYER_TYPE &fatherSize) {
    auto fatherCombo = new Combo(childCombo->layerNum);
    string fatherName;
    fatherSize = 0;
    if (childCombo->lamb > 1 and childCombo->lamb != childCombo->focusNum) {
        fatherCombo->copyCombo(childCombo);
        fatherCombo->lamb -= 1;
        fatherName = fatherCombo->toString();
        if (comboNameToId.find(fatherName) != comboNameToId.end()) {
            fatherIds[fatherSize] = comboNameToId[fatherName];
            fatherSize++;
        }
    }
    if (childCombo->focusNum != 0) {
        for (LAYER_TYPE i = 0; i < childCombo->focusNum; i++) {
            fatherCombo->copyCombo(childCombo);
            fatherCombo->focusLayers[i] = fatherCombo->layerNum + 1;
            fatherName = fatherCombo->toString();
//            LOG(INFO) << fmt::format("FatherName of {}: {}", idToCombo[childId].toString(), fatherName);
            if (comboNameToId.find(fatherName) != comboNameToId.end()) {
                fatherIds[fatherSize] = comboNameToId[fatherName];
                fatherSize++;
            }
        }
    }
    delete fatherCombo;
}

void ComboIndex::insert(Combo *childCombo) {
    LAYER_TYPE id = idToComboName.size();
    comboNameToId[childCombo->toString()] = idToComboName.size();
    idToComboName.emplace_back(childCombo->toString());
    tomb.emplace_back(false);
    comboNum++;
    fathersId.resize(fathersId.size() + 1);
    auto fatherIds = new LAYER_TYPE[childCombo->layerNum + 2];            // 每个combo最多也就layersNum个父亲
    LAYER_TYPE fatherSize;
    getOneHopFatherIds(childCombo, fatherIds, fatherSize);
    for (LAYER_TYPE i = 0; i < fatherSize; i++) {
        LAYER_TYPE fatherId = fatherIds[i];
        fathersId[id].emplace_back(fatherId);
    }
}

void ComboIndex::deleteLastInsert(Combo *childCombo) {
    idToComboName.pop_back();
    comboNameToId.erase(childCombo->toString());
}

bool ComboIndex::allFatherInIndex(Combo *childCombo) {
    auto fatherCombo = new Combo(childCombo->layerNum);
    string fatherName;
    if (childCombo->lamb > 1 and childCombo->lamb != childCombo->focusNum) {
        fatherCombo->copyCombo(childCombo);
        fatherCombo->lamb -= 1;
        fatherName = fatherCombo->toString();
        if (comboNameToId.find(fatherName) == comboNameToId.end()) {
            return false;
        }
    }
    if (childCombo->focusNum > 0) {
        for (LAYER_TYPE i = 0; i < childCombo->focusNum; i++) {
            fatherCombo->copyCombo(childCombo);
            fatherCombo->focusLayers[i] = fatherCombo->layerNum + 1;
            fatherName = fatherCombo->toString();
            if (comboNameToId.find(fatherName) == comboNameToId.end()) {
                return false;
            }
        }
    }
    delete fatherCombo;
    return true;
}

void ComboIndex::setTomb(LAYER_TYPE comboId) {
    if (comboId < tomb.size()) {
        tomb[comboId] = true;
    }
}

LAYER_TYPE ComboIndex::prunedFatherComboNum() {
    LAYER_TYPE ret = 0;
    for (auto item: tomb) {
        if (item) {
            ret++;
        }
    }
    return ret;
}

LAYER_TYPE ComboIndex::prunedSonComboNum(LAYER_TYPE layerNum) const {
    auto maxComboSize = LAYER_TYPE(pow(2, layerNum)) + layerNum * LAYER_TYPE(pow(2, layerNum - 1)) - 1;
    return maxComboSize - idToComboName.size();
}