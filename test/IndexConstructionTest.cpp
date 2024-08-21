//
// Created by 王润安 on 2024/4/23.
//

#include "../Algorithm/TrussOpt.h"
#include "testHeader.h"
#include "../constant.h"
#include "../Utils/logging.h"
#include "../Utils/Timer.h"

void testIndexConstruction(MLGWithSchema *mlg) {
    auto timer = new Timer();
    auto graphMemUsage = mlg->getMemUsage();
    auto maxComboSize = LAYER_TYPE(pow(2, mlg->layersNum)) + mlg->layersNum * LAYER_TYPE(pow(2, mlg->layersNum - 1));
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *ci = new ComboIndex();
    EDGE_TYPE totalSupernodeNum = 0;
    EDGE_TYPE totalSuperEdgeNum = 0;
    FoTrussDecompositionAll(mlg, trussness, ci);
    auto **index = new EqualTree *[ci->comboNum + 1];
    timer->startTimer();
    FoTrussIndexAll(mlg, index, ci, trussness);
    timer->endTimer();
    double totalIndexSize = 0;
    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        if (ci->tomb[i]) {
            continue;
        }
        totalIndexSize += index[i]->getMemUsage();
        totalSupernodeNum += index[i]->superNodeNum;
        for (const auto &item: index[i]->superNodeChild) {
            totalSuperEdgeNum += item.size();
        }

        delete index[i];
    }
    LOG(INFO) << fmt::format("Index Construction cost: {:.4f}s.", timer->getTimerSecond());
    LOG(INFO) << fmt::format("Index size: {:.4f}MB. Graph size: {:.4f}MB. {:.2f}%", totalIndexSize, graphMemUsage,
                             totalIndexSize / graphMemUsage * 100);
    LOG(INFO) << fmt::format("#SuperNodes = {}, #SuperEdges = {}", totalSupernodeNum, totalSuperEdgeNum);
    delete ci;
}