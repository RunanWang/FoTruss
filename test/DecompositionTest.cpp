//
// Created by 王润安 on 2024/3/22.
//
#include "../Algorithm/TrussOpt.h"
#include "testHeader.h"
#include "../constant.h"
#include "../Utils/logging.h"
#include "../Utils/Timer.h"


void testFoTrussDecomposition(MLGWithSchema *mlg, const filesystem::path &outPath) {
    auto timer1 = new Timer();
    timer1->startTimer();
    testFoTrussDecomposition(mlg, outPath, decompositionAlgoEnum(allD));
    timer1->endTimer();
    auto timer2 = new Timer();
    timer2->startTimer();
    testFoTrussDecomposition(mlg, outPath, decompositionAlgoEnum(plusD));
    timer2->endTimer();
    auto timer3 = new Timer();
    timer3->startTimer();
    testFoTrussDecomposition(mlg, outPath, decompositionAlgoEnum(basicD));
    timer3->endTimer();
    LOG(INFO) << fmt::format("Time for Decomposition* {:.4f}s, Decomposition+ {:.4f}s, and Decomposition {:.4f}s",
                             timer1->getTimerSecond(), timer2->getTimerSecond(), timer3->getTimerSecond());
}


void testFoTrussDecomposition(MLGWithSchema *mlg, const filesystem::path &outPath, int algorithm) {
    auto maxComboSize = LAYER_TYPE(pow(2, mlg->layersNum)) + mlg->layersNum * LAYER_TYPE(pow(2, mlg->layersNum - 1));
    auto result = new EDGE_TYPE *[maxComboSize];
    auto ci = new ComboIndex();
    if (algorithm == decompositionAlgoEnum(allD)) {
        LOG(INFO) << fmt::format("[Test] Time for Decomposition*");
        FoTrussDecompositionAll(mlg, result, ci);
    } else if (algorithm == decompositionAlgoEnum(plusD)) {
        LOG(INFO) << fmt::format("[Test] Time for Decomposition+");
        FoTrussDecompositionPlus(mlg, result, ci);
    } else {
        LOG(INFO) << fmt::format("[Test] Time for Decomposition");
        FoTrussDecompositionBasic(mlg, result, ci);
    }
    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        delete[] result[i];
    }
    delete[] result;
    delete ci;
}