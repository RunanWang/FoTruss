//
// Created by 王润安 on 2024/1/11.
//

#include "testHeader.h"
#include "../Algorithm/TrussOpt.h"
#include "../Algorithm/Algo.h"
#include "../Utils/logging.h"
#include "../Utils/Timer.h"
#include "../constant.h"
#include "stdlib.h"
#include "algorithm"

double
evalGroundTruthCommunityByDen(MLGWithSchema *mlg, EqualTree **index, ComboIndex *ci, NODE_TYPE query, double &bestden,
                              vector<NODE_TYPE> &gt) {
    double bestf1 = 0;
    bestden = 0;
    LAYER_TYPE bestComboId;
    vector<vector<NODE_TYPE>> qus;
    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        if (ci->tomb[i]) {
            continue;
        }
        queryNodeCommunity(mlg, index[i], query, getMaxTrussness(mlg, index[i], query), qus);
        for (const auto &qu: qus) {
            double f1 = f1Score(qu, gt);
            double den = density(qu, mlg, DEN_BETA);
            if (den > bestden) {
                bestComboId = i;
                bestf1 = f1;
                bestden = den;
            }
        }
    }
//    LOG(INFO) << fmt::format("Best community of vertex {} is {}", query, ci->idToComboName[bestComboId]);
    return bestf1;
}


double
evalGroundTruthCommunityByF1(MLGWithSchema *mlg, EqualTree **index, ComboIndex *ci, NODE_TYPE query, double &bestden,
                             vector<NODE_TYPE> &gt) {
    double bestf1 = 0;
    bestden = 0;
    EDGE_TYPE best_k = 0;
    LAYER_TYPE bestComboId;
    vector<vector<NODE_TYPE>> qus;
    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        if (ci->tomb[i]) {
            continue;
        }
        EDGE_TYPE k = getMaxTrussness(mlg, index[i], query);
        queryNodeCommunity(mlg, index[i], query, k, qus);
        for (const auto &qu: qus) {
            double f1 = f1Score(qu, gt);
            double den = density(qu, mlg, DEN_BETA);
            if (f1 > bestf1) {
                bestComboId = i;
                bestf1 = f1;
                bestden = den;
                best_k = k;
            }
        }
    }
    LOG(INFO) << fmt::format("Best community of vertex {} is {}-{}", query, ci->idToComboName[bestComboId], best_k);
    return bestf1;
}

double
evalGroundTruthCommunity(MLGWithSchema *mlg, EqualTree **index, ComboIndex *ci, NODE_TYPE query, double &bestden,
                         vector<NODE_TYPE> &gt) {
//    return evalGroundTruthCommunityByDen(mlg, index, ci, query, bestden, gt);
    return evalGroundTruthCommunityByF1(mlg, index, ci, query, bestden, gt);
//    return evalGroundTruthCommunityByF1Single(mlg, index, ci, query, bestden, gt);
}

void evalCommunityByDensity(MLGWithSchema *mlg, EqualTree **index, ComboIndex *ci, NODE_TYPE query, double &bestden,
                            double &bestdiam, Timer *timer) {
    bestden = 0;
    bestdiam = 0;
    LAYER_TYPE bestComboId = 1;
    vector<vector<NODE_TYPE>> qus;
    vector<NODE_TYPE> bestqu;
    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        if (ci->tomb[i]) {
            continue;
        }
        // LOG(INFO) << fmt::format("Evaluating query = {}, combo = {}", query, ci->idToComboName[i]);
        timer->startTimer();
        auto trussnessBound = getMaxTrussness(mlg, index[i], query);
        queryNodeCommunity(mlg, index[i], query, trussnessBound, qus);
//        queryNodeCommunity(mlg, index[i], query, 1, qu);
        timer->endTimer();
        for (const auto &qu: qus) {
            double den = density(qu, mlg, DEN_BETA);
//            LOG(INFO) << fmt::format("Combo = {}-{}, den = {}, community = [{}]", ci->idToComboName[i], trussnessBound,
//                                     den, fmt::join(qu.begin(), qu.end(), ", "));
            if (den > bestden) {
                bestComboId = i;
                bestden = den;
                bestqu.assign(qu.begin(), qu.end());
            }
        }
    }
    LOG(INFO) << fmt::format("Best community of vertex {} is {} with den = {}. Calculating diameter...", query,
                             ci->idToComboName[bestComboId], bestden);
//    bestdiam = double(diameterMultilayer(bestqu, mlg));
//    LOG(INFO) << fmt::format("Diameter = {:.4f}/{}", bestdiam, bestqu.size());
    bestdiam = 1;
}


void
evalCommunityByDensitySingle(MLGWithSchema *mlg, EqualTree **index, ComboIndex *ci, NODE_TYPE query, double &bestden,
                             double &bestdiam, Timer *timer) {
    bestden = 0;
    bestdiam = 0;
    LAYER_TYPE bestComboId = 1;
    vector<vector<NODE_TYPE>> qus;
    vector<NODE_TYPE> bestqu;
    for (LAYER_TYPE i = 1; i < ci->comboNum; i++) {
        if (ci->tomb[i]) {
            continue;
        }
        // LOG(INFO) << fmt::format("Evaluating query = {}, combo = {}", query, ci->idToComboName[i]);
        timer->startTimer();
        auto trussnessBound = getMaxTrussness(mlg, index[i], query);
        queryNodeCommunity(mlg, index[i], query, trussnessBound, qus);
//        queryNodeCommunity(mlg, index[i], query, 1, qu);
        timer->endTimer();
        for (const auto &qu: qus) {
            double den = density(qu, mlg, DEN_BETA);
//            LOG(INFO) << fmt::format("Combo = {}-{}, den = {}, community = [{}]", ci->idToComboName[i], trussnessBound,
//                                     den, fmt::join(qu.begin(), qu.end(), ", "));
            if (den > bestden) {
                bestComboId = i;
                bestden = den;
                bestqu.assign(qu.begin(), qu.end());
            }
        }
        break;
    }
    LOG(INFO) << fmt::format("Best community of vertex {} is {} with den = {}. Calculating diameter...", query,
                             ci->idToComboName[bestComboId], bestden);
    bestdiam = double(diameterMultilayer(bestqu, mlg));
    LOG(INFO) << fmt::format("Diameter = {:.4f}/{}", bestdiam, bestqu.size());
    if (bestqu.size() == 0) {
        bestdiam = 0;
    }
//    bestdiam = 1;
}


void evalCommunityByWDensity(MLGWithSchema *mlg, EqualTree **index, ComboIndex *ci, NODE_TYPE query, double &bestden,
                             double &bestdiam, Timer *timer) {
    bestden = 0;
    bestdiam = 0;
    vector<double> w = {1, 1, 1, 1, 4, 1, 1};
    LAYER_TYPE bestComboId = 1;
    vector<vector<NODE_TYPE>> qus;
    vector<NODE_TYPE> bestqu;
    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        if (ci->tomb[i]) {
            continue;
        }
        // LOG(INFO) << fmt::format("Evaluating query = {}, combo = {}", query, ci->idToComboName[i]);
        timer->startTimer();
        auto trussnessBound = getMaxTrussness(mlg, index[i], query);
        queryNodeCommunity(mlg, index[i], query, trussnessBound, qus);
//        queryNodeCommunity(mlg, index[i], query, 1, qu);
        timer->endTimer();
        for (const auto &qu: qus) {
            double den = weightedDensity(qu, mlg, w);
//            LOG(INFO) << fmt::format("Combo = {}-{}, den = {}, community = [{}]", ci->idToComboName[i], trussnessBound,
//                                     den, fmt::join(qu.begin(), qu.end(), ", "));
            if (den > bestden) {
                bestComboId = i;
                bestden = den;
                bestqu.assign(qu.begin(), qu.end());
            }
        }
    }
    LOG(INFO) << fmt::format("Best community of vertex {} is {} with den = {}. Calculating diameter...", query,
                             ci->idToComboName[bestComboId], bestden);
//    bestdiam = double(diameterMultilayer(bestqu, mlg));
//    LOG(INFO) << fmt::format("Diameter = {:.4f}/{}", bestdiam, bestqu.size());
    bestdiam = 1;
}


void testCommunityByDensityUsingRandomQueries(const filesystem::path &input_dir_path, const string &dataset_name,
                                              NODE_TYPE queryNum) {
    srand(1);
    auto path = filesystem::path(input_dir_path);
    auto graph_path = path.append(fmt::format("{}.txt", dataset_name));
    auto mlg = new MLGWithSchema();
    mlg->LoadFromFile(graph_path);
    auto maxComboSize = LAYER_TYPE(pow(2, mlg->layersNum)) + mlg->layersNum * LAYER_TYPE(pow(2, mlg->layersNum - 1));
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *ci = new ComboIndex();
    FoTrussDecompositionAll(mlg, trussness, ci);
    auto **index = new EqualTree *[ci->comboNum + 1];
    FoTrussIndexAll(mlg, index, ci, trussness);
    double sumDen = 0.0;
    double sumDiam = 0.0;
    double den, diam;
    auto timer = new Timer();
    set<NODE_TYPE> randomQueries;
    LOG(INFO) << fmt::format("Generating list of nodes satisfying condition {}...", ci->idToComboName[ci->rootComboId]);
    vector<NODE_TYPE> selectFromList;
    for (NODE_TYPE n = 0; n < mlg->nodesNum; n++) {
        if (getMaxTrussness(mlg, index[ci->rootComboId], n) >= 1) {
            selectFromList.emplace_back(n);
        }
    }
    LOG(INFO) << fmt::format("Selecting {} Nodes from {} candidates...", queryNum, selectFromList.size());
    for (int i = 0; i < queryNum; i++) {
        NODE_TYPE id = rand() % selectFromList.size();
        while (randomQueries.find(selectFromList[id]) != randomQueries.end()) {
            id = rand() % selectFromList.size();
        }
        randomQueries.insert(selectFromList[id]);
    }
    LOG(INFO) << fmt::format("Generated Queries = [{}]", fmt::join(randomQueries, ", "));
    for (auto q: randomQueries) {
        LOG(INFO) << fmt::format("Finding community for query vertex {}", q);
        evalCommunityByDensity(mlg, index, ci, q, den, diam, timer);
        sumDen += den;
        sumDiam += diam;
    }
    double avgDen = sumDen / double(queryNum);
    double avgDiam = sumDiam / double(queryNum);
    double avgQueryTime = timer->getTimerSecond() / double(queryNum);
    LOG(INFO) << fmt::format("Avg Density of {} is {:.4f}, Diameter is {}, cost time = {:.4f}s.", dataset_name, avgDen,
                             avgDiam, avgQueryTime);

    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        delete index[i];
    }
    delete ci;
}


void testCommunityByDensityUsingGivenQueries(const filesystem::path &input_dir_path, const string &dataset_name,
                                             const vector<NODE_TYPE> &queries) {
    auto path = filesystem::path(input_dir_path);
    auto graph_path = path.append(fmt::format("{}.txt", dataset_name));
    auto mlg = new MLGWithSchema();
    mlg->LoadFromFile(graph_path);
    auto maxComboSize = LAYER_TYPE(pow(2, mlg->layersNum)) + mlg->layersNum * LAYER_TYPE(pow(2, mlg->layersNum - 1));
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *ci = new ComboIndex();
    FoTrussDecompositionAll(mlg, trussness, ci);
    auto **index = new EqualTree *[ci->comboNum + 1];
    FoTrussIndexAll(mlg, index, ci, trussness);
    double sumDen = 0.0;
    double sumDiam = 0.0;
    NODE_TYPE havingDiam = 0;
    double den, diam;
    auto timer = new Timer();
    LOG(INFO) << fmt::format("Finding Communities for q in Queries = [{}]", fmt::join(queries, ", "));
    for (auto q: queries) {
        LOG(INFO) << fmt::format("Finding community for query vertex {}...", q);
//        evalCommunityByDensity(mlg, index, ci, q, den, diam, timer);
        evalCommunityByDensitySingle(mlg, index, ci, q, den, diam, timer);
        sumDen += den;
        sumDiam += diam;
        if (diam != 0) {
            havingDiam++;
        }
    }
    double avgDen = sumDen / double(queries.size());
    double avgDiam = sumDiam / double(havingDiam);
    double avgQueryTime = timer->getTimerSecond() / double(queries.size());
    LOG(INFO) << fmt::format("Avg Density of {} is {:.4f}, Diameter is {}, cost time = {:.4f}s.", dataset_name, avgDen,
                             avgDiam, avgQueryTime);

    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        delete index[i];
    }
    delete ci;
}


void testCommunityByDensityForSpecificNode(const filesystem::path &input_dir_path, const string &dataset_name,
                                           NODE_TYPE query) {
    auto path = filesystem::path(input_dir_path);
    auto graph_path = path.append(fmt::format("{}.txt", dataset_name));
    auto mlg = new MLGWithSchema();
    mlg->LoadFromFile(graph_path);
    auto maxComboSize = LAYER_TYPE(pow(2, mlg->layersNum)) + mlg->layersNum * LAYER_TYPE(pow(2, mlg->layersNum - 1));
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *ci = new ComboIndex();
    FoTrussDecompositionAll(mlg, trussness, ci);
    auto **index = new EqualTree *[ci->comboNum + 1];
    FoTrussIndexAll(mlg, index, ci, trussness);
    double den, diam;
    auto timer = new Timer();
    LOG(INFO) << fmt::format("Finding community for each query vertex...");
    evalCommunityByDensity(mlg, index, ci, query, den, diam, timer);

    LOG(INFO) << fmt::format("Avg Density of {} is {:.4f}, Diameter is {}, cost time = {:.4f}s.", dataset_name, den,
                             diam, timer->getTimerSecond());

    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        delete index[i];
    }
    delete ci;
}


void testRMCommunity(const filesystem::path &input_dir_path) {
    auto path = filesystem::path(input_dir_path);
    auto RM_path = path.append(fmt::format("RM.txt"));
    auto mlg = new MLGWithSchema();
    mlg->LoadFromFile(RM_path);
    auto maxComboSize = LAYER_TYPE(pow(2, mlg->layersNum)) + mlg->layersNum * LAYER_TYPE(pow(2, mlg->layersNum - 1));
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *ci = new ComboIndex();
    FoTrussDecompositionAll(mlg, trussness, ci);
    auto **index = new EqualTree *[ci->comboNum + 1];
    FoTrussIndexAll(mlg, index, ci, trussness);
    double sumF1 = 0.0;
    double sumDen = 0.0;
    double den;
    auto rmgt1 = new vector<NODE_TYPE>{12, 36, 38, 39, 40, 41, 43, 45, 46, 47, 48, 51, 54, 58, 60, 62, 64, 65, 68,
                                       69, 72, 73, 76, 80, 82, 84, 85};
    auto rmgt2 = new vector<NODE_TYPE>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
                                       23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 37, 42, 44, 49, 50, 52,
                                       53, 55, 56, 57, 59, 61, 63, 66, 67, 70, 71, 74, 75, 77, 78, 79, 81, 83, 86,
                                       87, 88, 89, 90, 91};
    auto timer = new Timer();
    timer->startTimer();
    LOG(INFO) << fmt::format("Finding community for each query vertex...");
    for (auto q: *rmgt1) {
        sumF1 += evalGroundTruthCommunity(mlg, index, ci, q, den, *rmgt1);
        sumDen += den;
    }
    for (auto q: *rmgt2) {
        sumF1 += evalGroundTruthCommunity(mlg, index, ci, q, den, *rmgt2);
        sumDen += den;
    }
    timer->endTimer();
    double avgF1 = sumF1 / 91.0;
    double avgDen = sumDen / 91.0;
    double avgQueryTime = timer->getTimerSecond() / 91.0;
    LOG(INFO) << fmt::format("Avg F1 of RM is {:.4f}, avg Den is {:.4f}, cost time = {:.4f}s.", avgF1, avgDen,
                             avgQueryTime);

    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        delete index[i];
    }
    delete ci;
}


void testTerroristCommunity(const filesystem::path &input_dir_path) {
    auto path = filesystem::path(input_dir_path);
    auto RM_path = path.append(fmt::format("terrorist.txt"));
    auto mlg = new MLGWithSchema();
    mlg->LoadFromFile(RM_path);
    auto maxComboSize = LAYER_TYPE(pow(2, mlg->layersNum)) + mlg->layersNum * LAYER_TYPE(pow(2, mlg->layersNum - 1));
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *ci = new ComboIndex();
    FoTrussDecompositionAll(mlg, trussness, ci);
    auto **index = new EqualTree *[ci->comboNum + 1];
    FoTrussIndexAll(mlg, index, ci, trussness);
    double sumF1 = 0.0;
    double sumDen = 0.0;
    double den;
    auto rmgt1 = new vector<NODE_TYPE>{3, 18, 24, 36, 38, 39, 44, 48, 50, 55, 56, 57, 63, 64};
    auto rmgt2 = new vector<NODE_TYPE>{1, 2, 9, 12, 14, 15, 16, 19, 29, 30, 35, 62, 66};
    auto rmgt3 = new vector<NODE_TYPE>{4, 20, 21, 27, 42, 47, 53, 60, 65};
    auto rmgt4 = new vector<NODE_TYPE>{5, 6, 7, 8, 10, 11, 13, 17, 22, 23, 25, 26, 28, 31, 33, 34, 40, 46, 49, 51, 52,
                                       54, 58, 59, 61, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79};
    auto rmgt5 = new vector<NODE_TYPE>{32, 37, 41, 43, 45};
    auto timer = new Timer();
    timer->startTimer();
    LOG(INFO) << fmt::format("Finding community for each query vertex...");
    for (auto q: *rmgt1) {
        sumF1 += evalGroundTruthCommunity(mlg, index, ci, q, den, *rmgt1);
        sumDen += den;
    }
    for (auto q: *rmgt2) {
        sumF1 += evalGroundTruthCommunity(mlg, index, ci, q, den, *rmgt2);
        sumDen += den;
    }
    for (auto q: *rmgt3) {
        sumF1 += evalGroundTruthCommunity(mlg, index, ci, q, den, *rmgt3);
        sumDen += den;
    }
    for (auto q: *rmgt4) {
        sumF1 += evalGroundTruthCommunity(mlg, index, ci, q, den, *rmgt4);
        sumDen += den;
    }
    for (auto q: *rmgt5) {
        sumF1 += evalGroundTruthCommunity(mlg, index, ci, q, den, *rmgt5);
        sumDen += den;
    }
    timer->endTimer();
    double avgF1 = sumF1 / 79.0;
    double avgDen = sumDen / 79.0;
    double avgQueryTime = timer->getTimerSecond() / 79.0;
    LOG(INFO)
            << fmt::format("Avg F1 of terrorist is {:.4f}, avg density is {:.4f}, cost time = {:.4f}s.", avgF1, avgDen,
                           avgQueryTime);

    for (LAYER_TYPE i = 0; i < ci->comboNum; i++) {
        delete index[i];
    }
    delete ci;
}
