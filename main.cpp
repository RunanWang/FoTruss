#include <filesystem>
#include <fmt/format.h>
#include "gflags/gflags.h"
#include "glog/logging.h"

#include "Utils/Timer.h"
#include "Utils/Memory.h"
#include "Utils/logging.h"
#include "DataStructure/MLGWithSchema.h"
#include "test/testHeader.h"
#include "Algorithm/TrussOpt.h"
#include "constant.h"

using namespace std;

DEFINE_string(output_path, "./output", "output dir");
DEFINE_string(input_path, "./Dataset", "input graph file");
DEFINE_string(dataset, "homo", "input dataset name");
DEFINE_string(routine, "run", "mlg w/wo schema");
DEFINE_string(algo, "index", "execute algorithm.");

void check_dir(const filesystem::path &path) {
    if (not filesystem::exists(path)) {
        filesystem::create_directory(path);
    }
}

void routineSchema() {
    auto output_path = filesystem::path(FLAGS_output_path);
    output_path.append(FLAGS_dataset);
    auto input_file_path = filesystem::path(FLAGS_input_path).append(fmt::format("{}.txt", FLAGS_dataset));
    // load graph
    LOG(INFO) << fmt::format("Reading Graph from {}...", input_file_path.string());
    auto timer = new Timer();
    timer->startTimer();
    auto mlg = new MLGWithSchema();
    if (not filesystem::exists(input_file_path)) {
        LOG(ERROR) << fmt::format("No Such File: {}", input_file_path.string());
        return;
    }
    mlg->LoadFromFile(input_file_path);
    timer->endTimer();
    LOG(INFO) << fmt::format("Loading Phase: {:.4f} Seconds.", timer->getTimerSecond());
    LOG(INFO) << fmt::format("Peak Memory Usage in loading = {:.4f} MB", GetPeakRSSInMB());
    if (FLAGS_algo == "decomposition") {
        LOG(INFO) << "### FoTruss Decomposition ###";
        auto outPath = output_path.append(fmt::format("{}-fotruss.txt", FLAGS_dataset));
        FoTrussDecompositionOpt(mlg, output_path);
    } else if (FLAGS_algo == "index") {
        LOG(INFO) << "### Constructing FoTruss Index ###";
        auto outPath = output_path.append(fmt::format("{}-fotruss.txt", FLAGS_dataset));
        FoTrussIndex(mlg, output_path);
    }
    delete timer;
    delete mlg;
}

void routineTest() {
    auto input_path = filesystem::path(FLAGS_input_path);
    if (FLAGS_algo == "sota-gt") {
        if (FLAGS_dataset == "RM") {
            LOG(INFO) << fmt::format("[Exp=SOTA-GT-RM] Generating density and f1 of community for RM.");
            testRMCommunity(input_path);
        } else if (FLAGS_dataset == "terrorist") {
            LOG(INFO) << fmt::format("[Exp=SOTA-GT-terrorist] Generating density and f1 of community for terrorist.");
            testTerroristCommunity(input_path);
        } else {
            LOG(ERROR) << fmt::format("Dataset {} dose not have ground truth community.", FLAGS_dataset);
        }

    } else if (FLAGS_algo == "sota-random") {
        LOG(INFO)
                << fmt::format("[Exp=SOTA-RANDOM] Generating density and diameter of community for {}.", FLAGS_dataset);
        testCommunityByDensityUsingRandomQueries(input_path, FLAGS_dataset, 100);
    } else if (FLAGS_algo == "sota-given") {
        LOG(INFO)
                << fmt::format("[Exp=SOTA-GIVEN] Generating density and diameter of community for {}.", FLAGS_dataset);
        vector<NODE_TYPE> queries = {310, 326, 659, 676, 772, 789, 812, 993, 1005, 1055, 1068, 1162, 1904, 1998, 2029,
                                     2170, 2214, 2270, 2612, 2631, 2717, 2757, 2968, 2976, 2985, 2991, 3157, 3184, 3236,
                                     3255, 3314, 3328, 3344, 3552, 3713, 3839, 3926, 4031, 4075, 4089, 4154, 4157, 4254,
                                     4277, 4331, 4553, 4707, 4812, 4975, 5231, 5249, 5274, 5520, 5980, 5994, 6362, 6387,
                                     6446, 6701, 7149, 7509, 7661, 8301, 8391, 8729, 8867, 8868, 8993, 9081, 9362, 9363,
                                     9393, 9558, 9696, 9849, 9968, 9999, 10368, 10789, 10793, 11048, 11115, 11913,
                                     12505, 12695, 12857, 12858, 13033, 13164, 13276, 13390, 13490, 13648, 14405, 14604,
                                     14872, 14899, 16370, 17478, 17631};
        testCommunityByDensityUsingGivenQueries(input_path, FLAGS_dataset, queries);
    } else if (FLAGS_algo == "sota-single") {
        NODE_TYPE queryNode = 287;
        LOG(INFO) << fmt::format(
                "[Exp=SOTA-SINGLE] Generating density and diameter of community for {}, with query = {}.",
                FLAGS_dataset, queryNode);
        testCommunityByDensityForSpecificNode(input_path, FLAGS_dataset, queryNode);
    } else if (FLAGS_algo == "sota-noindex") {
        LOG(INFO) << fmt::format("[Exp=SOTA-NOINDEX] Generating density of community for {}.", FLAGS_dataset);
        vector<NODE_TYPE> queries = {3255};
        searchCommByTrussnessForGivenQueries(input_path, FLAGS_dataset, queries);
    } else if (FLAGS_algo == "time-dec") {
        LOG(INFO) << fmt::format("[EXP=DecTime] Time for Dec*, Dec+, and Dec for {}", FLAGS_dataset);
        auto output_path = filesystem::path(FLAGS_output_path);
        output_path.append(FLAGS_dataset);
        auto outPath = output_path.append(fmt::format("{}-fotruss.txt", FLAGS_dataset));
        auto input_file_path = filesystem::path(FLAGS_input_path).append(fmt::format("{}.txt", FLAGS_dataset));
        auto mlg = new MLGWithSchema();
        mlg->LoadFromFile(input_file_path);
        testFoTrussDecomposition(mlg, output_path);
    } else if (FLAGS_algo == "index-con") {
        LOG(INFO) << fmt::format("[EXP=IndexConstruction] Index construction time, MemUsage and statistics for {}",
                                 FLAGS_dataset);
        auto input_file_path = filesystem::path(FLAGS_input_path).append(fmt::format("{}.txt", FLAGS_dataset));
        auto mlg = new MLGWithSchema();
        mlg->LoadFromFile(input_file_path);
        testIndexConstruction(mlg);
    } else if (FLAGS_algo == "fotruss-index-greedy") {
        LOG(INFO) << "### Constructing Part of FoTruss Index (GREEDY) ###";
        vector<NODE_TYPE> queries = {93, 3747, 6918, 8741, 17388};
        vector<NODE_TYPE> testQueries = {211, 393, 465, 642, 692, 720, 733, 932, 941, 992, 1445, 1530, 1720, 1727, 1829,
                                         1971, 1996, 2077, 2131, 2375, 2483, 2633, 2749, 3003, 3016, 3054, 3151, 3206,
                                         3228, 3249, 3275, 3345, 3402, 3434, 3448, 3818, 3820, 3994, 4042, 4083, 4272,
                                         4408, 4495, 4516, 4536, 4702, 4803, 4845, 4936, 4967, 5020, 5255, 5362, 5579,
                                         5618, 6174, 6314, 6470, 6471, 6509, 6779, 6998, 7287, 7342, 7535, 7575, 7690,
                                         7724, 7771, 7998, 8335, 8403, 8488, 8640, 8957, 9000, 9360, 9692, 9752, 9845,
                                         10028, 10046, 10301, 10776, 11093, 11113, 11459, 12861, 12899, 12959, 13703,
                                         13952, 14319, 14362, 17461};
        LOG(INFO) << "";
        LOG(INFO) << "###### FULL";
        LOG(INFO) << "";
        testCommunityByDensityUsingGivenQueries(input_path, FLAGS_dataset, testQueries);
        auto mlg = new MLGWithSchema();
        auto input_file_path = filesystem::path(FLAGS_input_path).append(fmt::format("{}.txt", FLAGS_dataset));
        mlg->LoadFromFile(input_file_path);
        vector<NODE_TYPE> comboSize = {mlg->layersNum / 8, mlg->layersNum / 4, mlg->layersNum / 2, mlg->layersNum,
                                       2 * mlg->layersNum};
        for (auto cSize: comboSize) {
            LOG(INFO) << "";
            LOG(INFO) << "###### size = " << fmt::format("{}", cSize);
            LOG(INFO) << "";
            vector<string> combos;
            vector<EqualTree *> index;
            FoTrussIndexGreedy(mlg, cSize, queries, testQueries, combos, index);
        }
    }
}

int main(int argc, char *argv[]) {
    // init
    FLAGS_alsologtostderr = true; //设置日志消息除了日志文件之外是否去标准输出
    FLAGS_log_prefix = true; //设置日志前缀是否应该添加到每行输出
    google::ParseCommandLineFlags(&argc, &argv, true);
    google::InitGoogleLogging(argv[0]);
    FLAGS_log_dir = FLAGS_output_path;

    if (FLAGS_routine == "run") routineSchema();
    else if (FLAGS_routine == "test") routineTest();
    else
        LOG(ERROR) << fmt::format("Routine {} not supported", FLAGS_routine);

    LOG(INFO) << fmt::format("Peak Memory Usage in total = {:.4f} MB", GetPeakRSSInMB());
    return 0;
}

