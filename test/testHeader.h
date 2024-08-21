//
// Created by 王润安 on 2023/12/11.
//

#ifndef MLGWORK_TESTHEADER_H
#define MLGWORK_TESTHEADER_H

#include "filesystem"
#include "../DataStructure/MLGWithSchema.h"

void testRMCommunity(const filesystem::path &input_dir_path);

void testTerroristCommunity(const filesystem::path &input_dir_path);

void testCommunityByDensityUsingRandomQueries(const filesystem::path &input_dir_path, const string &dataset_name,
                                              NODE_TYPE queryNum);

void testCommunityByDensityForSpecificNode(const filesystem::path &input_dir_path, const string &dataset_name,
                                           NODE_TYPE query);

void testCommunityByDensityUsingGivenQueries(const filesystem::path &input_dir_path, const string &dataset_name,
                                             const vector<NODE_TYPE> &queries);

void searchCommByTrussnessForGivenQueries(const filesystem::path &input_dir_path, const string &dataset_name,
                                          const vector<NODE_TYPE> &queries);

void testFoTrussDecomposition(MLGWithSchema *mlg, const filesystem::path &outPath, int algorithm);

void testFoTrussDecomposition(MLGWithSchema *mlg, const filesystem::path &outPath);

void testIndexConstruction(MLGWithSchema *mlg);

#endif //MLGWORK_TESTHEADER_H
