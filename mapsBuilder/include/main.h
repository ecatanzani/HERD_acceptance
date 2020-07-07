#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>
#include <vector>

#include "anyoption.h"

#pragma once

extern void evaluateEventMap(
    const std::string accPath,
    const std::string telPath,
    const std::string outPath,
    const bool verbose,
    AnyOption &opt);

extern void init_data_maps(
    std::vector<std::vector<float>> &pixel_dataMap, 
    const long npix);
void write_final_maps(
    std::vector<std::vector<float>> &pixel_dataMap, 
    const std::string outputPath, 
    AnyOption &opt);   

#endif