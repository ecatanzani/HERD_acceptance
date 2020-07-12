#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "anyoption.h"

#include "TF1.h"
#include "TH1D.h"

#pragma once

extern void evaluateEventMap(
    const std::string evtPath,
    const std::string telPath,
    const std::string outPath,
    const bool verbose,
    AnyOption &opt,
    const TF1 &fitFunc,
    const std::shared_ptr<TH1D> acceptance,
    const std::shared_ptr<TH1D> live_time);

extern void init_data_maps(
    std::vector<std::vector<float>> &pixel_dataMap, 
    const long npix);
void write_final_maps(
    std::vector<std::vector<float>> &pixel_dataMap, 
    const std::string outputPath, 
    AnyOption &opt);   

#endif