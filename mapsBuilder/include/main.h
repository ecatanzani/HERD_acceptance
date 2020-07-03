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
    const bool verbose);

#endif