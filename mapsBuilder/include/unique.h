#ifndef UNIQUE_H
#define UNIQUE_H

#include <string>

#include "anyoption.h"

extern const std::string uniqueOutFile(
    const std::string outputPath, 
    AnyOption &opt,
    const int binIdx);

#endif