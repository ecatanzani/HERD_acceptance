#include "unique.h"

#include <ctime>
#include <sstream>

const std::string uniqueOutFile(
    const std::string outputPath, 
    AnyOption &opt,
    const int binIdx)
{
    std::time_t ctime = std::time(0);
    std::stringstream fPath;
    if (opt.getValue("outputDir") || opt.getValue('d'))
        fPath << outputPath << "/mapsOutFile_energyBin_" << ctime << ".fits";
    else if (opt.getValue("output") || opt.getValue('o'))
        fPath << outputPath;
    else
        fPath << "mapsOutFile_energyBin_" << ctime << ".fits";

    return fPath.str();
}