#include "main.h"

int main(int argc, char** argv)
{
    AnyOption opt;

    opt.addUsage("Usage: ");
    opt.addUsage("");
    opt.addUsage(" -h  --help                                                   Prints this help");
    opt.addUsage(" -a  --acceptance     <path_to_acceptance_TFile>          (*) Input acceptance TFile");
    opt.addUsage(" -t  --telemetry      <path_to_telemetry_TFile>           (*) Input telemetry TFile");
    opt.addUsage(" -o  --output         <path_to_output_TFile>                  Output ROOT TFile");
    opt.addUsage(" -d  --outputDir      <path_to_output_TFile_dir>              Output ROOT TFile directory");
    opt.addUsage(" -v  --verbose                                                Verbose output");
    opt.addUsage("");

    opt.setFlag("help", 'h');
    opt.setOption("acceptance", 'a');
    opt.setOption("telemetry", 't');
    opt.setOption("output", 'o');
    opt.setOption("outputDir", 'd');
    opt.setFlag("verbose", 'v');

    opt.processCommandArgs(argc, argv);

    std::string accPath;
    std::string telPath;
    std::string outPath;
    bool verbose = false;

    if (!opt.hasOptions())
        opt.printUsage();

    if (opt.getFlag("help") || opt.getFlag('h'))
    {
        opt.printUsage();
        return 0;
    }
    if (opt.getValue("acceptance") || opt.getValue('a'))
        accPath = opt.getValue('a');
    if (opt.getValue("telemetry") || opt.getValue('t'))
        telPath = opt.getValue('t');
    if (opt.getValue("output") || opt.getValue('o'))
        outPath = opt.getValue('o');
    if (opt.getValue("outputDir") || opt.getValue('d'))
        outPath = opt.getValue('d');
    if (opt.getFlag("verbose") || opt.getFlag('v'))
        verbose = opt.getFlag('v');

    return 0;
}