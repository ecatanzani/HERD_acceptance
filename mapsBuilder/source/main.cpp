#include "main.h"
#include "fluxFit.h"
#include "acceptance_reader.h"
#include "compute_livetime.h"
#include "DAMPE_acquisition_rate.h"

int main(int argc, char** argv)
{
    AnyOption opt;

    opt.addUsage("Usage: ");
    opt.addUsage("");
    opt.addUsage(" -h  --help                                                   Prints this help");
    opt.addUsage(" -a  --acceptance     <path_to_acceptance_TFile>          (*) Input acceptance TFile");
    opt.addUsage(" -e  --event          <path_to_event_dist_TFile>          (*) Input angular event distribution TFile");
    opt.addUsage(" -f  --flux           <path_to_flux_TFile>                (*) Input flux TFile");
    opt.addUsage(" -t  --telemetry      <path_to_telemetry_TFile>           (*) Input telemetry TFile");
    opt.addUsage(" -r  --rate           <path_to_acq_rate_TFile>            (*) Input DAMPE acquisition rate TFile");
    opt.addUsage(" -o  --output         <path_to_output_TFile>                  Output ROOT TFile");
    opt.addUsage(" -d  --outputDir      <path_to_output_TFile_dir>              Output ROOT TFile directory");
    opt.addUsage(" -v  --verbose                                                Verbose output");
    opt.addUsage("");

    opt.setFlag("help", 'h');
    opt.setOption("acceptance", 'a');
    opt.setOption("event", 'e');
    opt.setOption("flux", 'f');
    opt.setOption("telemetry", 't');
    opt.setOption("rate", 'r');
    opt.setOption("output", 'o');
    opt.setOption("outputDir", 'd');
    opt.setFlag("verbose", 'v');

    opt.processCommandArgs(argc, argv);

    std::string accPath;
    std::string evtPath;
    std::string fluxPath;
    std::string telPath;
    std::string acqPath;
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
    if (opt.getValue("event") || opt.getValue('e'))
        evtPath = opt.getValue('e');
    if (opt.getValue("flux") || opt.getValue('f'))
        fluxPath = opt.getValue('f');
    if (opt.getValue("telemetry") || opt.getValue('t'))
        telPath = opt.getValue('t');
     if (opt.getValue("rate") || opt.getValue('r'))
        acqPath = opt.getValue('r');
    if (opt.getValue("output") || opt.getValue('o'))
        outPath = opt.getValue('o');
    if (opt.getValue("outputDir") || opt.getValue('d'))
        outPath = opt.getValue('d');
    if (opt.getFlag("verbose") || opt.getFlag('v'))
        verbose = opt.getFlag('v');

    // Extract all electron flux fit function
    auto fitFunc = GetFluxFittedFunc(fluxPath);

    // Extract HERD acceptance
    auto acceptance = GetAcceptanceHisto(accPath);
    
    // Extract DAMPE acquisition rate
    auto dampe_acquition_rate = GetDAMPEAcquisitionRate(acqPath);

    // Compute HERD livetime
    auto live_time = GetLiveTime(dampe_acquition_rate);
    
    
    evaluateEventMap(
        evtPath,
        telPath,
        outPath,
        verbose,
        opt,
        fitFunc,
        acceptance,
        live_time);
        
    return 0;
}