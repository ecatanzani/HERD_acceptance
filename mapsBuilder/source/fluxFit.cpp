#include "fluxFit.h"

TF1 GetFluxFittedFunc(const std::string fluxPath)
{
    TFile* inFile = TFile::Open(fluxPath.c_str(), "READ");
    if (inFile->IsZombie())
    {
        std::cerr << "\n\nError reading input file: " << fluxPath << std::endl;
        exit(100);
    }

    std::unique_ptr<TMultiGraph> mg = std::make_unique<TMultiGraph>();
    mg->SetTitle("All electron flux (AMS, DAMPE, CALET); Energy(GeV); Flux");
    mg->SetName("flux_mg");

    for (auto&& keyAsObj : *inFile->GetListOfKeys())
    {
        auto key = static_cast<TKey*>(keyAsObj);
        if (!strcmp(key->GetClassName(), "TGraphAsymmErrors"))
            mg->Add(static_cast<TGraphAsymmErrors *>(key->ReadObj()));
    }

    TF1 fitFuncE3("fitFuncE3", "pow(log10(x),3) * pow(10, [0] + [1]*log10(x) + [2]*pow(log10(x),2) +[3]*pow(log10(x),3) + [4]*pow(log10(x),4)) ", 5, 1e+4);
    fitFuncE3.SetNpx(10000);
    mg->Fit("fitFuncE3", "IR");

    TF1 fitFunc("fitFunc", "(pow(log10(x),3) * pow(10, [0] + [1]*log10(x) + [2]*pow(log10(x),2) +[3]*pow(log10(x),3) + [4]*pow(log10(x),4))) / pow(x,3)", 5, 1e+4);
    for (int nPar=0; nPar<fitFuncE3.GetNpar(); ++nPar)
        fitFunc.SetParameter(nPar, fitFuncE3.GetParameter(nPar));

    TFile outFile("flux_fit.root", "RECREATE");
    if (outFile.IsZombie())
    {
        std::cerr << "\n\nError writing output file" << std::endl;
        exit(100);
    }
    
    mg->Write();
    fitFuncE3.Write();
    fitFunc.Write();

    outFile.Close();

    return fitFunc;
}