#include "Utils/Plotting/PlottingManager.h"
#include <ROOT/RDataFrame.hxx>
#include "TFile.h"
#include "TChain.h"
#include <omp.h>

//This aims to take an input FD chain and produce a file containing 1d and 2d
//histograms with and without a reactor reweight applied
//to allow easy plotting, we need an input config telling us the binning to use
//in each parameter for and list the parameters to plot
//This is explicitly NOT a plotting script and never should be for performance reasons

M3::Plotting::PlottingManager *man;

void MakePosteriorHistograms(){

    TChain *inf = new TChain(man->getOption<std::string>("posteriorTreeName").c_str());
    auto infiles = man->getFileNames();
    for(auto file : infiles) inf->AddFile(file.c_str());
    int burnin = man->getOption<int>("burnin");
    YAML::Node woRCbinning = man->getOption("woRCparameterBinning");
    YAML::Node wRCbinning = man->getOption("wRCparameterBinning");
    auto reweights = man->getOption<std::vector<std::string>>("reweightBranches");
    auto parameter_names = man->getOption<std::vector<std::string>>("parameterNames");
    auto mass_order_parameter = man->getOption<std::string>("massOrderParameter");

    std::string reweight_string = reweights[0];
    for(unsigned int i=1; i<reweights.size(); i++) reweight_string+="*"+reweights[i];

    //what a declaration! basically a human-readable combination of parameters + mass ordering that map to a histogram
    //due to the lazy loading of RDataFrame these histograms aren't real histos until we do something like writing them
    //at which point it will do a single pass through the file and fill all histograms simultaneously
    std::map<std::tuple<std::string, std::string>, ROOT::RDF::RResultPtr<TH1D>> histograms_1d;
    std::map<std::tuple<std::string, std::string>, ROOT::RDF::RResultPtr<TH1D>> rc_histograms_1d;

    std::map<std::tuple<std::string, std::string, std::string>, ROOT::RDF::RResultPtr<TH2D>> histograms_2d;
    std::map<std::tuple<std::string, std::string, std::string>, ROOT::RDF::RResultPtr<TH2D>> rc_histograms_2d;

    ROOT::RDataFrame df(man->getOption<std::string>("posteriorTreeName").c_str(), infiles);

    //create a new column/temporary branch that stores the reactor weighting
    auto reweight = df.Define("total_weight", reweight_string);

    std::string burnin_filter = "step > "+std::to_string(burnin); //burnin cut
    std::map<std::string, std::string> mass_orderings = {{"MOMarg", ""},
                                                         {"IO", " && " + mass_order_parameter + "<0"},
                                                         {"NO", " && " + mass_order_parameter + ">0"}};

    for(auto [ordering, mass_filter] : mass_orderings){
        auto filter = burnin_filter + mass_filter;
        //basically similar to posteriors->Draw("theta_13", "step>XXXX && dm23>0")
        //creates a dataframe with just the steps we want for this specific configuration
        auto filtered_df = reweight.Filter(filter);

        //loop over parameters to set x axis parameter for 1d and 2d histograms
        for(auto parx : parameter_names){
            auto name = ordering+"_"+parx;
            auto title = ordering+" "+parx;
            auto tupleidx = std::make_tuple(ordering, parx);

            auto x = woRCbinning[parx]; //get the number of bins, bin lower and upper range
            //create essentially a promise to a TH1 that will be filled when we do something with this histogram later
            histograms_1d[tupleidx] = filtered_df.Histo1D({name.c_str(), title.c_str(), x[0].as<int>(), x[1].as<double>(), x[2].as<double>()}, parx.c_str());
            auto rcx = wRCbinning[parx]; //get the number of bins, bin lower and upper range
            rc_histograms_1d[tupleidx] = filtered_df.Histo1D({("rc_"+name).c_str(), ("RC "+title).c_str(), rcx[0].as<int>(), rcx[1].as<double>(), rcx[2].as<double>()}, parx.c_str(), "total_weight");

            //loop over parameters for the y axis of the 2d histograms
            for(auto pary : parameter_names){
                auto tupleidxy = std::make_tuple(ordering, parx, pary); //make a unique 'index' to the histogram
                name = ordering+"_"+parx+"_"+pary;
                title = ordering+" "+parx+" vs "+pary;

                auto y = woRCbinning[pary];
                histograms_2d[tupleidxy] = filtered_df.Histo2D({name.c_str(), title.c_str(), x[0].as<int>(), x[1].as<double>(), x[2].as<double>(),
                     y[0].as<int>(), y[1].as<double>(), y[2].as<double>()}, parx.c_str(), pary.c_str());

                auto rcy = wRCbinning[pary];
                rc_histograms_2d[tupleidxy] = filtered_df.Histo2D({("rc_"+name).c_str(), ("RC "+title).c_str(), rcx[0].as<int>(), rcx[1].as<double>(), rcx[2].as<double>(),
                     rcy[0].as<int>(), rcy[1].as<double>(), rcy[2].as<double>()}, parx.c_str(), pary.c_str(), "total_weight");
            }
        }
    }

    TFile *outf = new TFile(man->getOutputName().c_str(), "RECREATE");
    for(auto [key, histo] : histograms_1d) histo->Write();
    for(auto [key, histo] : rc_histograms_1d) histo->Write();
    for(auto [key, histo] : histograms_2d) histo->Write();
    for(auto [key, histo] : rc_histograms_2d) histo->Write();
    outf->Close();
}



int main(int argc, char **argv) {

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    const char* omp_env = std::getenv("OMP_NUM_THREADS");
    if (omp_env != nullptr) {
        int omp_threads = std::stoi(omp_env);
        ROOT::EnableImplicitMT(omp_threads);
    }
    else ROOT::EnableImplicitMT(16); //lets not go too crazy, limit to 16 threads if OMP_NUM_THREADS is not set
    
    SetMaCh3LoggerFormat();
    man = new M3::Plotting::PlottingManager();
    man->parseInputs(argc, argv);
    man->setExec("MakePosteriorHistograms");
    MakePosteriorHistograms();
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> runtime = std::chrono::duration_cast<std::chrono::nanoseconds>(end-start);
    MACH3LOG_INFO("Execution took {}s", runtime.count());

}
