#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <cmath>

#include "yaml-cpp/yaml.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLine.h"
#include "TMacro.h"
#include "TSystem.h"
#include "TTree.h"

/// @author David Riley

void UmbrellaPlotting() {
    // ============================================
    // Configuration: Number of bins for each variable
    // ============================================
    int nBins_dcp = 200;           // Number of bins for delta_cp
    int nBins_sin2th23 = 100;      // Number of bins for sin^2(theta_23)
    int nBins_sin2th13 = 100;
    int nBins_sin2th12 = 100;
    int nBins_delm23 = 3000;        // Number of bins for delta_m23^2
    int nBins_delm12 = 100;

    // Unified plotting ranges requested for all sin2th23 and delm23 views.
    double sin2th23_min = 0.45;
    double sin2th23_max = 0.57;
    double sin2th13_min = 0.015;
    double sin2th13_max = 0.03;
    double sin2th12_min = 0.24;
    double sin2th12_max = 0.36;
    double delm23_no_min = 0.00245;
    double delm23_no_max = 0.0026;
    double delm23_io_min = -delm23_no_max;
    double delm23_io_max = -delm23_no_min;
    double delm12_min = 6.5e-5;
    double delm12_max = 8.5e-5;

    auto normalizeHist = [&](TH1 *h) {
        if (!h) return;
        double integral = h->Integral();
        if (integral > 0.0) h->Scale(1.0 / integral);
    };

    auto printHistStats = [&](const std::string &label, TH1 *h) {
        if (!h) {
            std::cout << "[DEBUG] " << label << " : histogram is null" << std::endl;
            return;
        }
        std::cout << std::fixed << std::setprecision(6)
                  << "[DEBUG] " << label
                  << " entries=" << h->GetEntries()
                  << " effEntries=" << h->GetEffectiveEntries()
                  << " integral=" << h->Integral()
                  << " mean=" << h->GetMean()
                  << " rms=" << h->GetRMS()
                  << std::endl;
    };
    
    // ============================================
    // Open the ROOT file
    // ============================================
    //TFile *f = TFile::Open("hk_umbrella_stepscale_ON_width_0.05_mean_0.725_flipping_0_verbose_test.root");  // Update with your actual filename
    // TFile *f = TFile::Open("/home/ppe/d/driley/Projects/hk_umbrella/mach3-hk/build/output/hk_umbrella_stepscale_ON_width_0.01_mean_0.725_flipping_0_verbose_test_FLIPPING_DELM23_FIX.root");  // Update with your actual filename
    // TFile *f = TFile::Open("/home/ppe/d/driley/Projects/hk_umbrella/mach3-hk/build/output/hk_umbrella_stepscale_ON_width_0.01_mean_0.725_flipping_0_verbose_test_FLIPPING_DELM23_FIX_long.root"); 
    // TFile *f = TFile::Open("/home/ppe/d/driley/Projects/umbrella/utils_mcm/utils_umbrella_multicanonical/window_solver/umbrella_output_fullwidth.root");
    //TFile *f = TFile::Open("/project/6002456/driley/T2K_outputs/umbrella_testing/AsimovA22_wider_vonMises/window_solve_output/A22_wider_vonMises.root");
    
    // current inputs
    //TFile *f = TFile::Open("/home/driley/Projects-T2K-MaCh3/umbrella_t2k/mcm_utils/utils_umbrella_multicanonical/window_solver/A22_dense_windows_subset_test.root");
    //std::string outputDir = "outputs/dense_test/";
    //TFile *f = TFile::Open("/home/driley/Projects-T2K-MaCh3/umbrella_t2k/mcm_utils/utils_umbrella_multicanonical/window_solver/outputs/slurm_31456021/A22_dense_windows_subset_test.root");
    //std::string outputDir = "outputs/t2k_denser_windows_vmwideish/";
   
    //TFile *f = TFile::Open("/home/driley/Projects-T2K-MaCh3/umbrella_t2k/mcm_utils/utils_umbrella_multicanonical/window_solver/outputs/slurm_32871924/A22_gaussian_test.root");
    //std::string outputDir = "outputs/t2k_gaussian_but_accidentally_vonMises";
    
    //TFile *f = TFile::Open("/home/driley/Projects-T2K-MaCh3/umbrella_t2k/mcm_utils/utils_umbrella_multicanonical/window_solver/pihalf_fine_umbrella.root");
    //std::string outputDir = "outputs/hk_fine/";

    // Raw window file used to read MaCh3_Config YAML (contains Asimov metadata).
    //std::string mach3ConfigSourceFile = "/project/6002456/driley/HyperK_output_archive/umbrella_testing/dcp_mpihalf_gen_gaus_020626/burnin_window_0.root";
    // Combined plotting input.
    //TFile *f = TFile::Open("/project/6002456/driley/HyperK_output_archive/umbrella_testing/dcp_mpihalf_gen_gaus_020626/gen_gaussian_mpihalf.root");
    //std::string outputDir = "outputs/hk_gengaus_dcpmpihalf/";
    
    //std::string mach3ConfigSourceFile = "/project/6002456/driley/HyperK_output_archive/umbrella_testing/dcpzero_gen_gaus_020626/burnin_window_0.root";
    //TFile *f = TFile::Open("/project/6002456/driley/HyperK_output_archive/umbrella_testing/dcpzero_gen_gaus_020626/gen_gaussian_dcpzerp.root");
    //std::string outputDir = "outputs/hk_gengaus_dcpzero/";
    
    std::string mach3ConfigSourceFile = "/home/driley/projects/def-blairt2k/driley/HyperK_output_archive/umbrella_testing/hk_gen_gaus_dense_2205/archive_2905_sinth23flipping/burnin_window_1.root";
    TFile *f = TFile::Open("/home/driley/projects/def-blairt2k/driley/HyperK_output_archive/umbrella_testing/hk_gen_gaus_dense_2205/archive_2905_sinth23flipping/gen_gaussian_dense_flipon.root");
    std::string outputDir = "outputs/hk_gengaus_new_dense_flippingon/";

    //std::string mach3ConfigSourceFile = "/project/6002456/driley/HyperK_output_archive/umbrella_testing/hk_p_pihalf_IO/output_accidentally_NO/burnin_window_0.root";
    //TFile *f = TFile::Open("/project/6002456/driley/HyperK_output_archive/umbrella_testing/hk_p_pihalf_NO/accidental_NO_results/gen_gaussian_dcpzero_finer_NO_only.root");
    //std::string outputDir = "outputs/hk_p_pihalf_NO_only/";
    
    //std::string mach3ConfigSourceFile = "/project/6002456/driley/HyperK_output_archive/umbrella_testing/dcpzero_gen_gaus_finer_020626/burnin_window_0.root";
    //TFile *f = TFile::Open("/project/6002456/driley/HyperK_output_archive/umbrella_testing/dcpzero_gen_gaus_finer_020626/gen_gaussian_dcpzero_finer.root");
    //std::string outputDir = "outputs/hk_gengaus_dcpzero_finer/";

    // T2K regerence posterior input
    //std::string referenceFile = "/home/driley/projects/def-blairt2k/driley/T2K_outputs/reduce_AsimovA22Fit_070723_FinalChain.root";
    //std::string referenceTreeName = "osc_posteriors";
    bool use_reference_chain = false;
    bool use_reference_LLH_scan = true;
    
    // load in reference chain
    std::string referenceFile = "";
    std::string referenceTreeName = "";

    if(use_reference_chain) {
        referenceFile = "/home/driley/projects/def-blairt2k/driley/T2K_outputs/umbrella_testing/fixing_params/allfixed/all_oscParams_fixed_hadded.root";
        referenceTreeName = "posteriors";
    }

    // setup reference LLH scan input (works best when fit is only over the required variables)
    std::string referenceLLHScanFile = "";
    std::string referenceLLHScanDcpHistName = "";
    std::string referenceLLHScanDcp_NO_histName = "";
    std::string referenceLLHScanDcp_IO_histName = "";

    if(use_reference_LLH_scan) {
        referenceLLHScanFile = "/home/driley/Projects-T2K-MaCh3/umbrella_t2k/mcm_utils/target_1D_histos.root";
        referenceLLHScanDcpHistName = "dcp";
        referenceLLHScanDcp_NO_histName = "dcp_NO";
        referenceLLHScanDcp_IO_histName = "dcp_IO";
    }

    // begin loading in fit file
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open file!" << std::endl;
        return;
    }
    
    // Get the tree
    TTree *tree = (TTree*)f->Get("posteriors");  // Update with your tree name
    if (!tree) {
        std::cerr << "Error: Cannot find tree!" << std::endl;
        return;
    }

    // Read Asimov metadata from MaCh3_Config YAML in the raw chain file.
    // this is for plotting of asimov point 
    bool hasMaCh3Config = false;
    bool isAsimovChain = false;
    std::vector<double> asimovOscillationParameters;

    double asimovSin2th12 = std::numeric_limits<double>::quiet_NaN();
    double asimovSin2th23 = std::numeric_limits<double>::quiet_NaN();
    double asimovSin2th13 = std::numeric_limits<double>::quiet_NaN();
    double asimovDelm12 = std::numeric_limits<double>::quiet_NaN();
    double asimovDelm23 = std::numeric_limits<double>::quiet_NaN();
    double asimovDcp = std::numeric_limits<double>::quiet_NaN();
    double asimovBaseline = std::numeric_limits<double>::quiet_NaN();
    double asimovElectronDensity = std::numeric_limits<double>::quiet_NaN();

    TFile *fConfig = TFile::Open(mach3ConfigSourceFile.c_str());
    if (!fConfig || fConfig->IsZombie()) {
        std::cerr << "Error: Cannot open MaCh3 config source file: " << mach3ConfigSourceFile << std::endl;
        return;
    } else {
        TMacro *mach3ConfigMacro = dynamic_cast<TMacro *>(fConfig->Get("MaCh3_Config"));
        if (!mach3ConfigMacro) {
            std::cerr << "Error: MaCh3_Config not found in " << mach3ConfigSourceFile << std::endl;
            fConfig->Close();
            return;
        } else {
            hasMaCh3Config = true;
            const std::string tmpYamlPath = std::string(gSystem->TempDirectory()) + "/umbrella_plotting_MaCh3_Config.yaml";
            mach3ConfigMacro->SaveSource(tmpYamlPath.c_str());

            try {
                YAML::Node cfg = YAML::LoadFile(tmpYamlPath);
                YAML::Node general = cfg["General"];
                isAsimovChain = general["Asimov"].as<bool>();
                YAML::Node osc = general["OscillationParameters"];
                for (std::size_t i = 0; i < osc.size(); ++i) {
                    asimovOscillationParameters.push_back(osc[i].as<double>());
                }
            } catch (const std::exception &e) {
                std::cerr << "Error: Failed to parse MaCh3_Config YAML with yaml-cpp: " << e.what() << std::endl;
                fConfig->Close();
                return;
            }
        }
        fConfig->Close();
    }

    asimovSin2th12 = asimovOscillationParameters[0];
    asimovSin2th23 = asimovOscillationParameters[1];
    asimovSin2th13 = asimovOscillationParameters[2];
    asimovDelm12 = asimovOscillationParameters[3];
    asimovDelm23 = asimovOscillationParameters[4];
    asimovDcp = asimovOscillationParameters[5];
    asimovBaseline = asimovOscillationParameters[6];
    asimovElectronDensity = asimovOscillationParameters[7];

    std::cout << "[DEBUG] hasMaCh3Config=" << hasMaCh3Config
              << " isAsimovChain=" << isAsimovChain
              << " oscParamCount=" << asimovOscillationParameters.size() << std::endl;
    if (isAsimovChain && !asimovOscillationParameters.empty()) {
        std::cout << std::fixed << std::setprecision(6)
                  << "[DEBUG] Asimov named values: "
                  << "sin2th12=" << asimovSin2th12
                  << " sin2th23=" << asimovSin2th23
                  << " sin2th13=" << asimovSin2th13
                  << " delm12=" << asimovDelm12
                  << " delm23=" << asimovDelm23
                  << " dcp(rad)=" << asimovDcp
                  << " baseline=" << asimovBaseline
                  << " electronDensity=" << asimovElectronDensity
                  << std::endl;
    }

    auto calculateJarlskog = [&](double sin2th12, double sin2th23, double sin2th13, double dcp) {
        if (!std::isfinite(sin2th12) || !std::isfinite(sin2th23) || !std::isfinite(sin2th13) || !std::isfinite(dcp)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return std::sqrt(std::max(0.0, sin2th12) * std::max(0.0, 1.0 - sin2th12) *
                         std::max(0.0, sin2th23) * std::max(0.0, 1.0 - sin2th23)) *
               std::max(0.0, 1.0 - sin2th13) * std::sqrt(std::max(0.0, sin2th13)) *
               std::sin(dcp);
    };

    const double asimovJarlskog = calculateJarlskog(asimovSin2th12, asimovSin2th23, asimovSin2th13, asimovDcp);

    auto drawAsimovLineWithLegend = [&](double xValue, TLegend *existingLegend = nullptr) {
        if (!isAsimovChain || !std::isfinite(xValue) || !gPad) return;
        gPad->Update();
        const double yMin = gPad->GetUymin();
        const double yMax = gPad->GetUymax();
        TLine *asimovLine = new TLine(xValue, yMin, xValue, yMax);
        asimovLine->SetLineColor(kBlack);
        asimovLine->SetLineStyle(3);
        asimovLine->SetLineWidth(2);
        asimovLine->Draw("SAME");

        if (existingLegend) {
            existingLegend->AddEntry(asimovLine, "Asimov point", "l");
            existingLegend->Draw();
        } else {
            TLegend *asimovLegend = new TLegend(0.14, 0.83, 0.34, 0.91);
            asimovLegend->SetBorderSize(0);
            asimovLegend->SetFillStyle(0);
            asimovLegend->AddEntry(asimovLine, "Asimov point", "l");
            asimovLegend->Draw();
        }
    };

    auto drawAsimovMarkerWithLegend = [&](double xValue, double yValue, TLegend *existingLegend = nullptr) {
        if (!isAsimovChain || !std::isfinite(xValue) || !std::isfinite(yValue) || !gPad) return;
        gPad->Update();
        // small cross sized relative to the axis ranges
        const double xMin = gPad->GetUxmin();
        const double xMax = gPad->GetUxmax();
        const double yMin = gPad->GetUymin();
        const double yMax = gPad->GetUymax();
        const double dx = 0.012 * (xMax - xMin);
        const double dy = 0.012 * (yMax - yMin);

        TLine *l1 = new TLine(xValue - dx, yValue - dy, xValue + dx, yValue + dy);
        l1->SetLineColor(kRed);
        l1->SetLineWidth(2);
        l1->Draw("SAME");

        TLine *l2 = new TLine(xValue - dx, yValue + dy, xValue + dx, yValue - dy);
        l2->SetLineColor(kRed);
        l2->SetLineWidth(2);
        l2->Draw("SAME");

        if (existingLegend) {
            existingLegend->AddEntry(l1, "Asimov point", "l");
            existingLegend->Draw();
        } else {
            TLegend *lm = new TLegend(0.14, 0.83, 0.34, 0.91);
            lm->SetBorderSize(0);
            lm->SetFillStyle(0);
            lm->AddEntry(l1, "Asimov point", "l");
            lm->Draw();
        }
    };


    // Build reference histograms from osc_posteriors for direct overlay.
    bool hasReferenceDcp = false;
    bool hasReferenceSin2Th23 = false;
    bool hasReferenceDelm23 = false;
    TH1D *hRefDcpShape = nullptr;
    TH1D *hRefSin2Th23Shape = nullptr;
    TH1D *hRefDelm23Shape = nullptr;
    TFile *fRef = TFile::Open(referenceFile.c_str());
    if (!fRef || fRef->IsZombie()) {
        std::cerr << "Warning: Cannot open reference file: " << referenceFile << std::endl;
    } else {
        TTree *refTree = (TTree*)fRef->Get(referenceTreeName.c_str());
        if (!refTree) {
            std::cerr << "Warning: Cannot find reference tree: " << referenceTreeName << std::endl;
        } else if (!refTree->GetBranch("dcp") && !refTree->GetBranch("delta_cp")) {
            std::cerr << "Warning: dcp branch not found in reference tree." << std::endl;
        } else if (refTree->GetBranch("delta_cp")) {
            hRefDcpShape = new TH1D("hRefDcpShape", "", nBins_dcp, -3.1415, 3.1415);
            refTree->Draw("delta_cp>>hRefDcpShape", "", "goff");
            hRefDcpShape->SetDirectory(0);
            hasReferenceDcp = true;
            std::cout << "Loaded reference dcp histogram from " << referenceTreeName << std::endl;
	    } else {
            hRefDcpShape = new TH1D("hRefDcpShape", "", nBins_dcp, -3.1415, 3.1415);
            refTree->Draw("dcp>>hRefDcpShape", "", "goff");
            hRefDcpShape->SetDirectory(0);
            hasReferenceDcp = true;
            std::cout << "Loaded reference dcp histogram from " << referenceTreeName << std::endl;
        }

        if (!refTree->GetBranch("theta23") && !refTree->GetBranch("sin2th_23")) {
            std::cerr << "Warning: theta23 branch not found in reference tree." << std::endl;
        } else if (refTree->GetBranch("sin2th_23")) {
            hRefSin2Th23Shape = new TH1D("hRefSin2Th23Shape", "", nBins_sin2th23, sin2th23_min, sin2th23_max);
            refTree->Draw("sin2th_23>>hRefSin2Th23Shape", "", "goff");
            hRefSin2Th23Shape->SetDirectory(0);
            hasReferenceSin2Th23 = true;
            std::cout << "Loaded reference sin^{2}(theta23) histogram from " << referenceTreeName << std::endl;
	    
    	} else {
            hRefSin2Th23Shape = new TH1D("hRefSin2Th23Shape", "", nBins_sin2th23, sin2th23_min, sin2th23_max);
            refTree->Draw("theta23>>hRefSin2Th23Shape", "", "goff");
            hRefSin2Th23Shape->SetDirectory(0);
            hasReferenceSin2Th23 = true;
            std::cout << "Loaded reference sin^{2}(theta23) histogram from " << referenceTreeName << std::endl;
        }

        if (!refTree->GetBranch("dm23") && !refTree->GetBranch("delm2_23")) {
            std::cerr << "Warning: dm23 branch not found in reference tree." << std::endl;
        } else if (refTree->GetBranch("delm2_23")) {
            hRefDelm23Shape = new TH1D("hRefDelm23Shape", "", nBins_delm23, delm23_io_min, delm23_no_max);
            refTree->Draw("delm2_23>>hRefDelm23Shape", "", "goff");
            TH1D *hRefDelm23IO_dbg = new TH1D("hRefDelm23IO_dbg", "", nBins_delm23, delm23_io_min, delm23_io_max);
            TH1D *hRefDelm23NO_dbg = new TH1D("hRefDelm23NO_dbg", "", nBins_delm23, delm23_no_min, delm23_no_max);
            refTree->Draw("delm2_23>>hRefDelm23IO_dbg", "delm2_23 < 0", "goff");
            refTree->Draw("delm2_23>>hRefDelm23NO_dbg", "delm2_23 > 0", "goff");
            printHistStats("Reference delm23 IO (raw)", hRefDelm23IO_dbg);
            printHistStats("Reference delm23 NO (raw)", hRefDelm23NO_dbg);
            hRefDelm23Shape->SetDirectory(0);
            normalizeHist(hRefDelm23Shape);
            hasReferenceDelm23 = true;
            std::cout << "Loaded reference delm23 histogram from " << referenceTreeName << std::endl;
	    
    	}  else {
            hRefDelm23Shape = new TH1D("hRefDelm23Shape", "", nBins_delm23, delm23_io_min, delm23_no_max);
            refTree->Draw("dm23>>hRefDelm23Shape", "", "goff");
            TH1D *hRefDm23IO_dbg = new TH1D("hRefDm23IO_dbg", "", nBins_delm23, delm23_io_min, delm23_io_max);
            TH1D *hRefDm23NO_dbg = new TH1D("hRefDm23NO_dbg", "", nBins_delm23, delm23_no_min, delm23_no_max);
            refTree->Draw("dm23>>hRefDm23IO_dbg", "dm23 < 0", "goff");
            refTree->Draw("dm23>>hRefDm23NO_dbg", "dm23 > 0", "goff");
            printHistStats("Reference dm23 IO (raw)", hRefDm23IO_dbg);
            printHistStats("Reference dm23 NO (raw)", hRefDm23NO_dbg);
            hRefDelm23Shape->SetDirectory(0);
            normalizeHist(hRefDelm23Shape);
            hasReferenceDelm23 = true;
            std::cout << "Loaded reference dm23 histograms from " << referenceTreeName << std::endl;
        }
        fRef->Close();
    }
    
    // --------------------------------------------------
    // Load LLH-scan reference histograms (linear likelihoods)
    // --------------------------------------------------
    bool hasReferenceLLH_Dcp = false;
    bool hasReferenceLLH_Dcp_NO = false;
    bool hasReferenceLLH_Dcp_IO = false;
    TH1D *hRefLLH_Dcp = nullptr;
    TH1D *hRefLLH_Dcp_NO = nullptr;
    TH1D *hRefLLH_Dcp_IO = nullptr;

    if (use_reference_LLH_scan && !referenceLLHScanFile.empty()) {
        TFile *fLLH = TFile::Open(referenceLLHScanFile.c_str());
        if (!fLLH || fLLH->IsZombie()) {
            std::cerr << "Warning: Cannot open reference LLH-scan file: " << referenceLLHScanFile << std::endl;
        } else {
            // Generic getter that leaves a cloned, detached, normalized histogram
            auto loadAndNormalize1D = [&](const std::string &name, TH1D *&outHist, bool &outFlag) {
                TObject *obj = fLLH->Get(name.c_str());
                if (!obj) {
                    std::cerr << "Warning: LLH-scan histogram '" << name << "' not found in " << referenceLLHScanFile << std::endl;
                    return;
                }
                TH1 *htmp = dynamic_cast<TH1*>(obj);
                if (!htmp) {
                    std::cerr << "Warning: Object '" << name << "' is not a histogram." << std::endl;
                    return;
                }
                outHist = (TH1D*)htmp->Clone((std::string("hRefLLH_") + name + "_clone").c_str());
                outHist->SetDirectory(0);
                // LLH scans are likelihoods (not log-likelihood) so area-normalize for overlay
                double I = outHist->Integral();
                if (I > 0) outHist->Scale(1.0 / I);
                outFlag = true;
                printHistStats(std::string("Reference LLH ") + name, outHist);
            };

            if (!referenceLLHScanDcpHistName.empty()) loadAndNormalize1D(referenceLLHScanDcpHistName, hRefLLH_Dcp, hasReferenceLLH_Dcp);
            if (!referenceLLHScanDcp_NO_histName.empty()) loadAndNormalize1D(referenceLLHScanDcp_NO_histName, hRefLLH_Dcp_NO, hasReferenceLLH_Dcp_NO);
            if (!referenceLLHScanDcp_IO_histName.empty()) loadAndNormalize1D(referenceLLHScanDcp_IO_histName, hRefLLH_Dcp_IO, hasReferenceLLH_Dcp_IO);

            fLLH->Close();
        }
    }
    
    bool hasUmbrellaWeight = false;
    // check if the umbrella weight branch exists
    if (!tree->GetBranch("umbrella_weight")) {
        std::cerr << "Warning: umbrella_weight branch not found, plotting as unweighted umbrella" << std::endl;
    } else {
        hasUmbrellaWeight = true;
    }



    std::cout << "[DEBUG] Chain entries total=" << tree->GetEntries() << std::endl;
    TH1D *hChainDelm23IO_dbg = new TH1D("hChainDelm23IO_dbg", "", nBins_delm23, delm23_io_min, delm23_io_max);
    TH1D *hChainDelm23NO_dbg = new TH1D("hChainDelm23NO_dbg", "", nBins_delm23, delm23_no_min, delm23_no_max);
    tree->Draw("delm2_23>>hChainDelm23IO_dbg", "delm2_23 < 0", "goff");
    tree->Draw("delm2_23>>hChainDelm23NO_dbg", "delm2_23 > 0", "goff");
    printHistStats("Chain delm23 IO (raw)", hChainDelm23IO_dbg);
    printHistStats("Chain delm23 NO (raw)", hChainDelm23NO_dbg);

    if (hasUmbrellaWeight) {
        TH1D *hChainDelm23IO_w_dbg = new TH1D("hChainDelm23IO_w_dbg", "", nBins_delm23, delm23_io_min, delm23_io_max);
        TH1D *hChainDelm23NO_w_dbg = new TH1D("hChainDelm23NO_w_dbg", "", nBins_delm23, delm23_no_min, delm23_no_max);
        tree->Draw("delm2_23>>hChainDelm23IO_w_dbg", "(delm2_23 < 0) * umbrella_weight", "goff");
        tree->Draw("delm2_23>>hChainDelm23NO_w_dbg", "(delm2_23 > 0) * umbrella_weight", "goff");
        printHistStats("Chain delm23 IO (weighted)", hChainDelm23IO_w_dbg);
        printHistStats("Chain delm23 NO (weighted)", hChainDelm23NO_w_dbg);
    }

    // Create canvas
    TCanvas *c1 = new TCanvas("c1", "Delta CP vs Delta m23", 900, 700);
    c1->SetRightMargin(0.15);
    // Create 2D histogram: delta m23 vs delta CP
    TH2F *h2d = new TH2F("h2d", "#Delta m_{23}^{2} vs #delta_{CP};#delta_{CP};#Delta m_{23}^{2} (eV^{2})", 
                         nBins_dcp, -3.14159, 3.14159, 
                         nBins_delm23, delm23_io_min, delm23_no_max);
    // Draw from tree 
    tree->Draw("delm2_23:delta_cp>>h2d", "", "COLZ");
    normalizeHist(h2d);
    h2d->Draw("COLZ");
    // Styling
    h2d->SetStats(0);  // Remove stats box
    h2d->GetXaxis()->SetTitle("#delta_{CP}");
    h2d->GetYaxis()->SetTitle("#Delta m_{23}^{2} (eV^{2})");
    h2d->GetZaxis()->SetTitle("Normalized Entries");
    //gPad->SetLogz();  // Optional: log scale on z-axis for better visibility
    c1->Update();
    // draw Asimov marker and Save plot
    drawAsimovMarkerWithLegend(asimovDcp, asimovDelm23, nullptr);
    c1->SaveAs((outputDir + "delm23_vs_dcp.png").c_str());

    if (hasUmbrellaWeight){
        tree->Draw("delm2_23:delta_cp>>h2d", "umbrella_weight", "COLZ");
        normalizeHist(h2d);
        h2d->Draw("COLZ");
        c1->Update();
        drawAsimovMarkerWithLegend(asimovDcp, asimovDelm23, nullptr);
        c1->SaveAs((outputDir + "delm23_vs_dcp_weighted.png").c_str());
    }


    // Create side-by-side canvas for IO (left) and NO (right)
    TCanvas *c2 = new TCanvas("c2", "Delta m23 IO and NO", 1800, 700);
    c2->Divide(2, 1);
    
    // Left pad: IO
    c2->cd(1);
    gPad->SetRightMargin(0.15);
    // Create 2D histogram: dcp vs delm23 
    TH2F *h2d_IO = new TH2F("h2d_IO", "Inverted Ordering;#Delta m_{23}^{2} (eV^{2});#delta_{CP}", 
                             nBins_delm23, delm23_io_min, delm23_io_max,  // X-axis: delm23 (IO range)
                             nBins_dcp, -3.1415, 3.1415);               // Y-axis: dcp
    // Draw from tree 
    tree->Draw("delta_cp:delm2_23>>h2d_IO", "delm2_23 < 0", "COLZ");
    normalizeHist(h2d_IO);
    h2d_IO->Draw("COLZ");
    // Styling
    h2d_IO->SetStats(0);
    h2d_IO->GetXaxis()->SetTitle("#Delta m_{23}^{2} (eV^{2})");
    h2d_IO->GetYaxis()->SetTitle("#delta_{CP} / #pi");
    h2d_IO->GetZaxis()->SetTitle("Normalized Entries");
    //gPad->SetLogz();
    
    // Right pad: NO
    c2->cd(2);
    gPad->SetRightMargin(0.15);
    // Create 2D histogram: dcp vs delm23 
    TH2F *h2d_NO = new TH2F("h2d_NO", "Normal Ordering;#Delta m_{23}^{2} (eV^{2});#delta_{CP}", 
                             nBins_delm23, delm23_no_min, delm23_no_max,   // X-axis: delm23 (NO range)
                             nBins_dcp, -3.1415, 3.1415);               // Y-axis: dcp
    // Draw from tree 
    tree->Draw("delta_cp:delm2_23>>h2d_NO", "delm2_23 > 0", "COLZ");
    normalizeHist(h2d_NO);
    h2d_NO->Draw("COLZ");
    // Styling
    h2d_NO->SetStats(0);
    h2d_NO->GetXaxis()->SetTitle("#Delta m_{23}^{2} (eV^{2})");
    h2d_NO->GetYaxis()->SetTitle("#delta_{CP} / #pi");
    h2d_NO->GetZaxis()->SetTitle("Normalized Entries");
    //gPad->SetLogz();
    
    // Update canvas
    // draw Asimov markers per pad, update and save
    c2->cd(1);
    drawAsimovMarkerWithLegend(asimovDelm23, asimovDcp, nullptr);
    c2->cd(2);
    drawAsimovMarkerWithLegend(asimovDelm23, asimovDcp, nullptr);
    c2->Update();
    c2->SaveAs((outputDir + "delm23_vs_dcp_IO_NO.png").c_str());

    if (hasUmbrellaWeight){
        // Left pad: IO weighted
        c2->cd(1);
        tree->Draw("delta_cp:delm2_23>>h2d_IO", "(delm2_23 < 0) * umbrella_weight", "COLZ");
        normalizeHist(h2d_IO);
        h2d_IO->Draw("COLZ");
        
        // Right pad: NO weighted
        c2->cd(2);
        tree->Draw("delta_cp:delm2_23>>h2d_NO", "(delm2_23 > 0) * umbrella_weight", "COLZ");
        normalizeHist(h2d_NO);
        h2d_NO->Draw("COLZ");
        
        // Update canvas
        // draw Asimov markers per pad, update and save weighted
        c2->cd(1);
        drawAsimovMarkerWithLegend(asimovDelm23, asimovDcp, nullptr);
        c2->cd(2);
        drawAsimovMarkerWithLegend(asimovDelm23, asimovDcp, nullptr);
        c2->Update();
        c2->SaveAs((outputDir + "delm23_vs_dcp_IO_NO_weighted.png").c_str());
    }

    // Lambda to calculate HPD contour level for given confidence level (0-1)
    auto calculateHPDLevel = [](TH2F *hist, double confidenceLevel) -> double {
        std::vector<double> binContents;
        for (int i = 1; i <= hist->GetNbinsX(); ++i) {
            for (int j = 1; j <= hist->GetNbinsY(); ++j) {
                double content = hist->GetBinContent(i, j);
                if (content > 0) {
                    binContents.push_back(content);
                }
            }
        }
        
        if (binContents.empty()) return 0.0;
        
        std::sort(binContents.rbegin(), binContents.rend()); // Sort in descending order
        
        double totalSum = 0;
        for (double val : binContents) totalSum += val;
        
        double targetSum = confidenceLevel * totalSum;
        double runningSum = 0;
        
        for (double val : binContents) {
            runningSum += val;
            if (runningSum >= targetSum) {
                return val;
            }
        }
        return binContents.back();
    };

    // Lambda to create wrapped histogram for circular delta_cp coordinate
    auto createWrappedHistogram = [](TH2F *originalHist, const std::string &name) -> TH2F* {
        int nBinsX = originalHist->GetNbinsX();
        int nBinsY = originalHist->GetNbinsY();
        
        double xMin = originalHist->GetXaxis()->GetBinLowEdge(1);
        double xMax = originalHist->GetXaxis()->GetBinUpEdge(nBinsX);
        double yMin = originalHist->GetYaxis()->GetBinLowEdge(1);
        double yMax = originalHist->GetYaxis()->GetBinUpEdge(nBinsY);
        double yRange = yMax - yMin;
        
        // Extended histogram spanning 3 periods for proper boundary wrapping
        TH2F *wrapped = new TH2F(name.c_str(), "",
                                 nBinsX, xMin, xMax,
                                 nBinsY * 3, yMin - yRange, yMax + yRange);
        wrapped->SetDirectory(0);
        
        // Fill with periodically wrapped data
        for (int i = 1; i <= nBinsX; ++i) {
            for (int j = 1; j <= nBinsY; ++j) {
                double content = originalHist->GetBinContent(i, j);
                // Center period
                wrapped->SetBinContent(i, j + nBinsY, content);
                // Lower wrapped period
                wrapped->SetBinContent(i, j, content);
                // Upper wrapped period
                wrapped->SetBinContent(i, j + 2*nBinsY, content);
            }
        }
        
        return wrapped;
    };

    // Helper to compute one-sided bounds (do not draw) — integrates from side with more posterior
    auto computeOneSidedBoundsValues = [&](TH1 *hist, double targets[3], double outBounds[3], bool &integrateLeftToRight) {
        for (int i = 0; i < 3; ++i) outBounds[i] = std::numeric_limits<double>::quiet_NaN();
        integrateLeftToRight = true;
        if (!hist) return;
        double total = hist->Integral();
        if (total <= 0) return;
        int nBins = hist->GetNbinsX();
        double leftSum = 0.0;
        double rightSum = 0.0;
        for (int b = 1; b <= nBins; ++b) {
            double center = hist->GetBinCenter(b);
            double c = hist->GetBinContent(b);
            if (center < 0) leftSum += c; else rightSum += c;
        }
        integrateLeftToRight = (leftSum >= rightSum);

        for (int t = 0; t < 3; ++t) {
            double target = targets[t];
            double cum = 0.0;
            if (integrateLeftToRight) {
                for (int b = 1; b <= nBins; ++b) {
                    cum += hist->GetBinContent(b);
                    if (cum / total >= target) {
                        outBounds[t] = hist->GetBinLowEdge(b+1);
                        break;
                    }
                }
                if (!std::isfinite(outBounds[t])) outBounds[t] = hist->GetXaxis()->GetXmax();
            } else {
                for (int b = nBins; b >= 1; --b) {
                    cum += hist->GetBinContent(b);
                    if (cum / total >= target) {
                        outBounds[t] = hist->GetBinLowEdge(b);
                        break;
                    }
                }
                if (!std::isfinite(outBounds[t])) outBounds[t] = hist->GetXaxis()->GetXmin();
            }
        }
    };

    auto drawAbsDelm23Overlay = [&](const std::string &suffix,
                                    const std::string &ioSelection,
                                    const std::string &noSelection) {
        const int nBinsAbsDelm23 = 80;
        const int nBinsAbsDcp = 60;

        // Set reversed greyscale palette (grey to black) for this plot only
        gStyle->SetPalette(kGreyScale);
        TColor::InvertPalette();

        TCanvas *cAbs = new TCanvas(("c_abs_" + suffix).c_str(), "Delta m23 magnitude overlay", 900, 700);
        cAbs->SetRightMargin(0.15);
        cAbs->SetLeftMargin(0.12);
        cAbs->SetBottomMargin(0.12);

        TH2F *hAbsFrame = new TH2F(("hAbsFrame_" + suffix).c_str(),
                                   "|#Delta m_{23}^{2}| vs #delta_{CP};|#Delta m_{23}^{2}| (eV^{2});#delta_{CP} / #pi",
                                              nBinsAbsDelm23, delm23_no_min, delm23_no_max,
                                              nBinsAbsDcp, -3.1415, 3.1415);
        hAbsFrame->SetStats(0);
        hAbsFrame->Draw();

        TH2F *hAbsIO = new TH2F(("hAbsIO_" + suffix).c_str(), "",
                                          nBinsAbsDelm23, delm23_no_min, delm23_no_max,
                                          nBinsAbsDcp, -3.1415, 3.1415);
        TH2F *hAbsNO = new TH2F(("hAbsNO_" + suffix).c_str(), "",
                                          nBinsAbsDelm23, delm23_no_min, delm23_no_max,
                                          nBinsAbsDcp, -3.1415, 3.1415);

        std::string ioDrawExpr = std::string("delta_cp:TMath::Abs(delm2_23)>>") + hAbsIO->GetName();
        std::string noDrawExpr = std::string("delta_cp:TMath::Abs(delm2_23)>>") + hAbsNO->GetName();
        tree->Draw(ioDrawExpr.c_str(), ioSelection.c_str(), "goff");
        tree->Draw(noDrawExpr.c_str(), noSelection.c_str(), "goff");

        hAbsIO->SetDirectory(0);
        hAbsNO->SetDirectory(0);
        printHistStats("Abs delm23 IO", hAbsIO);
        printHistStats("Abs delm23 NO", hAbsNO);

        // Calculate HPD contour levels for 1-sigma (68.3%)
        double hpdLevel1SigmaNO = calculateHPDLevel(hAbsNO, 0.683);
        
        double hpdLevel1SigmaIO = calculateHPDLevel(hAbsIO, 0.683);

        // Draw density first
        hAbsNO->SetLineColor(kAzure + 2);
        hAbsNO->SetLineWidth(1);
        hAbsNO->Draw("COLZ");

        hAbsIO->SetLineColor(kOrange + 7);
        hAbsIO->SetLineWidth(1);
        hAbsIO->Draw("COLZ SAME");

        // Smooth the fine histograms for contours
        TH2F *hNO_contour = (TH2F *)hAbsNO->Clone("hNO_contour");
        TH2F *hIO_contour = (TH2F *)hAbsIO->Clone("hIO_contour");
        hNO_contour->SetDirectory(0);
        hIO_contour->SetDirectory(0);

        hNO_contour->Smooth(1, "k5a");
        hIO_contour->Smooth(1, "k5a");

        double noLevels[1] = {hpdLevel1SigmaNO};
        hNO_contour->SetContour(1, noLevels);
        hNO_contour->SetLineColor(kAzure + 2);
        hNO_contour->SetLineWidth(3);
        hNO_contour->Draw("CONT3 SAME");

        double ioLevels[1] = {hpdLevel1SigmaIO};
        hIO_contour->SetContour(1, ioLevels);
        hIO_contour->SetLineColor(kOrange + 7);
        hIO_contour->SetLineWidth(3);
        hIO_contour->Draw("CONT3 SAME");

        TLegend *legend = new TLegend(0.68, 0.73, 0.88, 0.91);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->AddEntry(hAbsNO, "NO 1#sigma", "l");
        legend->AddEntry(hAbsIO, "IO 1#sigma", "l");
        legend->Draw();
        // add Asimov marker (abs of delm for this plot) and save
        drawAsimovMarkerWithLegend(std::abs(asimovDelm23), asimovDcp, legend);
        cAbs->Update();
        cAbs->SaveAs((outputDir + "delm23_vs_dcp_abs_IO_NO" + suffix + ".png").c_str());
    };

    drawAbsDelm23Overlay("", "delm2_23 < 0", "delm2_23 > 0");

    if (hasUmbrellaWeight) {
        drawAbsDelm23Overlay("_weighted",
                             "(delm2_23 < 0) * umbrella_weight",
                             "(delm2_23 > 0) * umbrella_weight");
    }

    // Reset palette to default for all other plots
    gStyle->SetPalette(kBird);

    // Create canvas for Delta m23 NO vs sin^2(theta23)
    TCanvas *c4 = new TCanvas("c4", "Delta m23 NO vs sin^{2}(#theta_{23})", 900, 700);
    c4->SetRightMargin(0.15);
    // Create 2D histogram: delta m23 vs delta CP
    TH2F *h2d_sin = new TH2F("h2d_sin", "#Delta m_{23}^{2} vs sin^{2}(#theta_{23});sin^{2}(#theta_{23});#Delta m_{23}^{2} (eV^{2})",
        nBins_sin2th23, sin2th23_min, sin2th23_max, nBins_delm23, delm23_no_min, delm23_no_max);
    // Draw from tree (note: y:x format)
    tree->Draw("delm2_23:sin2th_23>>h2d_sin", "delm2_23 > 0", "COLZ");
    normalizeHist(h2d_sin);
    h2d_sin->Draw("COLZ");
    // Styling
    h2d_sin->SetStats(0);  // Remove stats box
    h2d_sin->GetXaxis()->SetTitle("sin^{2}(#theta_{23})");
    h2d_sin->GetYaxis()->SetTitle("#Delta m_{23}^{2} (eV^{2})");
    h2d_sin->GetZaxis()->SetTitle("Normalized Entries");
    //gPad->SetLogz();  // Optional: log scale on z-axis for better visibility
    // Update canvas
    // draw Asimov marker and Save plot
    drawAsimovMarkerWithLegend(asimovSin2th23, asimovDelm23, nullptr);
    c4->Update();
    c4->SaveAs((outputDir + "th23_vs_delm23_no.png").c_str());

    if (hasUmbrellaWeight){
        tree->Draw("delm2_23:sin2th_23>>h2d_sin", "(delm2_23 > 0) * umbrella_weight", "COLZ");
        normalizeHist(h2d_sin);
        h2d_sin->Draw("COLZ");
        drawAsimovMarkerWithLegend(asimovSin2th23, asimovDelm23, nullptr);
        c4->Update();
        c4->SaveAs((outputDir + "th23_vs_delm23_no_weighted.png").c_str());
    }

    // Create canvas for Delta m23 NO vs sin^2(theta23)
    TCanvas *c5 = new TCanvas("c5", "Delta m23 IO vs sin^{2}(#theta_{23})", 900, 700);
    c5->SetRightMargin(0.15);
    // Create 2D histogram: delta m23 vs delta CP
    TH2F *h2d_sin_IO = new TH2F("h2d_sin_IO", "#Delta m_{23}^{2} vs sin^{2}(#theta_{23});sin^{2}(#theta_{23});#Delta m_{23}^{2} (eV^{2})",
                                 nBins_sin2th23, sin2th23_min, sin2th23_max, 
                                 nBins_delm23, delm23_io_min, delm23_io_max);
    // Draw from tree (note: y:x format)
    tree->Draw("delm2_23:sin2th_23>>h2d_sin_IO", "delm2_23 < 0", "COLZ");
    normalizeHist(h2d_sin_IO);
    h2d_sin_IO->Draw("COLZ");
    // Styling
    h2d_sin_IO->SetStats(0);  // Remove stats box
    h2d_sin_IO->GetXaxis()->SetTitle("sin^{2}(#theta_{23})");
    h2d_sin_IO->GetYaxis()->SetTitle("#Delta m_{23}^{2} (eV^{2})");
    h2d_sin_IO->GetZaxis()->SetTitle("Normalized Entries");
    //gPad->SetLogz();  // Optional: log scale on z-axis for better visibility
    // Update canvas
    drawAsimovMarkerWithLegend(asimovSin2th23, asimovDelm23, nullptr);
    c5->Update();
    c5->SaveAs((outputDir + "th23_vs_delm23_io.png").c_str());

    if (hasUmbrellaWeight){
        tree->Draw("delm2_23:sin2th_23>>h2d_sin_IO", "(delm2_23 < 0)*umbrella_weight", "COLZ");
        normalizeHist(h2d_sin_IO);
        h2d_sin_IO->Draw("COLZ");
        drawAsimovMarkerWithLegend(asimovSin2th23, asimovDelm23, nullptr);
        c5->Update();
        c5->SaveAs((outputDir + "th23_vs_delm23_io_weighted.png").c_str());
    }

     // Create canvas for delta_cp vs sin^2(theta23)
    TCanvas *c6 = new TCanvas("c6", "Delta_cp vs sin^{2}(#theta_{23})", 900, 700);
    c6->SetRightMargin(0.15);
    // Create 2D histogram: delta m23 vs delta CP
    TH2F *h2d_delta_cp = new TH2F("h2d_delta_cp", "#delta_{CP} vs sin^{2}(#theta_{23});#delta_{CP};sin^{2}(#theta_{23})",
                                  nBins_dcp, -3.1415, 3.1415, 
                                  nBins_sin2th23, sin2th23_min, sin2th23_max);
    // Draw from tree 
    tree->Draw("sin2th_23:delta_cp>>h2d_delta_cp", "", "COLZ");
    normalizeHist(h2d_delta_cp);
    h2d_delta_cp->Draw("COLZ");
    // Styling
    h2d_delta_cp->SetStats(0);  // Remove stats box
    h2d_delta_cp->GetYaxis()->SetTitle("sin^{2}(#theta_{23})");
    h2d_delta_cp->GetXaxis()->SetTitle("#delta_{CP}");
    h2d_delta_cp->GetZaxis()->SetTitle("Normalized Entries");
    //gPad->SetLogz();  // Optional: log scale on z-axis for better visibility
    // Update canvas
    drawAsimovMarkerWithLegend(asimovDcp, asimovSin2th23, nullptr);
    c6->Update();
    c6->SaveAs((outputDir + "th23_vs_delta_cp.png").c_str());

    if (hasUmbrellaWeight){
        tree->Draw("sin2th_23:delta_cp>>h2d_delta_cp", "umbrella_weight", "COLZ");
        normalizeHist(h2d_delta_cp);
        h2d_delta_cp->Draw("COLZ");
        drawAsimovMarkerWithLegend(asimovDcp, asimovSin2th23, nullptr);
        c6->Update();
        c6->SaveAs((outputDir + "th23_vs_delta_cp_weighted.png").c_str());
    }

    // Create canvas for delta_cp vs sin^2(theta13)
    TCanvas *cth13 = new TCanvas("cth13", "Delta_cp vs sin^{2}(#theta_{13})", 900, 700);
    cth13->SetRightMargin(0.15);
    // Create 2D histogram: delta m23 vs delta CP
    TH2F *h2d_delta_cp_13 = new TH2F("h2d_delta_cp_13", "#delta_{CP} vs sin^{2}(#theta_{13});#delta_{CP};sin^{2}(#theta_{13})",
                                  nBins_dcp, -3.1415, 3.1415, 
                                  nBins_sin2th13, sin2th13_min, sin2th13_max);
    // Draw from tree 
    tree->Draw("sin2th_13:delta_cp>>h2d_delta_cp_13", "", "COLZ");
    normalizeHist(h2d_delta_cp_13);
    h2d_delta_cp_13->Draw("COLZ");
    // Styling
    h2d_delta_cp_13->SetStats(0);  // Remove stats box
    h2d_delta_cp_13->GetYaxis()->SetTitle("sin^{2}(#theta_{13})");
    h2d_delta_cp_13->GetXaxis()->SetTitle("#delta_{CP}");
    h2d_delta_cp_13->GetZaxis()->SetTitle("Normalized Entries");
    //gPad->SetLogz();  // Optional: log scale on z-axis for better visibility
    // Update canvas
    drawAsimovMarkerWithLegend(asimovDcp, asimovSin2th13, nullptr);
    cth13->Update();
    cth13->SaveAs((outputDir + "th13_vs_delta_cp.png").c_str());

    if (hasUmbrellaWeight){
        tree->Draw("sin2th_13:delta_cp>>h2d_delta_cp_13", "umbrella_weight", "COLZ");
        normalizeHist(h2d_delta_cp_13);
        h2d_delta_cp_13->Draw("COLZ");
        drawAsimovMarkerWithLegend(asimovDcp, asimovSin2th13, nullptr);
        cth13->Update();
        cth13->SaveAs((outputDir + "th13_vs_delta_cp_weighted.png").c_str());
    }

    // plot 1D histograms of delta_cp, delm23, theta23, theta13, theta12, delm12
    TCanvas *c7 = new TCanvas("c7", "1D Histograms", 3200, 1200);
    c7->Divide(4, 2);

    // delta_cp unweighted
    c7->cd(1);
    TH1F *h1d_delta_cp_unweighted = new TH1F("h1d_delta_cp_unweighted", "#delta_{CP};#delta_{CP};Entries", nBins_dcp, -3.1415, 3.1415);
    tree->Draw("delta_cp>>h1d_delta_cp_unweighted", "", "");
    normalizeHist(h1d_delta_cp_unweighted);
    h1d_delta_cp_unweighted->SetMinimum(0);
    h1d_delta_cp_unweighted->SetStats(0);
    h1d_delta_cp_unweighted->Draw("HIST");
    h1d_delta_cp_unweighted->Draw("E SAME");
    drawAsimovLineWithLegend(asimovDcp);

    c7->cd(2);
    TH1F *h1d_delm23_all_unweighted = new TH1F("h1d_delm23_all_unweighted", "#Delta m_{23}^{2} (full range);#Delta m_{23}^{2} (eV^{2});Normalized Entries", nBins_delm23, delm23_io_min, delm23_no_max);
    tree->Draw("delm2_23>>h1d_delm23_all_unweighted", "", "");
    normalizeHist(h1d_delm23_all_unweighted);
    h1d_delm23_all_unweighted->SetStats(0);
    h1d_delm23_all_unweighted->Draw("HIST");
    h1d_delm23_all_unweighted->Draw("E SAME");
    drawAsimovLineWithLegend(asimovDelm23);

    c7->cd(3);
    TH1F *h1d_delm23_unweighted = (TH1F*)h1d_delm23_all_unweighted->Clone("h1d_delm23_IO_unweighted");
    h1d_delm23_unweighted->SetTitle("#Delta m_{23}^{2} (IO)");
    h1d_delm23_unweighted->SetStats(0);
    h1d_delm23_unweighted->GetXaxis()->SetRangeUser(delm23_io_min, delm23_io_max);
    h1d_delm23_unweighted->Draw("HIST");
    h1d_delm23_unweighted->Draw("E SAME");
    drawAsimovLineWithLegend(asimovDelm23);

    c7->cd(4);
    TH1F *h1d_delm23_NO_unweighted = (TH1F*)h1d_delm23_all_unweighted->Clone("h1d_delm23_NO_unweighted");
    h1d_delm23_unweighted->SetTitle("#Delta m_{23}^{2} (NO)");
    h1d_delm23_NO_unweighted->SetStats(0);
    h1d_delm23_NO_unweighted->GetXaxis()->SetRangeUser(delm23_no_min, delm23_no_max);
    h1d_delm23_NO_unweighted->Draw("HIST");
    h1d_delm23_NO_unweighted->Draw("E SAME");
    drawAsimovLineWithLegend(asimovDelm23);

    c7->cd(5);
    TH1F *h1d_sinth23_unweighted = new TH1F("h1d_sinth23_unweighted", "sin^{2}(#theta_{23});sin^{2}(#theta_{23});Entries", nBins_sin2th23, sin2th23_min, sin2th23_max);
    tree->Draw("sin2th_23>>h1d_sinth23_unweighted", "", "");
    normalizeHist(h1d_sinth23_unweighted);
    h1d_sinth23_unweighted->SetStats(0);
    h1d_sinth23_unweighted->Draw("HIST");
    h1d_sinth23_unweighted->Draw("E SAME");
    drawAsimovLineWithLegend(asimovSin2th23);

    c7->cd(6);
    TH1F *h1d_sinth13_unweighted = new TH1F("h1d_sinth13_unweighted", "sin^{2}(#theta_{13});sin^{2}(#theta_{13});Entries", nBins_sin2th13, sin2th13_min, sin2th13_max);
    tree->Draw("sin2th_13>>h1d_sinth13_unweighted", "", "");
    normalizeHist(h1d_sinth13_unweighted);
    h1d_sinth13_unweighted->SetStats(0);
    h1d_sinth13_unweighted->Draw("HIST");
    h1d_sinth13_unweighted->Draw("E SAME");
    drawAsimovLineWithLegend(asimovSin2th13);

    c7->cd(7);
    TH1F *h1d_sinth12_unweighted = new TH1F("h1d_sinth12_unweighted", "sin^{2}(#theta_{12});sin^{2}(#theta_{12});Entries", nBins_sin2th12, sin2th12_min, sin2th12_max);
    tree->Draw("sin2th_12>>h1d_sinth12_unweighted", "", "");
    normalizeHist(h1d_sinth12_unweighted);
    h1d_sinth12_unweighted->SetStats(0);
    h1d_sinth12_unweighted->Draw("HIST");
    h1d_sinth12_unweighted->Draw("E SAME");
    drawAsimovLineWithLegend(asimovSin2th12);

    c7->cd(8);
    TH1F *h1d_delm12_unweighted = new TH1F("h1d_delm12_unweighted", "#Delta m_{12}^{2};#Delta m_{12}^{2} (eV^{2});Entries", nBins_delm12, delm12_min, delm12_max);
    tree->Draw("delm2_12>>h1d_delm12_unweighted", "", "");
    normalizeHist(h1d_delm12_unweighted);
    h1d_delm12_unweighted->SetStats(0);
    h1d_delm12_unweighted->Draw("HIST");
    h1d_delm12_unweighted->Draw("E SAME");
    drawAsimovLineWithLegend(asimovDelm12);

    c7->Update();
    c7->SaveAs((outputDir + "1D_histograms_unweighted.png").c_str());

    // Now create weighted version with sigma calculations
    TCanvas *c7_weighted = new TCanvas("c7_weighted", "1D Histograms Weighted", 3200, 1200);
    c7_weighted->Divide(4, 2);

    // delta_cp weighted
    c7_weighted->cd(1);
    TH1F *h1d_delta_cp = new TH1F("h1d_delta_cp", "#delta_{CP};#delta_{CP};Weighted Entries", nBins_dcp, -3.1415, 3.1415);
    if (hasUmbrellaWeight) {
        tree->Draw("delta_cp>>h1d_delta_cp", "umbrella_weight", "");
    } else {
        tree->Draw("delta_cp>>h1d_delta_cp", "", "");
        std::cout << "WARNING:::::No umbrella weights found, plotting unweighted delta_cp histogram." << std::endl;
    }
    normalizeHist(h1d_delta_cp);
    h1d_delta_cp->SetStats(0);
    h1d_delta_cp->Sumw2(); // Enable proper error calculation for weighted histograms
    h1d_delta_cp->Draw("HIST");
    h1d_delta_cp->Draw("E SAME");

    

    // add 1,3,5 sigma contour lines
    double sigma1 = 0.6827;
    double sigma3 = 0.9973;
    double sigma5 = 0.99994;

    double sigma1_upper;
    double sigma1_lower;
    double sigma3_upper;
    double sigma3_lower;
    double sigma5_upper;
    double sigma5_lower;

    // Find the HPD (Highest Posterior Density) bin
    int hpdBin = 1;
    double maxEntries = 0;
    for (int j = 1; j <= h1d_delta_cp->GetNbinsX(); j++) {
        if (h1d_delta_cp->GetBinContent(j) > maxEntries) {
            maxEntries = h1d_delta_cp->GetBinContent(j);
            hpdBin = j;
        }
    }
    
    std::cout << "HPD bin: " << hpdBin << " with " << maxEntries << " entries" << std::endl;
    
    // Create a copy of the histogram to work with
    TH1F *h1d_copy = (TH1F*)h1d_delta_cp->Clone("h1d_copy");
    
    // Track which bins have been included
    std::vector<bool> included(h1d_delta_cp->GetNbinsX() + 1, false);
    included[hpdBin] = true;
    
    // Helper function to find contiguous regions with a gap tolerance
    auto findContiguousRegions = [&](int gapTolerance = 2) {
        std::vector<std::pair<int, int>> regions; // pairs of (lower_bin, upper_bin)
        
        int start = -1;
        int gapCount = 0;
        
        for (int k = 1; k <= h1d_delta_cp->GetNbinsX(); k++) {
            if (included[k]) {
                if (start == -1) {
                    start = k; // start of a new region
                }
                gapCount = 0; // reset gap counter
            } else if (start != -1) {
                gapCount++;
                if (gapCount > gapTolerance) {
                    // end of region
                    regions.push_back({start, k - gapCount});
                    start = -1;
                    gapCount = 0;
                }
            }
        }
        
        // Close last region if it exists
        if (start != -1) {
            regions.push_back({start, h1d_delta_cp->GetNbinsX()});
        }
        
        return regions;
    };
    
    // loop over the bins from highest to lowest summing the bins until we reach the desired sigma levels
    double totalEntries = h1d_delta_cp->Integral();
    double sumEntries = h1d_delta_cp->GetBinContent(hpdBin);
    double level1 = 0;
    double level3 = 0;
    double level5 = 0;
    
    std::vector<std::pair<double, double>> sigma1_regions;
    std::vector<std::pair<double, double>> sigma3_regions;
    std::vector<std::pair<double, double>> sigma5_regions;
    
    for (int i = 0; i < h1d_delta_cp->GetNbinsX(); i++) {
        // find the bin with the highest entries that hasn't been included yet
        int maxBin = 0;
        double maxBinEntries = 0;
        for (int j = 1; j <= h1d_delta_cp->GetNbinsX(); j++) {
            if (!included[j] && h1d_copy->GetBinContent(j) > maxBinEntries) {
                maxBinEntries = h1d_copy->GetBinContent(j);
                maxBin = j;
            }
        }
        
        if (maxBin == 0) break; // no more bins to add
        
        sumEntries += maxBinEntries;
        included[maxBin] = true;
        double frac = sumEntries / totalEntries;
        
        std::cout << "Added bin " << maxBin << " with " << maxBinEntries << " entries. Fraction: " << frac << std::endl;
        
        if (frac >= sigma1 && level1 == 0) {
            level1 = maxBinEntries;
            std::cout << "1 sigma level: " << level1 << std::endl;
            // Find all contiguous regions
            auto regions = findContiguousRegions(2);
            std::cout << "1 sigma regions (" << regions.size() << " mode(s)):" << std::endl;
            for (const auto& region : regions) {
                double lower = h1d_delta_cp->GetBinLowEdge(region.first);
                double upper = h1d_delta_cp->GetBinLowEdge(region.second + 1);
                sigma1_regions.push_back({lower, upper});
                std::cout << "  [" << lower << ", " << upper << "]" << std::endl;
            }
        }
        if (frac >= sigma3 && level3 == 0) {
            level3 = maxBinEntries;
            std::cout << "3 sigma level: " << level3 << std::endl;
            // Find all contiguous regions
            auto regions = findContiguousRegions(2);
            std::cout << "3 sigma regions (" << regions.size() << " mode(s)):" << std::endl;
            for (const auto& region : regions) {
                double lower = h1d_delta_cp->GetBinLowEdge(region.first);
                double upper = h1d_delta_cp->GetBinLowEdge(region.second + 1);
                sigma3_regions.push_back({lower, upper});
                std::cout << "  [" << lower << ", " << upper << "]" << std::endl;
            }
        }
        if (frac >= sigma5 && level5 == 0) {
            level5 = maxBinEntries;
            std::cout << "5 sigma level: " << level5 << std::endl;
            // Find all contiguous regions
            auto regions = findContiguousRegions(2);
            std::cout << "5 sigma regions (" << regions.size() << " mode(s)):" << std::endl;
            for (const auto& region : regions) {
                double lower = h1d_delta_cp->GetBinLowEdge(region.first);
                double upper = h1d_delta_cp->GetBinLowEdge(region.second + 1);
                sigma5_regions.push_back({lower, upper});
                std::cout << "  [" << lower << ", " << upper << "]" << std::endl;
            }
        }
    }

    // Draw vertical lines at the sigma bounds for all regions
    double ymax = h1d_delta_cp->GetMaximum();
    
    // 1 sigma bounds (yellow/orange) - draw all regions
    for (const auto& region : sigma1_regions) {
        TLine *line_lower = new TLine(region.first, 0, region.first, ymax);
        line_lower->SetLineColor(kOrange+1);
        line_lower->SetLineStyle(2);
        line_lower->SetLineWidth(2);
        line_lower->Draw("same");
        
        TLine *line_upper = new TLine(region.second, 0, region.second, ymax);
        line_upper->SetLineColor(kOrange+1);
        line_upper->SetLineStyle(2);
        line_upper->SetLineWidth(2);
        line_upper->Draw("same");
    }

    // 3 sigma bounds (blue) - draw all regions
    for (const auto& region : sigma3_regions) {
        TLine *line_lower = new TLine(region.first, 0, region.first, ymax);
        line_lower->SetLineColor(kAzure+2);
        line_lower->SetLineStyle(2);
        line_lower->SetLineWidth(2);
        line_lower->Draw("same");
        
        TLine *line_upper = new TLine(region.second, 0, region.second, ymax);
        line_upper->SetLineColor(kAzure+2);
        line_upper->SetLineStyle(2);
        line_upper->SetLineWidth(2);
        line_upper->Draw("same");
    }

    // 5 sigma bounds (magenta) - draw all regions
    for (const auto& region : sigma5_regions) {
        TLine *line_lower = new TLine(region.first, 0, region.first, ymax);
        line_lower->SetLineColor(kMagenta+2);
        line_lower->SetLineStyle(2);
        line_lower->SetLineWidth(2);
        line_lower->Draw("same");
        
        TLine *line_upper = new TLine(region.second, 0, region.second, ymax);
        line_upper->SetLineColor(kMagenta+2);
        line_upper->SetLineStyle(2);
        line_upper->SetLineWidth(2);
        line_upper->Draw("same");
    }

    TLegend *legRef = nullptr;
    if (hasReferenceDcp && hRefDcpShape) {
        std::cout << "Overlaying reference delta_cp shape on weighted histogram" << std::endl;
        TH1D *hRefOverlay = (TH1D*)hRefDcpShape->Clone("hRefDcpShape_overlay_linear");
        hRefOverlay->SetDirectory(0);
        double refIntegral = hRefOverlay->Integral();
        if (refIntegral > 0) {
            hRefOverlay->Scale(h1d_delta_cp->Integral() / refIntegral);
        }
        hRefOverlay->SetLineColor(kMagenta + 2);
        hRefOverlay->SetLineWidth(3);
        hRefOverlay->SetLineStyle(1);
        hRefOverlay->SetFillStyle(0);
        hRefOverlay->Draw("HIST SAME");

        TLegend *legRef = new TLegend(0.15, 0.74, 0.48, 0.88);
        legRef->SetBorderSize(0);
        legRef->SetFillStyle(0);
        legRef->AddEntry(h1d_delta_cp, "Umbrella weighted", "l");
        legRef->AddEntry(hRefOverlay, "Reference dcp (shape, area-normalized)", "l");
        legRef->Draw();
    }
    drawAsimovLineWithLegend(asimovDcp, hasReferenceDcp ? legRef : nullptr);

    // If LLH-scan reference histograms were provided, overlay them (they were area-normalized above)
    if (hasReferenceLLH_Dcp && hRefLLH_Dcp) {
        std::cout << "Overlaying LLH-scan reference dcp histogram on weighted histogram (resampled to umbrella bins)" << std::endl;
        // Resample LLH histogram to match umbrella binning
        TH1D *hLLHResamp = new TH1D("hRefLLH_dcp_resamp", "", h1d_delta_cp->GetNbinsX(), h1d_delta_cp->GetXaxis()->GetXmin(), h1d_delta_cp->GetXaxis()->GetXmax());
        hLLHResamp->SetDirectory(0);
        for (int i = 1; i <= hLLHResamp->GetNbinsX(); ++i) {
            double xlow = hLLHResamp->GetBinLowEdge(i);
            double xhigh = hLLHResamp->GetBinLowEdge(i+1);
            double val = 0.0;
            int nOrig = hRefLLH_Dcp->GetNbinsX();
            for (int k = 1; k <= nOrig; ++k) {
                double origLow = hRefLLH_Dcp->GetXaxis()->GetBinLowEdge(k);
                double origHigh = origLow + hRefLLH_Dcp->GetBinWidth(k);
                double overlap = std::max(0.0, std::min(origHigh, xhigh) - std::max(origLow, xlow));
                if (overlap > 0.0) {
                    double origWidth = hRefLLH_Dcp->GetBinWidth(k);
                    double contrib = hRefLLH_Dcp->GetBinContent(k) * (overlap / origWidth);
                    val += contrib;
                }
            }
            hLLHResamp->SetBinContent(i, val);
        }
        // Now scale resampled LLH to umbrella histogram area
        double llhInt = hLLHResamp->Integral();
        if (llhInt > 0) hLLHResamp->Scale(h1d_delta_cp->Integral() / llhInt);
        hLLHResamp->SetLineColor(kBlack);
        hLLHResamp->SetLineWidth(3);
        hLLHResamp->SetLineStyle(2);
        hLLHResamp->SetFillStyle(0);
        hLLHResamp->Draw("HIST SAME");

        TLegend *legLLH = new TLegend(0.15, 0.60, 0.50, 0.74);
        legLLH->SetBorderSize(0);
        legLLH->SetFillStyle(0);
        legLLH->AddEntry(h1d_delta_cp, "Umbrella weighted", "l");
        legLLH->AddEntry(hLLHResamp, "LLH scan (resampled)", "l");
        legLLH->Draw();
    }

    c7_weighted->cd(2);
    TH1F *h1d_delm23_all = new TH1F("h1d_delm23_all", "#Delta m_{23}^{2} (full range);#Delta m_{23}^{2} (eV^{2});Normalized Entries", nBins_delm23, delm23_io_min, delm23_no_max);
    if (hasUmbrellaWeight) {
        tree->Draw("delm2_23>>h1d_delm23_all", "umbrella_weight", "");
        std::cout << "Applying umbrella weights to delm23 histogram" << std::endl;
    } else {
        tree->Draw("delm2_23>>h1d_delm23_all", "", "");
        std::cout << "WARNING:::::No umbrella weights found, plotting unweighted #Delta m_{23}^{2} histogram." << std::endl;
    }
    normalizeHist(h1d_delm23_all);
    h1d_delm23_all->SetStats(0);
    h1d_delm23_all->Draw("HIST");
    h1d_delm23_all->Draw("E SAME");

    TLegend *legRefDelm23All = nullptr;
    if (hasReferenceDelm23 && hRefDelm23Shape) {
        TH1D *hRefDelm23AllOverlay = (TH1D*)hRefDelm23Shape->Clone("hRefDelm23All_overlay_linear");
        hRefDelm23AllOverlay->SetDirectory(0);
        double refIntegral = hRefDelm23AllOverlay->Integral();
        if (refIntegral > 0) {
            hRefDelm23AllOverlay->Scale(h1d_delm23_all->Integral() / refIntegral);
        }
        hRefDelm23AllOverlay->SetLineColor(kMagenta + 2);
        hRefDelm23AllOverlay->SetLineWidth(3);
        hRefDelm23AllOverlay->SetLineStyle(1);
        hRefDelm23AllOverlay->SetFillStyle(0);
        hRefDelm23AllOverlay->Draw("HIST SAME");

        TLegend *legRefDelm23All = new TLegend(0.50, 0.74, 0.88, 0.88);
        legRefDelm23All->SetBorderSize(0);
        legRefDelm23All->SetFillStyle(0);
        legRefDelm23All->AddEntry(h1d_delm23_all, "Umbrella weighted", "l");
        legRefDelm23All->AddEntry(hRefDelm23AllOverlay, "Reference #Delta m_{23}^{2} (full range)", "l");
        legRefDelm23All->Draw();
    }
    drawAsimovLineWithLegend(asimovDelm23, hasReferenceDelm23 ? legRefDelm23All : nullptr);

    c7_weighted->cd(3);
    TH1F *h1d_delm23 = (TH1F*)h1d_delm23_all->Clone("h1d_delm23_IO");
    h1d_delm23->SetStats(0);
    h1d_delm23->Sumw2();
    h1d_delm23->GetXaxis()->SetRangeUser(delm23_io_min, delm23_io_max);
    h1d_delm23->Draw("HIST");
    h1d_delm23->Draw("E SAME");

    TLegend *legRefDelm23IO = nullptr;
    if (hasReferenceDelm23 && hRefDelm23Shape) {
        TH1D *hRefDelm23IOOverlay = (TH1D*)hRefDelm23Shape->Clone("hRefDelm23IO_overlay_linear");
        hRefDelm23IOOverlay->SetDirectory(0);
        hRefDelm23IOOverlay->GetXaxis()->SetRangeUser(delm23_io_min, delm23_io_max);
        hRefDelm23IOOverlay->SetLineColor(kMagenta + 2);
        hRefDelm23IOOverlay->SetLineWidth(3);
        hRefDelm23IOOverlay->SetLineStyle(1);
        hRefDelm23IOOverlay->SetFillStyle(0);
        hRefDelm23IOOverlay->Draw("HIST SAME");

        TLegend *legRefDelm23IO = new TLegend(0.50, 0.74, 0.88, 0.88);
        legRefDelm23IO->SetBorderSize(0);
        legRefDelm23IO->SetFillStyle(0);
        legRefDelm23IO->AddEntry(h1d_delm23, "Umbrella weighted", "l");
        legRefDelm23IO->AddEntry(hRefDelm23IOOverlay, "Reference #Delta m_{23}^{2} (IO)", "l");
        legRefDelm23IO->Draw();
    }
    drawAsimovLineWithLegend(asimovDelm23, hasReferenceDelm23 ? legRefDelm23IO : nullptr);

    c7_weighted->cd(4);
    TH1F *h1d_delm23_NO = (TH1F*)h1d_delm23_all->Clone("h1d_delm23_NO");
    h1d_delm23_NO->SetStats(0);
    h1d_delm23_NO->Sumw2();
    h1d_delm23_NO->GetXaxis()->SetRangeUser(delm23_no_min, delm23_no_max);
    h1d_delm23_NO->Draw("HIST");
    h1d_delm23_NO->Draw("E SAME");

    TLegend *legRefDelm23NO = nullptr;
    if (hasReferenceDelm23 && hRefDelm23Shape) {
        std::cout << "Drawing reference delm23 shape for NO" << std::endl;
        TH1D *hRefDelm23NOOverlay = (TH1D*)hRefDelm23Shape->Clone("hRefDelm23NO_overlay_linear");
        hRefDelm23NOOverlay->SetDirectory(0);
        hRefDelm23NOOverlay->GetXaxis()->SetRangeUser(delm23_no_min, delm23_no_max);
        hRefDelm23NOOverlay->SetLineColor(kMagenta + 2);
        hRefDelm23NOOverlay->SetLineWidth(3);
        hRefDelm23NOOverlay->SetLineStyle(1);
        hRefDelm23NOOverlay->SetFillStyle(0);
        hRefDelm23NOOverlay->Draw("HIST SAME");

        TLegend *legRefDelm23NO = new TLegend(0.50, 0.74, 0.88, 0.88);
        legRefDelm23NO->SetBorderSize(0);
        legRefDelm23NO->SetFillStyle(0);
        legRefDelm23NO->AddEntry(h1d_delm23_NO, "Umbrella weighted", "l");
        legRefDelm23NO->AddEntry(hRefDelm23NOOverlay, "Reference #Delta m_{23}^{2} (NO)", "l");
        legRefDelm23NO->Draw();
    }
    drawAsimovLineWithLegend(asimovDelm23, hasReferenceDelm23 ? legRefDelm23NO : nullptr);

    // Compute and print mass ordering preference (NO/IO) using the globally-normalized full-range histogram
    // Use `h1d_delm23_all` which has been normalized so its integral == 1, then sum bins over IO/NO ranges.
    if (h1d_delm23_all) {
        int binIO_low = h1d_delm23_all->GetXaxis()->FindBin(delm23_io_min);
        int binIO_high = h1d_delm23_all->GetXaxis()->FindBin(delm23_io_max);
        int binNO_low = h1d_delm23_all->GetXaxis()->FindBin(delm23_no_min);
        int binNO_high = h1d_delm23_all->GetXaxis()->FindBin(delm23_no_max);

        double fracIO = h1d_delm23_all->Integral(binIO_low, binIO_high);
        double fracNO = h1d_delm23_all->Integral(binNO_low, binNO_high);

        if (fracIO > 0) {
            double ratioNOtoIO = fracNO / fracIO;
            std::cout << "Mass ordering preference (NO/IO) = " << ratioNOtoIO
                      << "  (NO frac=" << fracNO << ", IO frac=" << fracIO << ")" << std::endl;
        } else {
            std::cout << "Mass ordering preference: IO fraction is zero; cannot compute ratio. "
                      << "NO_frac=" << fracNO << " IO_frac=" << fracIO << std::endl;
        }
    } else {
        std::cout << "Mass ordering preference: h1d_delm23_all is null; cannot compute." << std::endl;
    }

    c7_weighted->cd(5);
    TH1F *h1d_sinth23 = new TH1F("h1d_sinth23", "sin^{2}(#theta_{23});sin^{2}(#theta_{23});Weighted Entries", nBins_sin2th23, sin2th23_min, sin2th23_max);
    if (hasUmbrellaWeight) {
        tree->Draw("sin2th_23>>h1d_sinth23", "umbrella_weight", "");
    } else {
        tree->Draw("sin2th_23>>h1d_sinth23", "", "");
        std::cout << "WARNING:::::No umbrella weights found, plotting unweighted sin^2(theta23) histogram." << std::endl;
    }
    normalizeHist(h1d_sinth23);
    h1d_sinth23->SetStats(0);
    h1d_sinth23->Sumw2();
    h1d_sinth23->Draw("HIST");
    h1d_sinth23->Draw("E SAME");

    TLegend *legRefSin = nullptr;
    if (hasReferenceSin2Th23 && hRefSin2Th23Shape) {
        TH1D *hRefSinOverlay = (TH1D*)hRefSin2Th23Shape->Clone("hRefSin2Th23_overlay_linear");
        hRefSinOverlay->SetDirectory(0);
        double refIntegral = hRefSinOverlay->Integral();
        if (refIntegral > 0) {
            hRefSinOverlay->Scale(h1d_sinth23->Integral() / refIntegral);
        }
        hRefSinOverlay->SetLineColor(kMagenta + 2);
        hRefSinOverlay->SetLineWidth(3);
        hRefSinOverlay->SetLineStyle(1);
        hRefSinOverlay->SetFillStyle(0);
        hRefSinOverlay->Draw("HIST SAME");

        TLegend *legRefSin = new TLegend(0.50, 0.74, 0.88, 0.88);
        legRefSin->SetBorderSize(0);
        legRefSin->SetFillStyle(0);
        legRefSin->AddEntry(h1d_sinth23, "Umbrella weighted", "l");
        legRefSin->AddEntry(hRefSinOverlay, "Reference sin^{2}(#theta_{23})", "l");
        legRefSin->Draw();
    }
    drawAsimovLineWithLegend(asimovSin2th23, hasReferenceSin2Th23 ? legRefSin : nullptr);

    c7_weighted->cd(6);
    TH1F *h1d_sinth13 = new TH1F("h1d_sinth13", "sin^{2}(#theta_{13});sin^{2}(#theta_{13});Weighted Entries", nBins_sin2th13, sin2th13_min, sin2th13_max);
    if (hasUmbrellaWeight) {
        tree->Draw("sin2th_13>>h1d_sinth13", "umbrella_weight", "");
    } else {
        tree->Draw("sin2th_13>>h1d_sinth13", "", "");
        std::cout << "WARNING:::::No umbrella weights found, plotting unweighted sin^2(theta13) histogram." << std::endl;
    }
    normalizeHist(h1d_sinth13);
    h1d_sinth13->SetStats(0);
    h1d_sinth13->Sumw2();
    h1d_sinth13->Draw("HIST");
    h1d_sinth13->Draw("E SAME");
    drawAsimovLineWithLegend(asimovSin2th13);

    c7_weighted->cd(7);
    TH1F *h1d_sinth12 = new TH1F("h1d_sinth12", "sin^{2}(#theta_{12});sin^{2}(#theta_{12});Weighted Entries", nBins_sin2th12, sin2th12_min, sin2th12_max);
    if (hasUmbrellaWeight) {
        tree->Draw("sin2th_12>>h1d_sinth12", "umbrella_weight", "");
    } else {
        tree->Draw("sin2th_12>>h1d_sinth12", "", "");
        std::cout << "WARNING:::::No umbrella weights found, plotting unweighted sin^2(theta12) histogram." << std::endl;
    }
    normalizeHist(h1d_sinth12);
    h1d_sinth12->SetStats(0);
    h1d_sinth12->Sumw2();
    h1d_sinth12->Draw("HIST");
    h1d_sinth12->Draw("E SAME");
    drawAsimovLineWithLegend(asimovSin2th12);

    c7_weighted->cd(8);
    TH1F *h1d_delm12 = new TH1F("h1d_delm12", "#Delta m_{12}^{2};#Delta m_{12}^{2} (eV^{2});Weighted Entries", nBins_delm12, delm12_min, delm12_max);
    if (hasUmbrellaWeight) {
        tree->Draw("delm2_12>>h1d_delm12", "umbrella_weight", "");
    } else {
        tree->Draw("delm2_12>>h1d_delm12", "", "");
        std::cout << "WARNING:::::No umbrella weights found, plotting unweighted #Delta m_{12}^{2} histogram." << std::endl;
    }
    normalizeHist(h1d_delm12);
    h1d_delm12->SetStats(0);
    h1d_delm12->Sumw2();
    h1d_delm12->Draw("HIST");
    h1d_delm12->Draw("E SAME");
    drawAsimovLineWithLegend(asimovDelm12);
    
    c7_weighted->Update();
    c7_weighted->SaveAs((outputDir + "1D_histograms_weighted.png").c_str());
    
    // Create separate canvas for delta_cp with log scale
    TCanvas *c8 = new TCanvas("c8", "Delta CP Log Scale", 900, 700);
    c8->SetLogy();
    TH1F *h1d_delta_cp_log = (TH1F*)h1d_delta_cp->Clone("h1d_delta_cp_log");
    h1d_delta_cp_log->SetTitle("#delta_{CP} (Log Scale);#delta_{CP};Weighted Entries");
    h1d_delta_cp_log->Draw("HIST");
    h1d_delta_cp_log->Draw("E SAME");
    
    // Draw vertical lines at the sigma bounds for all regions on log scale plot
    double ymax_log = h1d_delta_cp_log->GetMaximum();
    double ymin_log = h1d_delta_cp_log->GetMinimum();
    if (ymin_log <= 0) ymin_log = 0.1; // avoid log(0)
    
    // 1 sigma bounds (user colors) - draw one-sided lines and HPD regions
    for (const auto& region : sigma1_regions) {
        TLine *line_lower = new TLine(region.first, ymin_log, region.first, ymax_log);
        line_lower->SetLineColor(kRed);
        line_lower->SetLineStyle(2);
        line_lower->SetLineWidth(2);
        line_lower->Draw("same");
        
        TLine *line_upper = new TLine(region.second, ymin_log, region.second, ymax_log);
        line_upper->SetLineColor(kRed);
        line_upper->SetLineStyle(2);
        line_upper->SetLineWidth(2);
        line_upper->Draw("same");
    }

    // 3 sigma bounds (turquoise/darker blue) - draw all regions
    for (const auto& region : sigma3_regions) {
        TLine *line_lower = new TLine(region.first, ymin_log, region.first, ymax_log);
        line_lower->SetLineColor(kAzure+2);
        line_lower->SetLineStyle(2);
        line_lower->SetLineWidth(2);
        line_lower->Draw("same");
        
        TLine *line_upper = new TLine(region.second, ymin_log, region.second, ymax_log);
        line_upper->SetLineColor(kAzure+2);
        line_upper->SetLineStyle(2);
        line_upper->SetLineWidth(2);
        line_upper->Draw("same");
    }

    // 5 sigma bounds (purple) - draw all regions
    for (const auto& region : sigma5_regions) {
        TLine *line_lower = new TLine(region.first, ymin_log, region.first, ymax_log);
        line_lower->SetLineColor(kMagenta+2);
        line_lower->SetLineStyle(2);
        line_lower->SetLineWidth(2);
        line_lower->Draw("same");
        
        TLine *line_upper = new TLine(region.second, ymin_log, region.second, ymax_log);
        line_upper->SetLineColor(kMagenta+2);
        line_upper->SetLineStyle(2);
        line_upper->SetLineWidth(2);
        line_upper->Draw("same");
    }

    TLegend *legRefLog = nullptr;
    if (hasReferenceDcp && hRefDcpShape) {
        TH1D *hRefOverlayLog = (TH1D*)hRefDcpShape->Clone("hRefDcpShape_overlay_log");
        hRefOverlayLog->SetDirectory(0);
        double refIntegral = hRefOverlayLog->Integral();
        if (refIntegral > 0) {
            hRefOverlayLog->Scale(h1d_delta_cp_log->Integral() / refIntegral);
        }
        hRefOverlayLog->SetLineColor(kMagenta + 2);
        hRefOverlayLog->SetLineWidth(3);
        hRefOverlayLog->SetLineStyle(1);
        hRefOverlayLog->SetFillStyle(0);
        hRefOverlayLog->Draw("HIST SAME");

        TLegend *legRefLog = new TLegend(0.50, 0.74, 0.88, 0.88);
        legRefLog->SetBorderSize(0);
        legRefLog->SetFillStyle(0);
        legRefLog->AddEntry(h1d_delta_cp_log, "Umbrella weighted", "l");
        legRefLog->AddEntry(hRefOverlayLog, "Reference dcp (shape, area-normalized)", "l");
        legRefLog->Draw();
    }
    if (hasReferenceLLH_Dcp && hRefLLH_Dcp) {
        // Resample LLH histogram to match log-plot binning
        TH1D *hLLHResampLog = new TH1D("hRefLLH_dcp_resamp_log", "", h1d_delta_cp_log->GetNbinsX(), h1d_delta_cp_log->GetXaxis()->GetXmin(), h1d_delta_cp_log->GetXaxis()->GetXmax());
        hLLHResampLog->SetDirectory(0);
        for (int i = 1; i <= hLLHResampLog->GetNbinsX(); ++i) {
            double xlow = hLLHResampLog->GetBinLowEdge(i);
            double xhigh = hLLHResampLog->GetBinLowEdge(i+1);
            double val = 0.0;
            int nOrig = hRefLLH_Dcp->GetNbinsX();
            for (int k = 1; k <= nOrig; ++k) {
                double origLow = hRefLLH_Dcp->GetXaxis()->GetBinLowEdge(k);
                double origHigh = origLow + hRefLLH_Dcp->GetBinWidth(k);
                double overlap = std::max(0.0, std::min(origHigh, xhigh) - std::max(origLow, xlow));
                if (overlap > 0.0) {
                    double origWidth = hRefLLH_Dcp->GetBinWidth(k);
                    double contrib = hRefLLH_Dcp->GetBinContent(k) * (overlap / origWidth);
                    val += contrib;
                }
            }
            hLLHResampLog->SetBinContent(i, val);
        }
        double llhI = hLLHResampLog->Integral();
        if (llhI > 0) hLLHResampLog->Scale(h1d_delta_cp_log->Integral() / llhI);
        hLLHResampLog->SetLineColor(kBlue + 2);
        hLLHResampLog->SetLineWidth(3);
        hLLHResampLog->SetLineStyle(2);
        hLLHResampLog->SetFillStyle(0);
        hLLHResampLog->Draw("HIST SAME");
        TLegend *legLLHLog = new TLegend(0.50, 0.60, 0.88, 0.74);
        legLLHLog->SetBorderSize(0);
        legLLHLog->SetFillStyle(0);
        legLLHLog->AddEntry(h1d_delta_cp_log, "Umbrella weighted", "l");
        legLLHLog->AddEntry(hLLHResampLog, "LLH scan (resampled & area-normalized)", "l");
        legLLHLog->Draw();
    }
    drawAsimovLineWithLegend(asimovDcp, hasReferenceDcp ? legRefLog : nullptr);
    
    c8->Update();
    c8->SaveAs((outputDir + "delta_cp_log_scale_weighted.png").c_str());

    // Create separate two-panel log-scale delta_cp plots marginalized over IO and NO.
    TCanvas *c9 = new TCanvas("c9", "Delta CP Log Scale Marginalized by Ordering", 1800, 700);
    c9->Divide(2, 1);

    std::string ioSelection = hasUmbrellaWeight ? "(delm2_23 < 0) * umbrella_weight" : "(delm2_23 < 0)";
    std::string noSelection = hasUmbrellaWeight ? "(delm2_23 > 0) * umbrella_weight" : "(delm2_23 > 0)";

    // IO panel
    c9->cd(1);
    gPad->SetLogy();
    TH1F *h1d_delta_cp_io_log = new TH1F("h1d_delta_cp_io_log", "#delta_{CP} marginalized over IO;#delta_{CP};Normalized Entries", nBins_dcp, -3.1415, 3.1415);
    tree->Draw("delta_cp>>h1d_delta_cp_io_log", ioSelection.c_str(), "");
    normalizeHist(h1d_delta_cp_io_log);
    h1d_delta_cp_io_log->SetStats(0);
    double ioMinPositive = 1.0;
    bool ioFoundPositive = false;
    for (int b = 1; b <= h1d_delta_cp_io_log->GetNbinsX(); ++b) {
        double y = h1d_delta_cp_io_log->GetBinContent(b);
        if (y > 0.0) {
            if (!ioFoundPositive || y < ioMinPositive) {
                ioMinPositive = y;
            }
            ioFoundPositive = true;
        }
    }
    h1d_delta_cp_io_log->SetMinimum(ioFoundPositive ? 0.5 * ioMinPositive : 1e-8);
    h1d_delta_cp_io_log->Draw("HIST");
    h1d_delta_cp_io_log->Draw("E SAME");
    drawAsimovLineWithLegend(asimovDcp);
    if (hasReferenceLLH_Dcp_IO && hRefLLH_Dcp_IO) {
        TH1D *hLLHIO = new TH1D("hRefLLH_dcp_io_resamp", "", h1d_delta_cp_io_log->GetNbinsX(), h1d_delta_cp_io_log->GetXaxis()->GetXmin(), h1d_delta_cp_io_log->GetXaxis()->GetXmax());
        hLLHIO->SetDirectory(0);
        for (int i = 1; i <= hLLHIO->GetNbinsX(); ++i) {
            double xlow = hLLHIO->GetBinLowEdge(i);
            double xhigh = hLLHIO->GetBinLowEdge(i+1);
            double val = 0.0;
            int nOrig = hRefLLH_Dcp_IO->GetNbinsX();
            for (int k = 1; k <= nOrig; ++k) {
                double origLow = hRefLLH_Dcp_IO->GetXaxis()->GetBinLowEdge(k);
                double origHigh = origLow + hRefLLH_Dcp_IO->GetBinWidth(k);
                double overlap = std::max(0.0, std::min(origHigh, xhigh) - std::max(origLow, xlow));
                if (overlap > 0.0) {
                    double origWidth = hRefLLH_Dcp_IO->GetBinWidth(k);
                    double contrib = hRefLLH_Dcp_IO->GetBinContent(k) * (overlap / origWidth);
                    val += contrib;
                }
            }
            hLLHIO->SetBinContent(i, val);
        }
        double li = hLLHIO->Integral();
        if (li > 0) hLLHIO->Scale(h1d_delta_cp_io_log->Integral() / li);
        hLLHIO->SetLineColor(kBlue + 2);
        hLLHIO->SetLineWidth(2);
        hLLHIO->SetLineStyle(2);
        hLLHIO->Draw("HIST SAME");
        TLegend *legLLHIO = new TLegend(0.50, 0.60, 0.88, 0.74);
        legLLHIO->SetBorderSize(0);
        legLLHIO->SetFillStyle(0);
        legLLHIO->AddEntry(h1d_delta_cp_io_log, "Umbrella weighted (IO)", "l");
        legLLHIO->AddEntry(hLLHIO, "LLH scan IO (resampled & area-normalized)", "l");
        legLLHIO->Draw();
    }

    // NO panel
    c9->cd(2);
    gPad->SetLogy();
    TH1F *h1d_delta_cp_no_log = new TH1F("h1d_delta_cp_no_log", "#delta_{CP} marginalized over NO;#delta_{CP};Normalized Entries", nBins_dcp, -3.1415, 3.1415);
    tree->Draw("delta_cp>>h1d_delta_cp_no_log", noSelection.c_str(), "");
    normalizeHist(h1d_delta_cp_no_log);
    h1d_delta_cp_no_log->SetStats(0);
    double noMinPositive = 1.0;
    bool noFoundPositive = false;
    for (int b = 1; b <= h1d_delta_cp_no_log->GetNbinsX(); ++b) {
        double y = h1d_delta_cp_no_log->GetBinContent(b);
        if (y > 0.0) {
            if (!noFoundPositive || y < noMinPositive) {
                noMinPositive = y;
            }
            noFoundPositive = true;
        }
    }
    h1d_delta_cp_no_log->SetMinimum(noFoundPositive ? 0.5 * noMinPositive : 1e-8);
    h1d_delta_cp_no_log->Draw("HIST");
    h1d_delta_cp_no_log->Draw("E SAME");
    drawAsimovLineWithLegend(asimovDcp);
    if (hasReferenceLLH_Dcp_NO && hRefLLH_Dcp_NO) {
        TH1D *hLLHNO = new TH1D("hRefLLH_dcp_no_resamp", "", h1d_delta_cp_no_log->GetNbinsX(), h1d_delta_cp_no_log->GetXaxis()->GetXmin(), h1d_delta_cp_no_log->GetXaxis()->GetXmax());
        hLLHNO->SetDirectory(0);
        for (int i = 1; i <= hLLHNO->GetNbinsX(); ++i) {
            double xlow = hLLHNO->GetBinLowEdge(i);
            double xhigh = hLLHNO->GetBinLowEdge(i+1);
            double val = 0.0;
            int nOrig = hRefLLH_Dcp_NO->GetNbinsX();
            for (int k = 1; k <= nOrig; ++k) {
                double origLow = hRefLLH_Dcp_NO->GetXaxis()->GetBinLowEdge(k);
                double origHigh = origLow + hRefLLH_Dcp_NO->GetBinWidth(k);
                double overlap = std::max(0.0, std::min(origHigh, xhigh) - std::max(origLow, xlow));
                if (overlap > 0.0) {
                    double origWidth = hRefLLH_Dcp_NO->GetBinWidth(k);
                    double contrib = hRefLLH_Dcp_NO->GetBinContent(k) * (overlap / origWidth);
                    val += contrib;
                }
            }
            hLLHNO->SetBinContent(i, val);
        }
        double ln = hLLHNO->Integral();
        if (ln > 0) hLLHNO->Scale(h1d_delta_cp_no_log->Integral() / ln);
        hLLHNO->SetLineColor(kBlue + 2);
        hLLHNO->SetLineWidth(2);
        hLLHNO->SetLineStyle(2);
        hLLHNO->Draw("HIST SAME");
        TLegend *legLLHNO = new TLegend(0.50, 0.60, 0.88, 0.74);
        legLLHNO->SetBorderSize(0);
        legLLHNO->SetFillStyle(0);
        legLLHNO->AddEntry(h1d_delta_cp_no_log, "Umbrella weighted (NO)", "l");
        legLLHNO->AddEntry(hLLHNO, "LLH scan NO (resampled & area-normalized)", "l");
        legLLHNO->Draw();
    }

    c9->Update();
    c9->SaveAs((outputDir + "delta_cp_log_scale_marginalized_io_no.png").c_str());

    // Jarlskog invariant plots (unitless)
    const double jarlskogMin = -0.06;
    const double jarlskogMax = 0.06;
    const int nBinsJarlskog = 160;

    TCanvas *c10 = new TCanvas("c10", "Jarlskog Invariant", 900, 700);
    TH1F *h1d_jarlskog_unweighted = new TH1F("h1d_jarlskog_unweighted", "Jarlskog invariant;Jarlskog invariant;Entries",
                                             nBinsJarlskog, jarlskogMin, jarlskogMax);
    tree->Draw("TMath::Sqrt(TMath::Max(0.,sin2th_12)*TMath::Max(0.,1.-sin2th_12)*TMath::Max(0.,sin2th_23)*TMath::Max(0.,1.-sin2th_23))*TMath::Max(0.,1.-sin2th_13)*TMath::Sqrt(TMath::Max(0.,sin2th_13))*TMath::Sin(delta_cp)>>h1d_jarlskog_unweighted",
               "", "");
    normalizeHist(h1d_jarlskog_unweighted);
    h1d_jarlskog_unweighted->SetStats(0);
    h1d_jarlskog_unweighted->Draw("HIST");
    h1d_jarlskog_unweighted->Draw("E SAME");
    drawAsimovLineWithLegend(asimovJarlskog);
    // One-sided confidence intervals for Jarlskog (unweighted)
    auto computeOneSidedBounds = [&](TH1 *hist, double targets[3], double outBounds[3]) {
        // targets: fractions (e.g., 0.6827, 0.9973, 0.99994)
        // outBounds: boundaries in x where cumulative from chosen side reaches target
        for (int i = 0; i < 3; ++i) outBounds[i] = std::numeric_limits<double>::quiet_NaN();
        if (!hist) return;
        double total = hist->Integral();
        if (total <= 0) return;
        int nBins = hist->GetNbinsX();
        // determine which side of zero holds more posterior (simple count of bin content)
        double leftSum = 0.0;
        double rightSum = 0.0;
        for (int b = 1; b <= nBins; ++b) {
            double center = hist->GetBinCenter(b);
            double c = hist->GetBinContent(b);
            if (center < 0) leftSum += c; else rightSum += c;
        }
        bool integrateLeftToRight = (leftSum >= rightSum);

        for (int t = 0; t < 3; ++t) {
            double target = targets[t];
            double cum = 0.0;
            if (integrateLeftToRight) {
                for (int b = 1; b <= nBins; ++b) {
                    cum += hist->GetBinContent(b);
                    if (cum / total >= target) {
                        // upper edge of bin b
                        outBounds[t] = hist->GetBinLowEdge(b+1);
                        break;
                    }
                }
                // if not reached, set to max
                if (!std::isfinite(outBounds[t])) outBounds[t] = hist->GetXaxis()->GetXmax();
            } else {
                for (int b = nBins; b >= 1; --b) {
                    cum += hist->GetBinContent(b);
                    if (cum / total >= target) {
                        // lower edge of bin b
                        outBounds[t] = hist->GetBinLowEdge(b);
                        break;
                    }
                }
                if (!std::isfinite(outBounds[t])) outBounds[t] = hist->GetXaxis()->GetXmin();
            }
        }
        // draw lines on current pad
        double ymax = hist->GetMaximum();
        // colors: 1sigma yellow-orangey, 3sigma turquoise/darker blue, 5sigma purple
        int colors[3] = {kOrange+1, kAzure+2, kMagenta+2};
        for (int t = 0; t < 3; ++t) {
            if (!std::isfinite(outBounds[t])) continue;
            TLine *line = nullptr;
            if (integrateLeftToRight) {
                // region is [xmin, outBounds[t]] -> draw vertical at outBounds[t]
                line = new TLine(outBounds[t], 0, outBounds[t], ymax);
            } else {
                // region is [outBounds[t], xmax] -> draw vertical at outBounds[t]
                line = new TLine(outBounds[t], 0, outBounds[t], ymax);
            }
            line->SetLineColor(colors[t]);
            line->SetLineStyle(2);
            line->SetLineWidth(2);
            line->Draw("same");
        }
    };
    // (no one-sided bounds for unweighted histograms)
    c10->Update();
    c10->SaveAs((outputDir + "jarlskog_unweighted.png").c_str());

    TCanvas *c10_log = new TCanvas("c10_log", "Jarlskog Invariant Log Scale", 900, 700);
    c10_log->SetLogy();
    TH1F *h1d_jarlskog_unweighted_log = (TH1F*)h1d_jarlskog_unweighted->Clone("h1d_jarlskog_unweighted_log");
    h1d_jarlskog_unweighted_log->SetTitle("Jarlskog invariant (Log Scale);Jarlskog invariant;Entries");
    h1d_jarlskog_unweighted_log->Draw("HIST");
    h1d_jarlskog_unweighted_log->Draw("E SAME");
    drawAsimovLineWithLegend(asimovJarlskog);
    // (no one-sided bounds for unweighted histograms)
    c10_log->Update();
    c10_log->SaveAs((outputDir + "jarlskog_unweighted_log.png").c_str());

    TCanvas *c10_weighted = new TCanvas("c10_weighted", "Jarlskog Invariant Weighted", 900, 700);
    TH1F *h1d_jarlskog_weighted = new TH1F("h1d_jarlskog_weighted", "Jarlskog invariant;Jarlskog invariant;Weighted Entries",
                                           nBinsJarlskog, jarlskogMin, jarlskogMax);
    if (hasUmbrellaWeight) {
        tree->Draw("TMath::Sqrt(TMath::Max(0.,sin2th_12)*TMath::Max(0.,1.-sin2th_12)*TMath::Max(0.,sin2th_23)*TMath::Max(0.,1.-sin2th_23))*TMath::Max(0.,1.-sin2th_13)*TMath::Sqrt(TMath::Max(0.,sin2th_13))*TMath::Sin(delta_cp)>>h1d_jarlskog_weighted",
                   "umbrella_weight", "");
    } else {
        tree->Draw("TMath::Sqrt(TMath::Max(0.,sin2th_12)*TMath::Max(0.,1.-sin2th_12)*TMath::Max(0.,sin2th_23)*TMath::Max(0.,1.-sin2th_23))*TMath::Max(0.,1.-sin2th_13)*TMath::Sqrt(TMath::Max(0.,sin2th_13))*TMath::Sin(delta_cp)>>h1d_jarlskog_weighted",
                   "", "");
        std::cout << "WARNING:::::No umbrella weights found, plotting unweighted Jarlskog histogram." << std::endl;
    }
    normalizeHist(h1d_jarlskog_weighted);
    h1d_jarlskog_weighted->SetStats(0);
    h1d_jarlskog_weighted->Sumw2();
    h1d_jarlskog_weighted->Draw("HIST");
    h1d_jarlskog_weighted->Draw("E SAME");
    drawAsimovLineWithLegend(asimovJarlskog);
    // One-sided confidence intervals for Jarlskog (weighted)
    double targets[3] = {sigma1, sigma3, sigma5};
    double bounds[3];
    computeOneSidedBounds(h1d_jarlskog_weighted, targets, bounds);
    // create legend entries for one-sided CI lines (weighted)
    double ymax_w = h1d_jarlskog_weighted->GetMaximum();
    TLine *legLine1 = nullptr;
    TLine *legLine3 = nullptr;
    TLine *legLine5 = nullptr;
    if (std::isfinite(bounds[0])) legLine1 = new TLine(bounds[0], 0, bounds[0], ymax_w);
    if (std::isfinite(bounds[1])) legLine3 = new TLine(bounds[1], 0, bounds[1], ymax_w);
    if (std::isfinite(bounds[2])) legLine5 = new TLine(bounds[2], 0, bounds[2], ymax_w);
    if (legLine1) { legLine1->SetLineColor(kOrange+1); legLine1->SetLineStyle(2); legLine1->SetLineWidth(2); }
    if (legLine3) { legLine3->SetLineColor(kAzure+2);  legLine3->SetLineStyle(2); legLine3->SetLineWidth(2); }
    if (legLine5) { legLine5->SetLineColor(kMagenta+2); legLine5->SetLineStyle(2); legLine5->SetLineWidth(2); }
    TLegend *legCI = new TLegend(0.15, 0.60, 0.50, 0.74);
    legCI->SetBorderSize(0);
    legCI->SetFillStyle(0);
    if (legLine1) legCI->AddEntry(legLine1, "1 #sigma (one-sided)", "l");
    if (legLine3) legCI->AddEntry(legLine3, "3 #sigma (one-sided)", "l");
    if (legLine5) legCI->AddEntry(legLine5, "5 #sigma (one-sided)", "l");
    legCI->Draw();
    c10_weighted->Update();
    c10_weighted->SaveAs((outputDir + "jarlskog_weighted.png").c_str());

    TCanvas *c10_weighted_log = new TCanvas("c10_weighted_log", "Jarlskog Invariant Weighted Log Scale", 900, 700);
    c10_weighted_log->SetLogy();
    TH1F *h1d_jarlskog_weighted_log = (TH1F*)h1d_jarlskog_weighted->Clone("h1d_jarlskog_weighted_log");
    h1d_jarlskog_weighted_log->SetTitle("Jarlskog invariant (Log Scale);Jarlskog invariant;Weighted Entries");
    h1d_jarlskog_weighted_log->Draw("HIST");
    h1d_jarlskog_weighted_log->Draw("E SAME");
    drawAsimovLineWithLegend(asimovJarlskog);
    // Draw previously computed one-sided bounds on the log canvas as well
    for (int t = 0; t < 3; ++t) {
        if (!std::isfinite(bounds[t])) continue;
        double ymax = h1d_jarlskog_weighted_log->GetMaximum();
        int colors_log[3] = {kOrange+1, kAzure+2, kMagenta+2};
        TLine *line = new TLine(bounds[t], 0, bounds[t], ymax);
        line->SetLineColor(colors_log[t]);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");
    }
    // legend for log canvas
    TLegend *legCI_log = new TLegend(0.15, 0.60, 0.50, 0.74);
    legCI_log->SetBorderSize(0);
    legCI_log->SetFillStyle(0);
    if (legLine1) legCI_log->AddEntry(legLine1, "1 #sigma (one-sided)", "l");
    if (legLine3) legCI_log->AddEntry(legLine3, "3 #sigma (one-sided)", "l");
    if (legLine5) legCI_log->AddEntry(legLine5, "5 #sigma (one-sided)", "l");
    legCI_log->Draw();
    c10_weighted_log->Update();
    c10_weighted_log->SaveAs((outputDir + "jarlskog_weighted_log.png").c_str());

    struct TriangleVar {
        std::string name;
        std::string title;
        std::string branch;
        int bins;
        double min;
        double max;
        double asimov;
    };

    std::vector<TriangleVar> triangleVars = {
        {"delta_cp", "#delta_{CP}", "delta_cp", nBins_dcp, -3.1415, 3.1415, asimovDcp},
        {"sin2th_23", "sin^{2}(#theta_{23})", "sin2th_23", nBins_sin2th23, sin2th23_min, sin2th23_max, asimovSin2th23},
        {"sin2th_13", "sin^{2}(#theta_{13})", "sin2th_13", nBins_sin2th13, sin2th13_min, sin2th13_max, asimovSin2th13},
        {"sin2th_12", "sin^{2}(#theta_{12})", "sin2th_12", nBins_sin2th12, sin2th12_min, sin2th12_max, asimovSin2th12},
        {"delm2_23", "#Delta m_{23}^{2} (eV^{2})", "delm2_23", nBins_delm23, delm23_io_min, delm23_no_max, asimovDelm23},
        {"delm2_12", "#Delta m_{12}^{2} (eV^{2})", "delm2_12", nBins_delm12, delm12_min, delm12_max, asimovDelm12}
    };

    auto drawAsimovLineNoLegend = [&](double xValue) {
        if (!isAsimovChain || !std::isfinite(xValue) || !gPad) return;
        gPad->Update();
        const double yMin = gPad->GetUymin();
        const double yMax = gPad->GetUymax();
        TLine *asimovLine = new TLine(xValue, yMin, xValue, yMax);
        asimovLine->SetLineColor(kBlack);
        asimovLine->SetLineStyle(3);
        asimovLine->SetLineWidth(2);
        asimovLine->Draw("SAME");
    };

    auto drawAsimovMarkerNoLegend = [&](double xValue, double yValue) {
        if (!isAsimovChain || !std::isfinite(xValue) || !std::isfinite(yValue) || !gPad) return;
        gPad->Update();
        const double xMin = gPad->GetUxmin();
        const double xMax = gPad->GetUxmax();
        const double yMin = gPad->GetUymin();
        const double yMax = gPad->GetUymax();
        const double dx = 0.012 * (xMax - xMin);
        const double dy = 0.012 * (yMax - yMin);

        TLine *l1 = new TLine(xValue - dx, yValue - dy, xValue + dx, yValue + dy);
        l1->SetLineColor(kRed);
        l1->SetLineWidth(2);
        l1->Draw("SAME");

        TLine *l2 = new TLine(xValue - dx, yValue + dy, xValue + dx, yValue - dy);
        l2->SetLineColor(kRed);
        l2->SetLineWidth(2);
        l2->Draw("SAME");
    };

    auto drawTrianglePlot = [&](const std::string &weightExpr, const std::string &suffix, int mo = 0) {
        if (mo == 1) {
            std::cout << "Drawing triangle plot for IO" << std::endl;
            triangleVars[4].min = delm23_no_min;
            triangleVars[4].max = delm23_no_max;
        } else if (mo == -1) {
            std::cout << "Drawing triangle plot for NO" << std::endl;
            triangleVars[4].min = delm23_io_min;
            triangleVars[4].max = delm23_io_max;
        } else {
            std::cout << "Drawing triangle plot for both MO" << std::endl;
        }
        
        const int nVars = static_cast<int>(triangleVars.size());
        TCanvas *cTri = new TCanvas(("c_triangle" + suffix).c_str(), "Triangle Plot", 1400, 1400);
        cTri->SetFillStyle(0);
        cTri->cd();

        // Use MCMCProcessor-style pad math for non-uniform margins
        const double TPm[4] = {.07, .07, .05, .05};
        const double Pm[2] = {.2, .1};

        const double TPw = 1.0 - TPm[0] - TPm[2];
        const double a_x = (Pm[0] * TPw) / (1.0 * nVars + Pm[0] * (1.0 - 1.0 * nVars));
        const double b_x = (TPw - a_x) / (1.0 * nVars);

        std::vector<double> X_Min(nVars), X_Max(nVars);
        X_Min[0] = TPm[0];
        X_Max[0] = X_Min[0] + a_x + b_x;
        for (int i = 1; i < nVars; ++i) {
            X_Min[i] = X_Max[i - 1];
            X_Max[i] = X_Min[i] + b_x;
        }

        const double TPh = 1.0 - TPm[1] - TPm[3];
        const double a_y = (Pm[1] * TPh) / (1.0 * nVars + Pm[1] * (1.0 - 1.0 * nVars));
        const double b_y = (TPh - a_y) / (1.0 * nVars);

        std::vector<double> Y_Min(nVars), Y_Max(nVars);
        Y_Min[nVars - 1] = TPm[1];
        Y_Max[nVars - 1] = Y_Min[nVars - 1] + a_y + b_y;
        for (int i = nVars - 2; i >= 0; --i) {
            Y_Min[i] = Y_Max[i + 1];
            Y_Max[i] = Y_Min[i] + b_y;
        }

        for (int row = 0; row < nVars; ++row) {
            for (int col = 0; col <= row; ++col) {
                const double x1 = X_Min[col];
                const double x2 = X_Max[col];
                const double y1 = Y_Min[row];
                const double y2 = Y_Max[row];

                TPad *pad = new TPad(Form("tri_pad_%s_%d_%d", suffix.c_str(), row, col), "", x1, y1, x2, y2);
                pad->SetFillStyle(0);
                pad->SetRightMargin(0.02);
                pad->SetTopMargin(0.02);
                pad->SetLeftMargin(col == 0 ? 0.26 : 0.06);
                pad->SetBottomMargin(row == nVars - 1 ? 0.18 : 0.06);
                pad->Draw();
                pad->cd();

                const TriangleVar &xVar = triangleVars[col];
                const TriangleVar &yVar = triangleVars[row];

                if (row == col) {
                    TH1F *h1 = new TH1F(Form("h_tri_1d_%s_%s", suffix.c_str(), xVar.name.c_str()),
                                        "", xVar.bins, xVar.min, xVar.max);
                    std::string drawExpr = xVar.branch + ">>" + h1->GetName();
                    tree->Draw(drawExpr.c_str(), weightExpr.c_str(), "goff");
                    normalizeHist(h1);
                    h1->SetStats(0);
                    h1->GetXaxis()->SetTitle(row == nVars - 1 ? xVar.title.c_str() : "");
                    h1->GetYaxis()->SetTitle(col == 0 ? "Entries" : "");
                    h1->GetXaxis()->SetLabelSize(row == nVars - 1 ? 0.07 : 0.0);
                    h1->GetYaxis()->SetLabelSize(col == 0 ? 0.07 : 0.0);
                    h1->GetXaxis()->SetTitleSize(row == nVars - 1 ? 0.08 : 0.0);
                    h1->GetYaxis()->SetTitleSize(col == 0 ? 0.08 : 0.0);
                    h1->Draw("HIST");
                    h1->Draw("E SAME");
                    drawAsimovLineNoLegend(xVar.asimov);
                } else {
                    TH2F *h2 = new TH2F(Form("h_tri_2d_%s_%s_%s", suffix.c_str(), yVar.name.c_str(), xVar.name.c_str()),
                                        "", xVar.bins, xVar.min, xVar.max, yVar.bins, yVar.min, yVar.max);
                    std::string drawExpr = yVar.branch + ":" + xVar.branch + ">>" + h2->GetName();
                    tree->Draw(drawExpr.c_str(), weightExpr.c_str(), "goff");
                    normalizeHist(h2);
                    h2->SetStats(0);
                    h2->GetXaxis()->SetTitle(row == nVars - 1 ? xVar.title.c_str() : "");
                    h2->GetYaxis()->SetTitle(col == 0 ? yVar.title.c_str() : "");
                    h2->GetXaxis()->SetLabelSize(row == nVars - 1 ? 0.07 : 0.0);
                    h2->GetYaxis()->SetLabelSize(col == 0 ? 0.07 : 0.0);
                    h2->GetXaxis()->SetTitleSize(row == nVars - 1 ? 0.08 : 0.0);
                    h2->GetYaxis()->SetTitleSize(col == 0 ? 0.08 : 0.0);
                    h2->Draw("COLZ");
                    drawAsimovMarkerNoLegend(xVar.asimov, yVar.asimov);
                }
                // return to canvas for next pad
                cTri->cd();
            }
        }

        cTri->Update();
        cTri->SaveAs((outputDir + "triangle" + suffix + ".png").c_str());
    };

    drawTrianglePlot("", "_unweighted", 0);
    if (hasUmbrellaWeight) {
        drawTrianglePlot("umbrella_weight", "_weighted", 0);
    } else {
        std::cout << "WARNING:::::No umbrella weights found, skipping weighted triangle plot." << std::endl;
    }

    drawTrianglePlot("(delm2_23 < 0)","_io_unweighted",-1);
    drawTrianglePlot("(delm2_23 > 0)","_no_unweighted",1);
    if (hasUmbrellaWeight) {
        drawTrianglePlot("(delm2_23 < 0) * umbrella_weight","_io_weighted",-1);
        drawTrianglePlot("(delm2_23 > 0) * umbrella_weight","_no_weighted",1);
    } else {
        std::cout << "WARNING:::::No umbrella weights found, skipping weighted triangle plots marginalized by ordering." << std::endl;
    }

}
