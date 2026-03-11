// C++
#include <algorithm>
#include <iomanip>

// MACH3 PLOTTING
#include "PlottingUtils/PlottingManager.h"

//this file has lots of usage of the ROOT plotting interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

/// @file PlotLLH.cpp
/// @ingroup MaCh3Plotting
/// @author Ewan Miller

// some options for the plots
double ratioPlotSplit;
double yTitleOffset;
double sampleLabelThreshold;
int lineWidth;
bool totalOnSplitPlots;
bool sameAxis;

double ratioLabelScaling;

/// @warning KS: keep raw pointer or ensure manual delete of PlotMan. If spdlog in automatically deleted before PlotMan then destructor has some spdlog and this could cause segfault
MaCh3Plotting::PlottingManager *PlotMan;

/// @brief TPad is SUPER FRAGILE, it is safer to just make raw pointer, while ROOT behave weirdly with smart pointers
void SetTPads(TPad*& LLHPad, TPad*& ratioPad)
{
  int originalErrorLevel = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  if (PlotMan->getPlotRatios())
  {
    LLHPad = new TPad("LLHPad", "LLHPad", 0.0, ratioPlotSplit, 1.0, 1.0);
    LLHPad->SetBottomMargin(0.0);
    ratioPad = new TPad("ratioPad", "ratioPad", 0.0, 0.0, 1.0, ratioPlotSplit);
    ratioPad->SetTopMargin(0.0);
    ratioPad->SetBottomMargin(0.3);
  } else {
    LLHPad = new TPad("AllSampPad", "AllSampPad", 0.0, 0.0, 1.0, 1.0);
    ratioPad = new TPad("AllSampRatioPad", "AllSampRatioPad", 0.0, 0.0, 0.0, 0.0);
  }
  LLHPad->AppendPad();
  ratioPad->AppendPad();
  gErrorIgnoreLevel = originalErrorLevel;
}

void getSplitSampleStack(int fileIdx, std::string parameterName, TH1D LLH_allSams,
                         std::vector<float> &cumSums, std::vector<bool> &drawLabel,
                         THStack *sampleStack, TLegend *splitSamplesLegend,
                         float baselineLLH_main = 0.00001) 
  {
  std::vector<std::string> sampNames = PlotMan->input().getTaggedSamples(PlotMan->getOption<std::vector<std::string>>("sampleTags"));
  size_t nSamples = sampNames.size();

  cumSums.resize(nSamples);
  drawLabel.resize(nSamples);

  float LLH_main_integ = LLH_allSams.Integral();
  float cumSum = 0.0;
  int nBins = LLH_allSams.GetNbinsX();

  MACH3LOG_DEBUG("getting split sample THStack for {} known samples", nSamples);
  for (uint i = 0; i < nSamples; i++)
  {
    std::string sampName = sampNames[i];
    
    MACH3LOG_DEBUG("  on sample {}/{}: {}", i, nSamples, sampName);

    TH1D *LLH_indivSam =
        new TH1D(PlotMan->input().getSampleSpecificLLHScan_TH1D(fileIdx, parameterName, sampName));
    LLH_indivSam->SetName(Form("%i_%s_%s", fileIdx, parameterName.c_str(), sampName.c_str()));
    LLH_indivSam->SetBit(kCanDelete);

    // the hist for this sample didn't exist so move on
    if (LLH_indivSam->GetNbinsX() == 1)
    {
      MACH3LOG_DEBUG("    sample hist had only 1 bin - assuming it doesn't exist");
      delete LLH_indivSam;
      drawLabel[i] = false;
      cumSums[i] = -999.9;
      continue;
    }

    LLH_indivSam->SetStats(0);
    LLH_indivSam->SetLineColor(TColor::GetColorPalette(
        floor(static_cast<float>(i) * TColor::GetNumberOfColors() / static_cast<float>(nSamples))));
    LLH_indivSam->SetFillColor(TColor::GetColorPalette(
        floor(static_cast<float>(i) * TColor::GetNumberOfColors() / static_cast<float>(nSamples))));
    sampleStack->Add(LLH_indivSam);
    splitSamplesLegend->AddEntry(LLH_indivSam, PlotMan->style().prettifySampleName(sampName).c_str(), "lf");

    float lastBinLLH = LLH_indivSam->GetBinContent(nBins);
    cumSum += lastBinLLH;

    MACH3LOG_DEBUG("    Last bin LLH = {} :: cumulative LLH = {}", lastBinLLH, cumSum);
    MACH3LOG_DEBUG("    LLH fraction = {} / {} = {}", LLH_indivSam->Integral(), LLH_main_integ, LLH_indivSam->Integral() / LLH_main_integ);

    cumSums[i] = cumSum;
    
    // dont draw a label if the likelihood contribution is less than threshold%
    if ((LLH_indivSam->Integral() / LLH_main_integ > sampleLabelThreshold) &&
        (LLH_indivSam->Integral() / baselineLLH_main > sampleLabelThreshold))
    {
      drawLabel[i] = true;
    } else {
      drawLabel[i] = false;
    }
    MACH3LOG_DEBUG("    drawLabel = {}", drawLabel.back()); 
    MACH3LOG_DEBUG("");
  }
}

// handy dandy helper function for drawing and nicely formatting a stack of ratio plots
void drawRatioStack(THStack *ratioCompStack) {
  // do this so 1.0 is in the middle of the plot vertically
  double stackMax = ratioCompStack->GetMaximum("NOSTACK");
  double stackMin = ratioCompStack->GetMinimum("NOSTACK");

  double stackLim = std::max(std::abs(1.0 - stackMax), std::abs(1.0 - stackMin));

  ratioCompStack->SetMinimum(1.0 - 1.05 * stackLim);
  ratioCompStack->SetMaximum(1.0 + 1.05 * stackLim);

  // draw it
  ratioCompStack->Draw(Form("NOSTACK%s", PlotMan->getDrawOptions().c_str()));

  // make it look a bit nicer
  ratioCompStack->GetXaxis()->SetLabelSize(ratioLabelScaling *
                                           ratioCompStack->GetXaxis()->GetLabelSize());
  ratioCompStack->GetXaxis()->SetTitleSize(ratioLabelScaling *
                                           ratioCompStack->GetXaxis()->GetLabelSize());
  ratioCompStack->GetXaxis()->SetTitle("Parameter Variation");

  ratioCompStack->GetYaxis()->SetLabelSize(ratioLabelScaling *
                                           ratioCompStack->GetYaxis()->GetLabelSize());
  ratioCompStack->GetYaxis()->SetTitleSize(ratioLabelScaling *
                                           ratioCompStack->GetYaxis()->GetLabelSize());
  ratioCompStack->GetYaxis()->SetTitleOffset(yTitleOffset);
  ratioCompStack->GetYaxis()->SetNdivisions(5, 2, 0);
}

void makeLLHScanComparisons(const std::string& paramName,
                            const std::string& LLHType,
                            const std::string& outputFileName,
                            const std::unique_ptr<TCanvas>& canv) {
  // will use these to make comparisons of likelihoods
  auto compStack = std::make_unique<THStack>("LLH_Stack", "");
  auto ratioCompStack = std::make_unique<THStack>("LLH_Ratio_Stack", "");
  auto legend = std::make_unique<TLegend>(0.3, 0.6, 0.7, 0.8);

  TPad* LLHPad = nullptr;
  TPad* ratioPad = nullptr;
  SetTPads(LLHPad, ratioPad);

  // get the sample reweight hist from the main file
  TH1D LLH_main = PlotMan->input().getLLHScan_TH1D(0, paramName, LLHType);
  LLH_main.SetStats(0);

  LLH_main.SetLineColor(kBlack);
  compStack->Add(&LLH_main);
  legend->AddEntry(&LLH_main, PlotMan->getFileLabel(0).c_str(), "l");

  int nBins = LLH_main.GetNbinsX();

  // go through the other files
  for (unsigned int extraFileIdx = 1; extraFileIdx < PlotMan->input().getNInputFiles(); extraFileIdx++)
  {
    int originalErrorLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;

    TH1D *compHist = new TH1D(PlotMan->input().getLLHScan_TH1D(extraFileIdx, paramName, LLHType));
    compHist->SetName(Form("LLHScan_%s_%s_%d", paramName.c_str(), LLHType.c_str(), extraFileIdx));
    compHist->SetBit(kCanDelete); // <- will allow this to be deleted by root once done plotting
    gErrorIgnoreLevel = originalErrorLevel;
    if (compHist->GetNbinsX() == 0)
      continue;

    // make them look different to each other
    compHist->SetLineColor(
        TColor::GetColorPalette(floor(static_cast<float>(extraFileIdx) * TColor::GetNumberOfColors() /
                                      static_cast<float>(PlotMan->input().getNInputFiles()))));
    compHist->SetLineStyle(2 + extraFileIdx % 9);
    compHist->SetLineWidth(lineWidth);

    TH1D *divHist = static_cast<TH1D*>(compHist->Clone(Form("RatioHist_%i", extraFileIdx)));
    divHist->SetBit(kCanDelete);

    if (PlotMan->getPlotRatios())
    {
      // get the ratio hist
      divHist->Divide(compHist, &LLH_main);
      ratioCompStack->Add(divHist);
    }

    // add it to the comparison stack and legend
    compStack->Add(compHist);
    legend->AddEntry(compHist, PlotMan->getFileLabel(extraFileIdx).c_str(), "l");
  }

  // draw the log likelihoods
  canv->cd();
  canv->Draw();
  LLHPad->Draw();
  LLHPad->cd();
  if (PlotMan->getDrawGrid())
    LLHPad->SetGrid();
  compStack->Draw(Form("NOSTACK%s", PlotMan->getDrawOptions().c_str()));
  if (!PlotMan->getPlotRatios())
    compStack->GetXaxis()->SetTitle("Parameter Variation");
  compStack->GetYaxis()->SetTitle(Form("-2LLH_{%s}", LLHType.c_str()));
  compStack->SetTitle(PlotMan->style().prettifyParamName(paramName).c_str());
  compStack->GetYaxis()->SetTitleOffset(yTitleOffset);
  legend->Draw();

  // add the ratio plot if specified
  if (PlotMan->getPlotRatios())
  {
    canv->cd();
    ratioPad->Draw();
    ratioPad->cd();
    if (PlotMan->getDrawGrid())
      ratioPad->SetGrid();

    drawRatioStack(ratioCompStack.get());

    // add horizontal line at 1 for reference
    TLine line = TLine();
    line.SetLineColor(kBlack);
    line.SetLineWidth(lineWidth);
    line.DrawLine(LLH_main.GetBinLowEdge(1), 1.0, LLH_main.GetBinLowEdge(nBins + 1), 1.0);
  }

  // save to the output file
  canv->SaveAs(outputFileName.c_str());
  delete LLHPad;
  delete ratioPad;
}

void makeSplitSampleLLHScanComparisons(const std::string& paramName,
                                       const std::string& outputFileName,
                                       const std::unique_ptr<TCanvas>& canv) {
  MACH3LOG_DEBUG(" Making split sample LLH comparison");
  // split the canvas if plotting ratios
  TPad* LLHPad = nullptr;
  TPad* ratioPad = nullptr;
  SetTPads(LLHPad, ratioPad);
  canv->Clear();
  canv->Draw();

  canv->Divide(PlotMan->getNFiles());

  // get the sample hist from the main file
  TH1D LLH_main = PlotMan->input().getLLHScan_TH1D(0, paramName, "sample");
  if (LLH_main.GetNbinsX() == 1)
  {
    MACH3LOG_DEBUG("  Main LLH had only 1 bin, assuming it doesn't exist");
    return;
  }

  auto baseSplitSamplesStack = std::make_unique<THStack>(
      paramName.c_str(), Form("%s - %s", paramName.c_str(), PlotMan->getFileLabel(0).c_str()));

  auto baseSplitSamplesLegend = std::make_unique<TLegend>(0.37, 0.475, 0.63, 0.9);

  if (PlotMan->getDrawGrid())
    canv->cd(1)->SetGrid();
  else
    canv->cd(1);

  canv->cd(1)->SetLeftMargin(0.15);
  std::vector<float> cumSums;
  std::vector<bool> drawLabel;

  getSplitSampleStack(0, paramName, LLH_main, cumSums, drawLabel, baseSplitSamplesStack.get(),
                      baseSplitSamplesLegend.get());
  baseSplitSamplesStack->Draw(PlotMan->getDrawOptions().c_str());
  // KS: Not sure why but need to plot after it's drawn otherwise ROOT throws segfault...
  baseSplitSamplesStack->GetYaxis()->SetTitle("-2LLH_{sam}");
  if (PlotMan->getPlotRatios() == false) baseSplitSamplesStack->GetXaxis()->SetTitle("Parameter Variation");
  gPad->Modified();
  gPad->Update();

  if (totalOnSplitPlots)
  {
    LLH_main.SetLineWidth(1); // undo SetLineWidth that was done above
    LLH_main.Draw(Form("same%s", PlotMan->getDrawOptions().c_str()));
    baseSplitSamplesLegend->AddEntry(&LLH_main, "All Samples", "l");
  }

  baseSplitSamplesLegend->Draw();

  auto label = std::make_unique<TLatex>();
  // format the label
  label->SetTextAlign(11);
  label->SetTextAngle(-55);
  label->SetTextSize(0.012);

  // need to draw the labels after other stuff or they don't show up
  for (uint i = 0; i < PlotMan->input().getTaggedSamples(PlotMan->getOption<std::vector<std::string>>("sampleTags")).size(); i++)
  {
    MACH3LOG_DEBUG("  Will I draw the label for sample {}??", i);
    std::string sampName = PlotMan->input().getTaggedSamples(PlotMan->getOption<std::vector<std::string>>("sampleTags"))[i];
    if (!drawLabel[i])
    { 
      MACH3LOG_DEBUG("   - Not drawing label");
      continue;
    }
    MACH3LOG_DEBUG("   - Drawing label");

    label->DrawLatex(LLH_main.GetBinLowEdge(LLH_main.GetNbinsX() + 1), cumSums[i],
                     Form("#leftarrow%s", PlotMan->style().prettifySampleName(sampName).c_str()));
    MACH3LOG_DEBUG("  I drew the label!");
  }

  // now we plot the comparison file plots
  for (unsigned int extraFileIdx = 1; extraFileIdx < PlotMan->getNFiles(); extraFileIdx++)
  {
    MACH3LOG_DEBUG("  - Adding plot for additional file {}", extraFileIdx);
    canv->cd(1 + extraFileIdx);

    std::vector<float> extraCumSums;
    std::vector<bool> extraDrawLabel;

    THStack *splitSamplesStack =
        new THStack(paramName.c_str(),
                    Form("%s - %s", paramName.c_str(), PlotMan->getFileLabel(extraFileIdx).c_str()));

    auto splitSamplesLegend = std::make_unique<TLegend>(0.37, 0.475, 0.63, 0.9);
    splitSamplesStack->SetBit(kCanDelete);
    splitSamplesLegend->SetBit(kCanDelete);

    int originalErrorLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kFatal;
    TH1D compLLH_main = PlotMan->input().getLLHScan_TH1D(extraFileIdx, paramName, "sample");
    compLLH_main.SetName((compLLH_main.GetName() + std::string("_file_") + extraFileIdx).Data());
    gErrorIgnoreLevel = originalErrorLevel;
    if (compLLH_main.GetNbinsX() == 1)
    {
      delete splitSamplesStack;
      continue;
    }

    // if on the same y axis, also check that the contribution of the sample compared to the
    // baseline LLH integral is above the threshold otherwise the labels might get very crowded if
    // the comparisson LLH is much smaller than the baseline one
    if (sameAxis)
      getSplitSampleStack(extraFileIdx, paramName, compLLH_main, extraCumSums, extraDrawLabel,
                          splitSamplesStack, splitSamplesLegend.get(), LLH_main.Integral());
    else
      getSplitSampleStack(extraFileIdx, paramName, compLLH_main, extraCumSums, extraDrawLabel,
                          splitSamplesStack, splitSamplesLegend.get());

    // go to the pad for the histograms
    if (PlotMan->getDrawGrid())
      LLHPad->SetGrid();
    LLHPad->Draw();
    LLHPad->cd();

    splitSamplesStack->Draw(PlotMan->getDrawOptions().c_str());
    if (sameAxis)
      splitSamplesStack->SetMaximum(baseSplitSamplesStack->GetMaximum());

    if (totalOnSplitPlots)
    {
      compLLH_main.Draw(Form("same%s", PlotMan->getDrawOptions().c_str()));
      compLLH_main.SetLineColor(kBlack);
      splitSamplesLegend->AddEntry(&compLLH_main, "All Samples", "l");
    }
    splitSamplesLegend->Draw();

    // KS: Not sure why but need to plot after it's drawn otherwise ROOT throws segfault...
    baseSplitSamplesStack->GetYaxis()->SetTitle("-2LLH_{sam}");
    if (PlotMan->getPlotRatios() == false) baseSplitSamplesStack->GetXaxis()->SetTitle("Parameter Variation");
    gPad->Modified();
    gPad->Update();

    // need to draw the labels after other stuff or they dont show up
    for (uint i = 0; i < PlotMan->input().getTaggedSamples(PlotMan->getOption<std::vector<std::string>>("sampleTags")).size(); i++)
    {
      std::string sampName = PlotMan->input().getTaggedSamples(PlotMan->getOption<std::vector<std::string>>("sampleTags"))[i];
      if (!drawLabel[i])
        continue;
      label->DrawLatex(compLLH_main.GetBinLowEdge(compLLH_main.GetNbinsX() + 1), extraCumSums[i],
                       Form("#leftarrow%s", PlotMan->style().prettifySampleName(sampName).c_str()));
    }

    if (PlotMan->getPlotRatios())
    {
      THStack *splitSamplesStackRatios = new THStack(paramName.c_str(), "");

      TList *baselineHistList = baseSplitSamplesStack->GetHists();
      TList *compHistList = splitSamplesStack->GetHists();
      for (uint sampleIdx = 0; sampleIdx < PlotMan->input().getTaggedSamples(PlotMan->getOption<std::vector<std::string>>("sampleTags")).size(); sampleIdx++)
      {
        TH1D *divHist = new TH1D(
            Form("%s_%s_splitDiv_%i", paramName.c_str(), PlotMan->getFileLabel(extraFileIdx).c_str(),
                 sampleIdx),
            Form("%s_%s_splitDiv", paramName.c_str(), PlotMan->getFileLabel(extraFileIdx).c_str()),
            compLLH_main.GetNbinsX(), compLLH_main.GetBinLowEdge(1),
            compLLH_main.GetBinLowEdge(compLLH_main.GetNbinsX() + 1));
        divHist->Divide(static_cast<TH1D*>(compHistList->At(sampleIdx)),
                        static_cast<TH1D*>(baselineHistList->At(sampleIdx)));
        splitSamplesStackRatios->Add(divHist);
        divHist->SetLineColor((static_cast<TH1D*>(compHistList->At(sampleIdx))->GetLineColor()));
      }

      canv->cd(2 + extraFileIdx);
      if (PlotMan->getDrawGrid())
        ratioPad->SetGrid();
      ratioPad->Draw();
      ratioPad->cd();
      drawRatioStack(splitSamplesStackRatios);
    }
  }
  canv->SaveAs(outputFileName.c_str());
  delete LLHPad;
}

int PlotLLH() {
  PlotMan->style().setPalette(PlotMan->getOption<std::string>("colorPalette"));

  // make the canvas and other plotting stuff
  auto canv = std::make_unique<TCanvas>("canv", "", 1024, 1024);
  gStyle->SetOptTitle(2);
  // split up the canvas so can have side by side plots, one for each file
  auto splitSamplesCanv = std::make_unique<TCanvas>("splitSampCanv", "", 4096 * PlotMan->getNFiles(), 4096);

  canv->SaveAs((PlotMan->getOutputName("_Sample") + "[").c_str());
  canv->SaveAs((PlotMan->getOutputName("_Penalty") + "[").c_str());
  canv->SaveAs((PlotMan->getOutputName("_Total") + "[").c_str());

  if (PlotMan->getSplitBySample())
    canv->SaveAs((PlotMan->getOutputName("_bySample") + "[").c_str());

  for( std::string par: PlotMan->getOption<std::vector<std::string>>("parameterTags")) std::cout << par << ", ";
  // loop over the spline parameters
  for (std::string paramName : PlotMan->input().getTaggedParameters(PlotMan->getOption<std::vector<std::string>>("parameterTags")))
  {
    MACH3LOG_DEBUG("working on parameter {}", paramName);
    // ###############################################################
    // First lets do just the straight up likelihoods from all samples
    // ###############################################################
    makeLLHScanComparisons(paramName, "sample", PlotMan->getOutputName("_Sample"), canv);
    makeLLHScanComparisons(paramName, "penalty", PlotMan->getOutputName("_Penalty"), canv);
    makeLLHScanComparisons(paramName, "total", PlotMan->getOutputName("_Total"), canv);
    // #########################################
    // ## now lets make plots split by sample ##
    // #########################################
    if (PlotMan->getSplitBySample())
    {
      makeSplitSampleLLHScanComparisons(paramName, PlotMan->getOutputName("_bySample"), splitSamplesCanv);
    }
  }

  canv->SaveAs((PlotMan->getOutputName("_Sample") + "]").c_str());
  canv->SaveAs((PlotMan->getOutputName("_Penalty") + "]").c_str());
  canv->SaveAs((PlotMan->getOutputName("_Total") + "]").c_str());
  if (PlotMan->getSplitBySample())
    canv->SaveAs((PlotMan->getOutputName("_bySample") + "]").c_str());

  if(PlotMan != nullptr) delete PlotMan;

  return 0;
}

int main(int argc, char **argv) {
  SetMaCh3LoggerFormat();
  MaCh3Utils::MaCh3Welcome();

  PlotMan = new MaCh3Plotting::PlottingManager();
  PlotMan->parseInputs(argc, argv);

  PlotMan->setExec("PlotLLH");

  ratioPlotSplit = PlotMan->getOption<double>("ratioPlotSplit");
  yTitleOffset = PlotMan->getOption<double>("yTitleOffset");
  sampleLabelThreshold = PlotMan->getOption<double>("sampleLabelThreshold");
  lineWidth = PlotMan->getOption<int>("lineWidth");
  totalOnSplitPlots = PlotMan->getOption<bool>("totalOnSplitPlots");
  sameAxis = PlotMan->getOption<bool>("sameAxis");

  // scale the ratio plot labels by this much to make them same size as the normal plot
  ratioLabelScaling = (1.0 / ratioPlotSplit - 1.0);

  if(PlotMan->getPlotRatios() && PlotMan->input().getNInputFiles() == 1){
    MACH3LOG_ERROR("Can't plot ratio with single file...");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  return PlotLLH();
}
