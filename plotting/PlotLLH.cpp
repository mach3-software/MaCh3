// C++
#include <algorithm>
#include <iomanip>

// MACH3 PLOTTING
#include "plottingUtils/plottingManager.h"

//this file has lots of usage of the ROOT plotting interface that only takes floats, turn this warning off for this CU for now
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#pragma GCC diagnostic ignored "-Wconversion"

// some options for the plots
double ratioPlotSplit;
double yTitleOffset;
double sampleLabelThreshold;
int lineWidth;
bool totalOnSplitPlots;
bool sameAxis;

double ratioLabelScaling;

MaCh3Plotting::PlottingManager *man;

void getSplitSampleStack(int fileIdx, std::string parameterName, TH1D LLH_allSams,
                         std::vector<float> &cumSums, std::vector<bool> &drawLabel,
                         THStack *sampleStack, TLegend *splitSamplesLegend,
                         float baselineLLH_main = 0.00001) 
  {
  std::vector<std::string> sampNames = man->input().getTaggedSamples(man->getOption<std::vector<std::string>>("sampleTags"));
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
        new TH1D(man->input().getSampleSpecificLLHScan_TH1D(fileIdx, parameterName, sampName));
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
    splitSamplesLegend->AddEntry(LLH_indivSam, man->style().prettifySampleName(sampName).c_str(),
                                 "lf");

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
    } else
    {
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
  ratioCompStack->Draw(Form("NOSTACK%s", man->getDrawOptions().c_str()));

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

void makeLLHScanComparisons(std::string paramName, std::string LLHType, std::string outputFileName,
                            TCanvas *canv, TPad *LLHPad, TPad *ratioPad) {
  // will use these to make comparisons of likelihoods
  THStack *compStack = new THStack("LLH_Stack", "");
  THStack *ratioCompStack = new THStack("LLH_Ratio_Stack", "");
  auto legend = std::make_unique<TLegend>(0.3, 0.6, 0.7, 0.8);

  // get the sample reweight hist from the main file
  TH1D LLH_main = man->input().getLLHScan_TH1D(0, paramName, LLHType);
  LLH_main.SetStats(0);

  LLH_main.SetLineColor(kBlack);
  compStack->Add(&LLH_main);
  legend->AddEntry(&LLH_main, man->getFileLabel(0).c_str(), "l");

  int nBins = LLH_main.GetNbinsX();

  // go through the other files
  for (unsigned int extraFileIdx = 1; extraFileIdx < man->input().getNInputFiles(); extraFileIdx++)
  {
    TH1D *compHist = new TH1D(man->input().getLLHScan_TH1D(extraFileIdx, paramName, LLHType));
    compHist->SetBit(kCanDelete); // <- will allow this to be deleted by root once done plotting
    if (compHist->GetNbinsX() == 0)
      continue;

    // make them look different to each other
    compHist->SetLineColor(
        TColor::GetColorPalette(floor(static_cast<float>(extraFileIdx) * TColor::GetNumberOfColors() /
                                      static_cast<float>(man->input().getNInputFiles()))));
    compHist->SetLineStyle(2 + extraFileIdx % 9);
    compHist->SetLineWidth(lineWidth);

    TH1D *divHist = static_cast<TH1D*>(compHist->Clone(Form("RatioHist_%i", extraFileIdx)));
    divHist->SetBit(kCanDelete);

    if (man->getPlotRatios())
    {
      // get the ratio hist
      divHist->Divide(compHist, &LLH_main);
      ratioCompStack->Add(divHist);
    }

    // add it to the comparisson hstack and legend
    compStack->Add(compHist);
    legend->AddEntry(compHist, man->getFileLabel(extraFileIdx).c_str(), "l");
  }

  // draw the log likelihoods
  canv->cd();
  canv->Draw();
  LLHPad->Draw();
  LLHPad->cd();
  if (man->getDrawGrid())
    LLHPad->SetGrid();
  compStack->Draw(Form("NOSTACK%s", man->getDrawOptions().c_str()));
  if (!man->getPlotRatios())
    compStack->GetXaxis()->SetTitle("Parameter Variation");
  compStack->GetYaxis()->SetTitle(Form("-2LLH_{%s}", LLHType.c_str()));
  compStack->SetTitle(man->style().prettifyParamName(paramName).c_str());
  compStack->GetYaxis()->SetTitleOffset(yTitleOffset);
  legend->Draw();

  // add the ratio plot if specified
  if (man->getPlotRatios())
  {
    canv->cd();
    ratioPad->Draw();
    ratioPad->cd();
    if (man->getDrawGrid())
      ratioPad->SetGrid();

    drawRatioStack(ratioCompStack);

    // add horizontal line at 1 for reference
    TLine line = TLine();
    line.SetLineColor(kBlack);
    line.SetLineWidth(lineWidth);
    line.DrawLine(LLH_main.GetBinLowEdge(1), 1.0, LLH_main.GetBinLowEdge(nBins + 1), 1.0);
  }

  // save to the output file
  canv->SaveAs(outputFileName.c_str());

  // will use these to make comparisons of likelihoods
  delete compStack;
  delete ratioCompStack;
}

void makeSplitSampleLLHScanComparisons(std::string paramName, std::string outputFileName,
                                       TCanvas *canv, TPad *LLHPad, TPad *ratioPad) {
  MACH3LOG_DEBUG(" Making split sample LLH comparison");
  canv->Clear();
  canv->Draw();

  canv->Divide(man->getNFiles());

  // get the sample hist from the main file
  TH1D LLH_main = man->input().getLLHScan_TH1D(0, paramName, "sample");
  if (LLH_main.GetNbinsX() == 1)
  {
    MACH3LOG_DEBUG("  Main LLH had only 1 bin, assuming it doesn't exist");
    return;
  }

  THStack *baseSplitSamplesStack = new THStack(
      paramName.c_str(), Form("%s - %s", paramName.c_str(), man->getFileLabel(0).c_str()));
  TLegend *baseSplitSamplesLegend = new TLegend(0.37, 0.475, 0.63, 0.9);

  if (man->getDrawGrid())
    canv->cd(1)->SetGrid();
  else
    canv->cd(1);

  std::vector<float> cumSums;
  std::vector<bool> drawLabel;

  getSplitSampleStack(0, paramName, LLH_main, cumSums, drawLabel, baseSplitSamplesStack,
                      baseSplitSamplesLegend);

  baseSplitSamplesStack->Draw(man->getDrawOptions().c_str());

  if (totalOnSplitPlots)
  {
    LLH_main.SetLineWidth(1); // undo SetLineWidth that was done above
    LLH_main.Draw(Form("same%s", man->getDrawOptions().c_str()));
    baseSplitSamplesLegend->AddEntry(&LLH_main, "All Samples", "l");
  }

  baseSplitSamplesLegend->Draw();

  TLatex *label = new TLatex;
  // format the label
  label->SetTextAlign(11);
  label->SetTextAngle(-55);
  label->SetTextSize(0.012);

  // need to draw the labels after other stuff or they dont show up
  for (uint i = 0; i < man->input().getTaggedSamples(man->getOption<std::vector<std::string>>("sampleTags")).size(); i++)
  {
    MACH3LOG_DEBUG("  Will I draw the label for sample {}??", i);
    std::string sampName = man->input().getTaggedSamples(man->getOption<std::vector<std::string>>("sampleTags"))[i];
    if (!drawLabel[i])
    { 
      MACH3LOG_DEBUG("   - Not drawing label");
      continue;
    }
    MACH3LOG_DEBUG("   - Drawing label");

    label->DrawLatex(LLH_main.GetBinLowEdge(LLH_main.GetNbinsX() + 1), cumSums[i],
                     Form("#leftarrow%s", man->style().prettifySampleName(sampName).c_str()));
    MACH3LOG_DEBUG("  I drew the label!");
  }

  // now we plot the comparisson file plots
  for (unsigned int extraFileIdx = 1; extraFileIdx < man->getNFiles(); extraFileIdx++)
  {
    MACH3LOG_DEBUG("  - Adding plot for additional file {}", extraFileIdx);
    canv->cd(1 + extraFileIdx);

    std::vector<float> extraCumSums;
    std::vector<bool> extraDrawLabel;

    THStack *splitSamplesStack =
        new THStack(paramName.c_str(),
                    Form("%s - %s", paramName.c_str(), man->getFileLabel(extraFileIdx).c_str()));
    TLegend *splitSamplesLegend = new TLegend(0.37, 0.475, 0.63, 0.9);

    splitSamplesStack->SetBit(kCanDelete);
    splitSamplesLegend->SetBit(kCanDelete);

    TH1D compLLH_main = man->input().getLLHScan_TH1D(extraFileIdx, paramName, "sample");
    if (compLLH_main.GetNbinsX() == 1)
    {
      delete splitSamplesStack;
      delete splitSamplesLegend;
      continue;
    }

    // if on the same y axis, also check that the contribution of the sample compared to the
    // baseline LLH integral is above the threshold otherwise the labels might get very crowded if
    // the comparisson LLH is much smaller than the baseline one
    if (sameAxis)
      getSplitSampleStack(extraFileIdx, paramName, compLLH_main, extraCumSums, extraDrawLabel,
                          splitSamplesStack, splitSamplesLegend, LLH_main.Integral());
    else
      getSplitSampleStack(extraFileIdx, paramName, compLLH_main, extraCumSums, extraDrawLabel,
                          splitSamplesStack, splitSamplesLegend);

    // go to the pad for the histograms
    if (man->getDrawGrid())
      LLHPad->SetGrid();
    LLHPad->Draw();
    LLHPad->cd();

    splitSamplesStack->Draw(man->getDrawOptions().c_str());
    splitSamplesStack->GetXaxis()->SetTitle("Parameter Variation");
    splitSamplesStack->GetYaxis()->SetTitle("-2LLH_{sam}");
    if (sameAxis)
      splitSamplesStack->SetMaximum(baseSplitSamplesStack->GetMaximum());

    if (totalOnSplitPlots)
    {
      compLLH_main.Draw(Form("same%s", man->getDrawOptions().c_str()));
      compLLH_main.SetLineColor(kBlack);
      splitSamplesLegend->AddEntry(&compLLH_main, "All Samples", "l");
    }
    splitSamplesLegend->Draw();

    // need to draw the labels after other stuff or they dont show up
    for (uint i = 0; i < man->input().getTaggedSamples(man->getOption<std::vector<std::string>>("sampleTags")).size(); i++)
    {
      std::string sampName = man->input().getTaggedSamples(man->getOption<std::vector<std::string>>("sampleTags"))[i];
      if (!drawLabel[i])
        continue;
      label->DrawLatex(compLLH_main.GetBinLowEdge(compLLH_main.GetNbinsX() + 1), extraCumSums[i],
                       Form("#leftarrow%s", man->style().prettifySampleName(sampName).c_str()));
    }

    if (man->getPlotRatios())
    {
      THStack *splitSamplesStackRatios = new THStack(paramName.c_str(), "");

      TList *baselineHistList = baseSplitSamplesStack->GetHists();
      TList *compHistList = splitSamplesStack->GetHists();

      for (uint sampleIdx = 0; sampleIdx < man->input().getTaggedSamples(man->getOption<std::vector<std::string>>("sampleTags")).size(); sampleIdx++)
      {
        TH1D *divHist = new TH1D(
            Form("%s_%s_splitDiv_%i", paramName.c_str(), man->getFileLabel(extraFileIdx).c_str(),
                 sampleIdx),
            Form("%s_%s_splitDiv", paramName.c_str(), man->getFileLabel(extraFileIdx).c_str()),
            compLLH_main.GetNbinsX(), compLLH_main.GetBinLowEdge(1),
            compLLH_main.GetBinLowEdge(compLLH_main.GetNbinsX() + 1));
        divHist->Divide(static_cast<TH1D*>(compHistList->At(sampleIdx)),
                        static_cast<TH1D*>(baselineHistList->At(sampleIdx)));
        splitSamplesStackRatios->Add(divHist);
        divHist->SetLineColor((static_cast<TH1D*>(compHistList->At(sampleIdx))->GetLineColor()));
      }

      canv->cd(2 + extraFileIdx);
      if (man->getDrawGrid())
        ratioPad->SetGrid();
      ratioPad->Draw();
      ratioPad->cd();

      drawRatioStack(splitSamplesStackRatios);
    }
  }

  canv->SaveAs(outputFileName.c_str());

  delete baseSplitSamplesStack;
  delete baseSplitSamplesLegend;
}

int PlotLLH() {

  man->style().setPalette(man->getOption<std::string>("colorPalette"));

  // make the canvas and other plotting stuff
  TCanvas *canv = new TCanvas("canv", "", 1024, 1024);

  // split up the canvas so can have side by side plots, one for each file
  TCanvas *splitSamplesCanv = new TCanvas("splitSampCanv", "", 4096 * man->getNFiles(), 4096);

  // split the canvas if plotting ratios
  TPad *LLHPad, *ratioPad;
  if (man->getPlotRatios())
  {
    LLHPad = new TPad("LLHPad", "LLHPad", 0.0, ratioPlotSplit, 1.0, 1.0);
    LLHPad->SetBottomMargin(0.0);
    ratioPad = new TPad("ratioPad", "ratioPad", 0.0, 0.0, 1.0, ratioPlotSplit);
    ratioPad->SetTopMargin(0.0);
    ratioPad->SetBottomMargin(0.3);
  } else
  {
    LLHPad = new TPad("AllSampPad", "AllSampPad", 0.0, 0.0, 1.0, 1.0);
    ratioPad = new TPad("AllSampRatioPad", "AllSampRatioPad", 0.0, 0.0, 0.0, 0.0);
  }

  canv->SaveAs((man->getOutputName("_Sample") + "[").c_str());
  canv->SaveAs((man->getOutputName("_Penalty") + "[").c_str());
  canv->SaveAs((man->getOutputName("_Total") + "[").c_str());

  if (man->getSplitBySample())
    canv->SaveAs((man->getOutputName("_bySample") + "[").c_str());

  for( std::string par: man->getOption<std::vector<std::string>>("parameterTags")) std::cout << par << ", ";
  // loop over the spline parameters
  for (std::string paramName : man->input().getTaggedParameters(man->getOption<std::vector<std::string>>("parameterTags")))
  {
    MACH3LOG_DEBUG("working on parameter {}", paramName);
    // ###############################################################
    // First lets do just the straight up likelihoods from all samples
    // ###############################################################

    makeLLHScanComparisons(paramName, "sample", man->getOutputName("_Sample"), canv, LLHPad,
                           ratioPad);
    makeLLHScanComparisons(paramName, "penalty", man->getOutputName("_Penalty"), canv, LLHPad,
                           ratioPad);
    makeLLHScanComparisons(paramName, "total", man->getOutputName("_Total"), canv, LLHPad,
                           ratioPad);

    // #########################################
    // ## now lets make plots split by sample ##
    // #########################################
    if (man->getSplitBySample())
    {
      makeSplitSampleLLHScanComparisons(paramName, man->getOutputName("_bySample"),
                                        splitSamplesCanv, LLHPad, ratioPad);
    }
  }

  canv->SaveAs((man->getOutputName("_Sample") + "]").c_str());
  canv->SaveAs((man->getOutputName("_Penalty") + "]").c_str());
  canv->SaveAs((man->getOutputName("_Total") + "]").c_str());
  if (man->getSplitBySample())
    canv->SaveAs((man->getOutputName("_bySample") + "]").c_str());

  delete canv;

  return 0;
}

int main(int argc, char **argv) {
  SetMaCh3LoggerFormat();

  man = new MaCh3Plotting::PlottingManager();
  man->parseInputs(argc, argv);

  man->setExec("PlotLLH");

  ratioPlotSplit = man->getOption<double>("ratioPlotSplit");
  yTitleOffset = man->getOption<double>("yTitleOffset");
  sampleLabelThreshold = man->getOption<double>("sampleLabelThreshold");
  lineWidth = man->getOption<int>("lineWidth");
  totalOnSplitPlots = man->getOption<bool>("totalOnSplitPlots");
  sameAxis = man->getOption<bool>("sameAxis");

  // scale the ratio plot labels by this much to make them same size as the normal plot
  ratioLabelScaling = (1.0 / ratioPlotSplit - 1.0);

  return PlotLLH();
}
