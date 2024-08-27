#include "styleManager.h"

namespace MaCh3Plotting {
StyleManager::StyleManager(std::string styleConfigName) {
  _styleConfig = YAML::LoadFile(styleConfigName);
}

std::string StyleManager::prettifyName(const std::string &origName, const std::string &nameType) const {
  std::string prettyName = origName;

  YAML::Node prettyNames = _styleConfig["PrettyNames"][nameType];

  if (prettyNames[origName])
  {
    prettyName = prettyNames[origName].as<std::string>();
  }

  return prettyName;
}

void StyleManager::setPalette(int rootPlotStyle) const {
  // set the colour palette to one of the root palettes
  gStyle->SetPalette(rootPlotStyle);
}

void StyleManager::setPalette(std::string configStyleName) const {
  // set the colour palette to one of the palettes defined in PlottingConfig.yaml

  // get the definition of the provided style from the config file
  YAML::Node palettes = _styleConfig["ColorPallettes"];
  std::vector<std::vector<double>> paletteDef =
      palettes[configStyleName].as<std::vector<std::vector<double>>>();

  const Int_t NCont = (Int_t)(paletteDef[0][0]);

  std::vector<double> stopVec = paletteDef[1];
  std::vector<double> redsVec = paletteDef[2];
  std::vector<double> greensVec = paletteDef[3];
  std::vector<double> bluesVec = paletteDef[4];

  // get the number of colour stops and check all vectors are same size
  const Int_t NRGBs = stopVec.size();
  if ((Int_t)redsVec.size() != NRGBs || (Int_t)greensVec.size() != NRGBs ||
      (Int_t)bluesVec.size() != NRGBs)
  {
    MACH3LOG_ERROR("invalid colour palettet defined in style config file: {}");
    MACH3LOG_ERROR("RGB arrays dont all have the same size, please fix that");
  }

  // root only likes arrays so convert to those
  Double_t stops[NRGBs];
  Double_t red[NRGBs];
  Double_t green[NRGBs];
  Double_t blue[NRGBs];
  for (int i = 0; i < NRGBs; i++)
  {
    stops[i] = stopVec[i];
    red[i] = redsVec[i];
    green[i] = greensVec[i];
    blue[i] = bluesVec[i];
  }

  // now actually set the palette
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

void StyleManager::setTH1Style(TH1 *hist, std::string styleName) const {
  // get the definition of the provided style from the config file
  YAML::Node TH1Styles = _styleConfig["TH1Styles"];
  YAML::Node styleDef = TH1Styles[styleName];

  if (styleDef["MarkerColor"])
  {
    hist->SetMarkerColor(styleDef["MarkerColor"].as<int>());
  }
  if (styleDef["MarkerStyle"])
  {
    hist->SetMarkerStyle(styleDef["MarkerStyle"].as<int>());
  }
  if (styleDef["FillColor"])
  {
    hist->SetFillColor(styleDef["FillColor"].as<int>());
  }
  if (styleDef["FillStyle"])
  {
    hist->SetFillStyle(styleDef["FillStyle"].as<int>());
  }
  if (styleDef["LineColor"])
  {
    hist->SetLineColor(styleDef["LineColor"].as<int>());
  }
  if (styleDef["LineStyle"])
  {
    hist->SetLineStyle(styleDef["LineStyle"].as<int>());
  }
}

} // namespace MaCh3Plotting