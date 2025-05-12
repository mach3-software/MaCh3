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

  const Int_t NCont = Int_t(paletteDef[0][0]);

  std::vector<double> stopVec = paletteDef[1];
  std::vector<double> redsVec = paletteDef[2];
  std::vector<double> greensVec = paletteDef[3];
  std::vector<double> bluesVec = paletteDef[4];

  // get the number of colour stops and check all vectors are same size
  const size_t NRGBs = stopVec.size();
  if (redsVec.size() != NRGBs || greensVec.size() != NRGBs ||
      bluesVec.size() != NRGBs)
  {
    MACH3LOG_ERROR("invalid colour palettet defined in style config file: {}");
    MACH3LOG_ERROR("RGB arrays dont all have the same size, please fix that");
  }

  // now actually set the palette
  TColor::CreateGradientColorTable(int(NRGBs), stopVec.data(), redsVec.data(), greensVec.data(), bluesVec.data(), NCont);
  gStyle->SetNumberContours(NCont);
}

void StyleManager::setTH1Style(TH1 *hist, std::string styleName) const {
  // get the definition of the provided style from the config file
  YAML::Node TH1Styles = _styleConfig["TH1Styles"];
  YAML::Node styleDef = TH1Styles[styleName];

  if (styleDef["MarkerColor"])
  {
    hist->SetMarkerColor(styleDef["MarkerColor"].as<Color_t>());
  }
  if (styleDef["MarkerStyle"])
  {
    hist->SetMarkerStyle(styleDef["MarkerStyle"].as<Color_t>());
  }
  if (styleDef["FillColor"])
  {
    hist->SetFillColor(styleDef["FillColor"].as<Color_t>());
  }
  if (styleDef["FillStyle"])
  {
    hist->SetFillStyle(styleDef["FillStyle"].as<Color_t>());
  }
  if (styleDef["LineColor"])
  {
    hist->SetLineColor(styleDef["LineColor"].as<Color_t>());
  }
  if (styleDef["LineStyle"])
  {
    hist->SetLineStyle(styleDef["LineStyle"].as<Color_t>());
  }
}

} // namespace MaCh3Plotting