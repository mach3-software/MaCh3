#include "styleManager.h"

namespace MaCh3Plotting {
StyleManager::StyleManager(std::string styleConfigName) {
  _styleConfig = M3OpenConfig(styleConfigName);
}

std::string StyleManager::prettifyName(const std::string &origName, const std::string &nameType) const {
  YAML::Node prettyNames = _styleConfig["PrettyNames"][nameType];
  auto prettyName = GetFromManager<std::string>(prettyNames[origName], origName, __FILE__, __LINE__);

  return prettyName;
}

void StyleManager::setPalette(const int rootPlotStyle) const {
  // set the colour palette to one of the root palettes
  gStyle->SetPalette(rootPlotStyle);
}

void StyleManager::setPalette(const std::string& configStyleName) const {
  // set the colour palette to one of the palettes defined in PlottingConfig.yaml

  // get the definition of the provided style from the config file
  YAML::Node palettes = _styleConfig["ColorPallettes"];
  YAML::Node styleDef;
  if (palettes) {
    styleDef = palettes[configStyleName];
  }

  auto paletteDef = GetFromManager<std::vector<std::vector<double>>>(
    styleDef,
    std::vector<std::vector<double>>{
      {4.0},                         // NCont
      {0.0   , 0.33  , 0.66 , 1.0},  // stops
      {0.0   , 0.0   , 0.0  , 0.0},  // Reds
      {0.0   , 0.0   , 1.0  , 0.0},  // Greens
      {0.0   , 1.0   , 0.0  , 0.0}   // Blues
    }, __FILE__, __LINE__);
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
    MACH3LOG_ERROR("invalid colour palette defined in style config file: {}");
    MACH3LOG_ERROR("RGB arrays don't all have the same size, please fix that");
  }

  // now actually set the palette
  TColor::CreateGradientColorTable(int(NRGBs), stopVec.data(), redsVec.data(), greensVec.data(), bluesVec.data(), NCont);
  gStyle->SetNumberContours(NCont);
}

void StyleManager::setTH1Style(TH1 *hist, const std::string& styleName) const {
  // get the definition of the provided style from the config file
  YAML::Node TH1Styles = _styleConfig["TH1Styles"];
  YAML::Node styleDef ;
  if (TH1Styles) {
    styleDef = TH1Styles[styleName];
  }

  hist->SetMarkerColor(GetFromManager<Color_t>(styleDef["MarkerColor"], kRed, __FILE__, __LINE__));
  hist->SetMarkerStyle(GetFromManager<Color_t>(styleDef["MarkerStyle"], 7, __FILE__, __LINE__));
  hist->SetFillColor(GetFromManager<Color_t>(styleDef["FillColor"], kRed, __FILE__, __LINE__));
  hist->SetFillStyle(GetFromManager<Color_t>(styleDef["FillStyle"], 3003, __FILE__, __LINE__));
  hist->SetLineColor(GetFromManager<Color_t>(styleDef["LineColor"], kRed, __FILE__, __LINE__));
  hist->SetLineStyle(GetFromManager<Color_t>(styleDef["LineStyle"], 1, __FILE__, __LINE__));
}

} // namespace MaCh3Plotting
