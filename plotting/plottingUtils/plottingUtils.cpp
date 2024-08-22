#include "plottingUtils.h"

namespace MaCh3Plotting {

/// It will go through your provided graph and make a histogram binning by taking the midpoints of
/// all the graph points x values, then fill the histogram with the graphs y values. For the first
/// and last points it will extend the binning out of the graph bounds using the width between the
/// outermost and second outermost points. This can be useful if e.g. you want to draw cumulative
/// stacks of LLH scans.
TH1D TGraphToTH1D(TGraph graph, std::string newName, std::string newTitle) {
  std::string name;
  std::string title;

  if (newName == "")
    name = graph.GetName();
  else
    name = newName;

  if (newTitle == "")
    title = graph.GetTitle();
  else
    title = newTitle;

  int nPoints = graph.GetN();
  Double_t pointsX[nPoints];
  Double_t pointsY[nPoints];

  // Get the points out
  Double_t x, y;

  for (int pointId = 0; pointId < nPoints; pointId++)
  {
    graph.GetPoint(pointId, x, y);
    pointsX[pointId] = (Double_t)x;
    pointsY[pointId] = (Double_t)y;
  }

  // get the bin edges
  Double_t binEdges[nPoints + 1];
  binEdges[0] = pointsX[0] - (pointsX[1] - pointsX[0]) / 2.0;
  binEdges[nPoints] = pointsX[nPoints - 1] + (pointsX[nPoints - 1] - pointsX[nPoints - 2]) / 2.0;

  for (int pointId = 1; pointId < nPoints; pointId++)
  {
    // take the midpoint of the two surrounding points
    binEdges[pointId] = (pointsX[pointId] + pointsX[pointId - 1]) / 2.0;
  }

  TH1D retHist = TH1D(name.c_str(), title.c_str(), nPoints, binEdges);

  for (int binId = 0; binId < nPoints; binId++)
  {
    retHist.SetBinContent(binId + 1, pointsY[binId]);
  }

  return retHist;
}

} // namespace MaCh3Plotting