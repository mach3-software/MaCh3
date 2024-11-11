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
  if(nPoints < 2){
    MACH3LOG_ERROR("Too few points in the graph.");
    throw MaCh3Exception(__FILE__,__LINE__);
  }
  std::vector<double> pointsX(nPoints);
  std::vector<double> pointsY(nPoints);

  // Get the points out
  Double_t x, y;

  for (int pointId = 0; pointId < nPoints; pointId++)
  {
    graph.GetPoint(pointId, x, y);
    pointsX[pointId] = x;
    pointsY[pointId] = y;
  }

  // get the bin edges
  std::vector<double> binEdges(nPoints + 1);
  binEdges[0] = pointsX[0] - (pointsX[1] - pointsX[0]) / 2.0;
  binEdges[nPoints] = pointsX[nPoints - 1] + (pointsX[nPoints - 1] - pointsX[nPoints - 2]) / 2.0;

  for (int pointId = 1; pointId < nPoints; pointId++)
  {
    // take the midpoint of the two surrounding points
    binEdges[pointId] = (pointsX[pointId] + pointsX[pointId - 1]) / 2.0;
  }

  TH1D retHist = TH1D(name.c_str(), title.c_str(), nPoints, binEdges.data());

  for (int binId = 0; binId < nPoints; binId++)
  {
    retHist.SetBinContent(binId + 1, pointsY[binId]);
  }

  return retHist;
}


std::vector<std::vector<double>> TGraphToVector(TGraph graph) {

  int nPoints = graph.GetN();
  std::vector<std::vector<double>> ret(2);
  std::vector<double> pointsX(nPoints);
  std::vector<double> pointsY(nPoints);

  // Get the points out
  Double_t x, y;

  for (int pointId = 0; pointId < nPoints; pointId++)
  {
    graph.GetPoint(pointId, x, y);
    pointsX[pointId] = x;
    pointsY[pointId] = y;
  }

  ret[0] = pointsX;
  ret[1] = pointsY;

  return ret;
}


std::vector<std::vector<double>> TGraphToVector(TGraph2D graph) {

  int nPoints = graph.GetN();
  std::vector<std::vector<double>> ret(3);
  std::vector<double> pointsX(nPoints);
  std::vector<double> pointsY(nPoints);
  std::vector<double> pointsZ(nPoints);

  // Get the points out
  Double_t x, y, z;

  for (int pointId = 0; pointId < nPoints; pointId++)
  {
    graph.GetPoint(pointId, x, y, z);
    pointsX[pointId] = x;
    pointsY[pointId] = y;
    pointsZ[pointId] = z;
  }

  ret[0] = pointsX;
  ret[1] = pointsY;
  ret[2] = pointsZ;

  return ret;
}

} // namespace MaCh3Plotting 
