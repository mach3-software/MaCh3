#include "GenericBinningTools.h"

#include "spdlog/fmt/bundled/ranges.h"

namespace {
auto get_bin_edges_from_TAxis(TAxis const &ax) {
  std::vector<double> bin_edges = {ax.GetBinLowEdge(1)};
  for (int bi = 0; bi < ax.GetNbins(); ++bi) {
    bin_edges.push_back(ax.GetBinUpEdge(bi + 1));
  }
  return bin_edges;
}
} // namespace

std::unique_ptr<TH1> GetGenericBinningTH1(samplePDFFDBase &sample,
                                          std::string const &hname,
                                          std::string const &htitle,
                                          bool divide_by_ND_hypervolume) {

  auto const flathist = sample.get1DHist();
  auto const &binning = sample.generic_binning;

  if (!binning.GetNDimensions()) {
    MACH3LOG_ERROR("GetGenericBinningTH1 passed a samplePDFFDBase that isn't "
                   "using generic binning");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  std::unique_ptr<TH1> hout;

  int naxes = binning.Axes.size();
  if (naxes == 1) {
    // this is a bit silly, but ROOT doesn't have constructors that take TAxis
    // instances. So, really, who's the silly one?
    auto bin_edges = get_bin_edges_from_TAxis(binning.Axes[0]);

    hout = std::make_unique<TH1D>(hname.c_str(), htitle.c_str(),
                                  bin_edges.size() - 1, bin_edges.data());
  } else {
    hout = std::make_unique<TH1D>(
        hname.c_str(), htitle.c_str(), flathist->GetXaxis()->GetNbins(), -0.5,
        double(flathist->GetXaxis()->GetNbins()) - 0.5);
  }

  for (int gbi = 0; gbi < flathist->GetXaxis()->GetNbins(); ++gbi) {
    double bw = divide_by_ND_hypervolume
                    ? 1.0 / binning.GetGlobalBinHyperVolume(gbi)
                    : 1.0;
    hout->SetBinContent(gbi + 1, flathist->GetBinContent(gbi + 1) * bw);
    hout->SetBinError(gbi + 1, flathist->GetBinError(gbi + 1) * bw);
  }

  hout->GetXaxis()->SetTitle(binning.Axes[0].GetTitle());

  if (htitle.size()) {
    hout->SetTitle(htitle.c_str());
  }

  hout->SetDirectory(nullptr);
  return hout;
}

std::unique_ptr<TH2> GetGenericBinningTH2(samplePDFFDBase &sample,
                                          std::string const &hname,
                                          std::string const &htitle,
                                          bool divide_by_ND_hypervolume) {

  auto const flathist = sample.get1DHist();
  auto const &binning = sample.generic_binning;

  if (!binning.GetNDimensions()) {
    MACH3LOG_ERROR("GetGenericBinningTH2 passed a samplePDFFDBase that isn't "
                   "using generic binning");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  int naxes = binning.Axes.size();
  if (naxes != 2) {
    MACH3LOG_ERROR("Can only use GetGenericBinningTH2 on a 2D histogram, but "
                   "this one is {}D",
                   naxes);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // this is a bit silly, but ROOT doesn't have constructors that take TAxis
  // instances. So, really, who's the silly one?
  auto xbin_edges = get_bin_edges_from_TAxis(binning.Axes[0]);
  auto ybin_edges = get_bin_edges_from_TAxis(binning.Axes[1]);

  auto hout = std::make_unique<TH2D>(hname.c_str(), htitle.c_str(),
                                     xbin_edges.size() - 1, xbin_edges.data(),
                                     ybin_edges.size() - 1, ybin_edges.data());

  for (int gbi = 0; gbi < flathist->GetXaxis()->GetNbins(); ++gbi) {
    auto axis_binning = binning.DecomposeGlobalBinNumber(gbi);
    double bw = divide_by_ND_hypervolume
                    ? 1.0 / binning.GetGlobalBinHyperVolume(gbi)
                    : 1.0;
    hout->SetBinContent(axis_binning[0] + 1, axis_binning[1] + 1,
                        flathist->GetBinContent(gbi + 1) * bw);
    hout->SetBinError(axis_binning[0] + 1, axis_binning[1] + 1,
                      flathist->GetBinError(gbi + 1) * bw);
  }

  hout->GetXaxis()->SetTitle(binning.Axes[0].GetTitle());
  hout->GetYaxis()->SetTitle(binning.Axes[1].GetTitle());

  if (htitle.size()) {
    hout->SetTitle(htitle.c_str());
  }

  hout->SetDirectory(nullptr);
  return hout;
}

std::unique_ptr<TH3> GetGenericBinningTH3(samplePDFFDBase &sample,
                                          std::string const &hname,
                                          std::string const &htitle,
                                          bool divide_by_ND_hypervolume) {

  auto const flathist = sample.get1DHist();
  auto const &binning = sample.generic_binning;

  if (!binning.GetNDimensions()) {
    MACH3LOG_ERROR("GetGenericBinningTH3 passed a samplePDFFDBase that isn't "
                   "using generic binning");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  int naxes = binning.Axes.size();
  if (naxes != 3) {
    MACH3LOG_ERROR("Can only use GetGenericBinningTH3 on a 3D histogram, but "
                   "this one is {}D",
                   naxes);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // this is a bit silly, but ROOT doesn't have constructors that take TAxis
  // instances. So, really, who's the silly one?
  auto xbin_edges = get_bin_edges_from_TAxis(binning.Axes[0]);
  auto ybin_edges = get_bin_edges_from_TAxis(binning.Axes[1]);
  auto zbin_edges = get_bin_edges_from_TAxis(binning.Axes[2]);

  auto hout = std::make_unique<TH3D>(hname.c_str(), htitle.c_str(),
                                     xbin_edges.size() - 1, xbin_edges.data(),
                                     ybin_edges.size() - 1, ybin_edges.data(),
                                     zbin_edges.size() - 1, zbin_edges.data());

  for (int gbi = 0; gbi < flathist->GetXaxis()->GetNbins(); ++gbi) {
    auto axis_binning = binning.DecomposeGlobalBinNumber(gbi);
    double bw = divide_by_ND_hypervolume
                    ? 1.0 / binning.GetGlobalBinHyperVolume(gbi)
                    : 1.0;
    hout->SetBinContent(axis_binning[0] + 1, axis_binning[1] + 1,
                        axis_binning[2] + 1,
                        flathist->GetBinContent(gbi + 1) * bw);
    hout->SetBinError(axis_binning[0] + 1, axis_binning[1] + 1,
                      axis_binning[2] + 1, flathist->GetBinError(gbi + 1) * bw);
  }

  hout->GetXaxis()->SetTitle(binning.Axes[0].GetTitle());
  hout->GetYaxis()->SetTitle(binning.Axes[1].GetTitle());
  hout->GetZaxis()->SetTitle(binning.Axes[2].GetTitle());

  if (htitle.size()) {
    hout->SetTitle(htitle.c_str());
  }

  hout->SetDirectory(nullptr);
  return hout;
}

std::unique_ptr<TH1>
GetGenericBinningTH1Slice(samplePDFFDBase &sample,
                          std::vector<double> slice_definition,
                          std::string const &hname, std::string const &htitle,
                          bool divide_by_ND_hypervolume) {

  auto const flathist = sample.get1DHist();
  auto const &binning = sample.generic_binning;

  if (binning.GetNDimensions() < 2) {
    MACH3LOG_ERROR(
        "GetGenericBinningTH1Slice passed a samplePDFFDBase that isn't "
        "using generic binning of dimension >2, dimension = {}",
        binning.GetNDimensions());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (slice_definition.size() != binning.GetNDimensions()) {
    MACH3LOG_ERROR(
        "GetGenericBinningTH1Slice operating on a generic binning of {} "
        "dimensions, but passed a slice_definition with {} entries.",
        slice_definition.size(), binning.GetNDimensions());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  int nslice = 0;
  size_t slice_ax = 0;
  for (size_t axi = 0; axi < binning.Axes.size(); ++axi) {
    if (slice_definition[axi] == kSliceAx) {
      nslice++;
      slice_ax = axi;
    }
  }

  if (nslice != 1) {
    MACH3LOG_ERROR(
        "GetGenericBinningTH1Slice operating on the slice_definition: {} "
        "and expected to find one entry == kSliceAx(={}), specifying the axis "
        "to preserve, but found {}.",
        slice_definition, kSliceAx, nslice);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  std::stringstream cut_title;
  for (size_t axi = 0; axi < binning.Axes.size(); ++axi) {
    if (axi == slice_ax) {
      continue;
    }

    auto const &ax = binning.Axes[axi];
    int ax_bin = ax.FindFixBin(slice_definition[axi]);
    cut_title << fmt::format("{} < {} < {}{}", ax.GetBinLowEdge(ax_bin),
                             ax.GetTitle(), ax.GetBinUpEdge(ax_bin),
                             ((axi + 1) == binning.Axes.size()) ? "" : ", ");
  }

  auto const &xax = binning.Axes[slice_ax];

  // this is a bit silly, but ROOT doesn't have constructors that take TAxis
  // instances. So, really, who's the silly one?
  auto bin_edges = get_bin_edges_from_TAxis(xax);

  auto hout = std::make_unique<TH1D>(hname.c_str(), htitle.c_str(),
                                     bin_edges.size() - 1, bin_edges.data());

  for (int sbi = 0; sbi < xax.GetNbins(); ++sbi) {

    slice_definition[slice_ax] = xax.GetBinCenter(sbi + 1);

    int gbi = binning.GetGlobalBinNumber(slice_definition);
    double bw = divide_by_ND_hypervolume
                    ? 1.0 / binning.GetGlobalBinHyperVolume(gbi)
                    : 1.0;
    hout->SetBinContent(sbi + 1, flathist->GetBinContent(gbi + 1) * bw);
    hout->SetBinError(sbi + 1, flathist->GetBinError(gbi + 1) * bw);
  }

  hout->GetXaxis()->SetTitle(xax.GetTitle());

  if (htitle.size()) {
    hout->SetTitle(htitle.c_str());
  }

  // add the slice definition to the title
  std::string hout_title = hout->GetTitle();
  hout_title += cut_title.str();
  hout->SetTitle(hout_title.c_str());

  hout->SetDirectory(nullptr);
  return hout;
}

std::unique_ptr<TH2>
GetGenericBinningTH2Slice(samplePDFFDBase &sample,
                          std::vector<double> slice_definition,
                          std::string const &hname, std::string const &htitle,
                          bool divide_by_ND_hypervolume) {

  auto const flathist = sample.get1DHist();
  auto const &binning = sample.generic_binning;

  if (binning.GetNDimensions() < 3) {
    MACH3LOG_ERROR(
        "GetGenericBinningTH2Slice passed a samplePDFFDBase that isn't "
        "using a generic binning of dimension >3 dimension = {}",
        binning.GetNDimensions());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (slice_definition.size() != binning.GetNDimensions()) {
    MACH3LOG_ERROR(
        "GetGenericBinningTH2Slice operating on a generic binning of {} "
        "dimensions, but passed a slice_definition with {} entries.",
        slice_definition.size(), binning.GetNDimensions());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  int nslice = 0;
  size_t xslice_ax = std::numeric_limits<size_t>::max();
  size_t yslice_ax = std::numeric_limits<size_t>::max();
  for (size_t axi = 0; axi < binning.Axes.size(); ++axi) {
    if (slice_definition[axi] == kSliceAx) {
      nslice++;
      if (xslice_ax == std::numeric_limits<size_t>::max()) {
        xslice_ax = axi;
      } else {
        yslice_ax = axi;
      }
    }
  }

  if (nslice != 2) {
    MACH3LOG_ERROR(
        "GetGenericBinningTH2Slice operating on the slice_definition: {} "
        "and expected to find two entries == kSliceAx(={}), specifying the "
        "axes to preserve, but found {}.",
        slice_definition, kSliceAx, nslice);
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  std::stringstream cut_title;
  for (size_t axi = 0; axi < binning.Axes.size(); ++axi) {
    if ((axi == xslice_ax) || (axi == yslice_ax)) {
      continue;
    }

    auto const &ax = binning.Axes[axi];
    int ax_bin = ax.FindFixBin(slice_definition[axi]);
    cut_title << fmt::format("{:.4g} < {} < {:.4g}{}", ax.GetBinLowEdge(ax_bin),
                             ax.GetTitle(), ax.GetBinUpEdge(ax_bin),
                             ((axi + 1) == binning.Axes.size()) ? "" : ", ");
  }

  auto const &xax = binning.Axes[xslice_ax];
  auto const &yax = binning.Axes[yslice_ax];

  // this is a bit silly, but ROOT doesn't have constructors that take TAxis
  // instances. So, really, who's the silly one?
  auto xbin_edges = get_bin_edges_from_TAxis(xax);
  auto ybin_edges = get_bin_edges_from_TAxis(yax);

  auto hout = std::make_unique<TH2D>(hname.c_str(), htitle.c_str(),
                                     xbin_edges.size() - 1, xbin_edges.data(),
                                     ybin_edges.size() - 1, ybin_edges.data());

  for (int xsbi = 0; xsbi < xax.GetNbins(); ++xsbi) {
    for (int ysbi = 0; ysbi < yax.GetNbins(); ++ysbi) {

      slice_definition[xslice_ax] = xax.GetBinCenter(xsbi + 1);
      slice_definition[yslice_ax] = yax.GetBinCenter(ysbi + 1);

      int gbi = binning.GetGlobalBinNumber(slice_definition);
      double bw = divide_by_ND_hypervolume
                      ? 1.0 / binning.GetGlobalBinHyperVolume(gbi)
                      : 1.0;

      hout->SetBinContent(xsbi + 1, ysbi + 1,
                          flathist->GetBinContent(gbi + 1) * bw);
      hout->SetBinError(xsbi + 1, ysbi + 1,
                        flathist->GetBinError(gbi + 1) * bw);
    }
  }

  hout->GetXaxis()->SetTitle(xax.GetTitle());
  hout->GetYaxis()->SetTitle(yax.GetTitle());

  if (htitle.size()) {
    hout->SetTitle(htitle.c_str());
  }

  // add the slice definition to the title
  std::string hout_title = hout->GetTitle();
  hout_title += cut_title.str();
  hout->SetTitle(hout_title.c_str());

  hout->SetDirectory(nullptr);
  return hout;
}

std::vector<std::unique_ptr<TH1>> GetGenericBinningTH1Slices(
    samplePDFFDBase &sample, size_t slice_ax, std::string const &hname_pat,
    std::string const &htitle, bool divide_by_ND_hypervolume) {
  auto const &binning = sample.generic_binning;

  if (binning.GetNDimensions() < 2) {
    MACH3LOG_ERROR(
        "GetGenericBinningTH1XSlices passed a samplePDFFDBase that isn't "
        "using a generic binning of dimension >1 dimension = {}",
        binning.GetNDimensions());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  std::vector<std::unique_ptr<TH1>> out;

  // loop through all global bins
  for (int gbi = 0; gbi < binning.nbins_per_slice.back(); ++gbi) {
    // only want one output per slice, so skip all global bins where ax_bin !=
    // 0
    auto axis_bins = binning.DecomposeGlobalBinNumber(gbi);
    if (axis_bins[slice_ax] != 0) {
      continue;
    }

    std::stringstream ss;
    ss << fmt::format("ax{}slice_", slice_ax);
    // build the slice definition for this bin, could do it directly with bin
    // numbers, but this uses the existing public interface and isn't much extra
    // work.
    std::vector<double> slice_definition(axis_bins.size());
    for (int axi = 0; axi < binning.Axes.size(); ++axi) {
      slice_definition[axi] =
          binning.Axes[axi].GetBinCenter(axis_bins[axi] + 1);
      ss << (axi == slice_ax) ? fmt::format("along{}", axi)
                              : fmt::format("bin{}", axis_bins[axi]);
    }
    slice_definition[slice_ax] = kSliceAx;

    out.push_back(GetGenericBinningTH1Slice(
        sample, slice_definition, fmt::format("{}_{}", hname_pat, ss.str()),
        htitle, divide_by_ND_hypervolume));
  }

  return out;
}
