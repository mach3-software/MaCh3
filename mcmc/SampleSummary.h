#pragma once

// MaCh3 includes
#include "samplePDF/samplePDFBase.h"
#include "mcmc/StatisticalUtils.h"
#include "mcmc/MCMCProcessor.h"

namespace M3 {
  /// @brief KS: Different Information Criterion tests mostly based Gelman paper
  enum kInfCrit {
    kBIC,      //!< Bayesian Information Criterion
    kDIC,      //!< Deviance Information Criterion
    kWAIC,     //!< Watanabe-Akaike information criterion
    kInfCrits  //!< This only enumerates
  };
}
// *******************
/// @brief Class to calculate pvalue produce posterior predictive and many fancy Bayesian stuff \cite gelman1996posterior
/// @details For more information, visit the [Wiki](https://github.com/mach3-software/MaCh3/wiki/10.-Posterior-Predictive,-p%E2%80%90value-etc.).
/// @author Clarence Wret
/// @author Kamil Skwarczynski
class SampleSummary {
// *******************
  public:
    /// @brief Constructor
    /// @param n_Samples total number of samples
    /// @param Filename name of output file
    /// @param sample pointer to sample PDF object
    /// @param nChainSteps number of steps in a chain, 0 indicate prior predictive was used
    SampleSummary(const int n_Samples, const std::string &Filename, samplePDFBase* const sample, const int nSteps);
    /// @brief Destructor
    ~SampleSummary();

    /// @brief KS: Add data histograms
    /// @param DataHist Histogram with data even rates for each sample
    void AddData(std::vector<TH2Poly*> &DataHist);
    /// @brief KS: Add prior histograms
    void AddNominal(std::vector<TH2Poly*> &NominalHist, std::vector<TH2Poly*> &W2Nom);
    /// @brief KS: Add histograms with throws
    void AddThrow(std::vector<TH2Poly*> &MCHist, std::vector<TH2Poly*> &W2Hist, const double LLHPenalty = 0.0, const double Weight = 1.0, const int DrawNumber = 0);
    /// @brief KS: Add histograms for each mode
    void AddThrowByMode(std::vector<std::vector<TH2Poly*>> &SampleVector_ByMode);

    /// @brief KS: Write results into root file
    void Write();

    /// @brief KS: Set likelihood type
    inline void SetLikelihood(const TestStatistic TestStat){ likelihood = TestStat;};
    /// @brief Set number of model params used for BIC
    inline void SetNModelParams(const int nPars){ nModelParams = nPars;};

  private:
    /// @brief Finalise the distributions from the thrown samples
    inline void MakePredictive();

    /// @brief KS: Prepare output tree and necessary variables
    inline void PrepareOutput();

    /// @brief Helper functions to calculate likelihoods using TH2Poly, will modify MC hist tittle to include LLH
    /// @param Data histogram with data distribution for a single sample
    /// @param MC histogram with MC distribution for a single sample
    /// @param W2 histogram with W2 distribution for a single sample
    inline void CalcLLH(TH2Poly * const & Data, TH2Poly * const & MC, TH2Poly * const & W2);
    /// @brief Helper functions to calculate likelihoods using TH1D, will modify MC hist tittle to include LLH
    /// @param Data histogram with data distribution for a single sample
    /// @param MC histogram with MC distribution for a single sample
    /// @param W2 histogram with W2 distribution for a single sample
    inline void CalcLLH(TH1D * const & Data, TH1D * const & MC, TH1D * const & W2);

    /// @brief Helper functions to calculate likelihoods using TH2Poly
    /// @param Data histogram with data distribution for a single sample
    /// @param MC histogram with MC distribution for a single sample
    /// @param W2 histogram with W2 distribution for a single sample
    inline double GetLLH(TH2Poly * const & Data, TH2Poly * const & MC, TH2Poly * const & W2);
    /// @brief Helper functions to calculate likelihoods using TH1D
    /// @param Data histogram with data distribution for a single sample
    /// @param MC histogram with MC distribution for a single sample
    /// @param W2 histogram with W2 distribution for a single sample
    inline double GetLLH(TH1D * const & Data, TH1D * const & MC, TH1D * const & W2);

    /// @brief KS: In Barlow Beeston we have Beta Parameters which scale generated MC
    inline void PlotBetaParameters();

    /// @brief KS: Study how correlated are sample or kinematic bins
    inline void StudyKinematicCorrelations();

    /// @brief Make the cut LLH histogram
    inline void MakeCutLLH();
    // Make the 1D cut distribution and give the 1D p-value
    inline void MakeCutLLH1D(TH1D *Histogram, double llh_ref = -999);
    /// @brief Make the 2D cut distribution and give the 2D p-value
    inline void MakeCutLLH2D(TH2D *Histogram);
    /// @brief Make the 1D Event Rate Hist
    inline void MakeCutEventRate(TH1D *Histogram, const double DataRate);
    /// @brief Make the fluctuated histograms (2D and 1D) for the chi2s
    /// Essentially taking the MCMC draws and calculating their LLH to the Posterior predictive distribution
    /// And additionally taking the data histogram and calculating the LLH to the predictive distribution
    /// Additionally we calculate the chi2 of the draws (fluctuated) of  the MC with the prior/posterior predictive and plot it vs the chi2 from the draws of MCMC and the data
    inline void MakeChi2Hists();

    /// @brief Check the length of samples agrees
    inline bool CheckSamples(const int Length);

    /// @brief Helper to project TH2D onto axis
    inline TH1D* ProjectHist(TH2D* Histogram, const bool ProjectX);
    /// @brief Helper to project TH2Poly onto axis
    inline TH1D* ProjectPoly(TH2Poly* Histogram, const bool ProjectX, const int selection, const bool MakeErrorHist = false);

    /// @brief Make Poisson fluctuation of TH1D hist
    inline void MakeFluctuatedHistogram(TH1D *FluctHist, TH1D* PolyHist);

    /// @brief Make Poisson fluctuation of TH2Poly hist
    inline void MakeFluctuatedHistogram(TH2Poly *FluctHist, TH2Poly* PolyHist);

    /// @brief Information Criterion
    inline void StudyInformationCriterion(M3::kInfCrit Criterion);

    /// @brief Study Bayesian Information Criterion (BIC)
    /// @cite Gelman2014
    inline void StudyBIC();

    /// @brief KS: Get the Deviance Information Criterion (DIC)
    /// @cite Spiegelhalter2002
    /// @cite BRugsDIC
    inline void StudyDIC();

    /// @brief KS: Get the Watanabe-Akaike information criterion (WAIC)
    /// @cite Gelman2014
    /// @cite Hartig2024WAIC
    inline void StudyWAIC();

    /// Random number generator
    std::unique_ptr<TRandom3> rnd;
    /// KS: Hacky flag to let us know if this is first toy
    bool first_pass;

    /// KS: We have two methods for Poissonian fluctuation
    bool StandardFluctuation;

    /// Vector of vectors which holds the loaded MC histograms
    std::vector<std::vector<TH2Poly*>> MCVector;
    /// Vector of vectors which holds the loaded W2 histograms
    std::vector<std::vector<TH2Poly*>> W2MCVector;
    /// Vector of vectors which holds the loaded MC histograms for each mode
    std::vector<std::vector<std::vector<TH2Poly*>>> MCVectorByMode;

    /// Vector to hold the penalty term
    std::vector<double> LLHPenaltyVector;
    /// Vector holding weight
    std::vector<double> WeightVector;

    /// Number of samples
    int nSamples;

    /// name for each sample
    std::vector<std::string> SampleNames;

    /// The posterior predictive for the whole selection: this gets built after adding in the toys. Now an array of Th1ds, 1 for each poly bin, for each sample
    std::vector<std::vector<std::unique_ptr<TH1D>>> PosteriorHist;
    /// The posterior predictive for the whole selection: this gets built after adding in the toys. Now an array of Th1ds, 1 for each poly bin, for each sample for W2
    std::vector<std::vector<std::unique_ptr<TH1D>>> w2Hist;

    /// Posterior predictive but for X projection but as a violin plot
    std::vector<TH2D*> ViolinHists_ProjectX;
    /// Posterior predictive but for Y projection but as a violin plot
    std::vector<TH2D*> ViolinHists_ProjectY;
    
    /// The data histogram for the selection
    std::vector<TH2Poly*> DataHist;
    /// The data histogram for the selection X projection
    std::vector<TH1D*> DataHist_ProjectX;
    /// The data histogram for the selection Y projection
    std::vector<TH1D*> DataHist_ProjectY;
    /// The nominal histogram for the selection
    std::vector<TH2Poly*> NominalHist;
    /// Pointer to the w2 histograms (for nominal values).
    std::vector<TH2Poly*> W2NomHist;
    /// Pointer to the w2 histograms (for mean values).
    std::vector<TH2Poly*> W2MeanHist;
    /// Pointer to the w2 histograms (for mode values).
    std::vector<TH2Poly*> W2ModeHist;

    /// The histogram containing the lnL for each throw
    std::unique_ptr<TH1D> lnLHist;
    /// The lnLhist for the draw vs MC fluctuated
    std::unique_ptr<TH1D> lnLHist_drawfluc;
    /// The lnLhist for the draw vs draw fluctuated
    std::unique_ptr<TH1D> lnLHist_drawflucdraw;
    /// The lnLhist for the draw vs data
    std::unique_ptr<TH1D> lnLHist_drawdata;
    /// The 2D lnLhist, showing (draw vs data) and (draw vs fluct), anything above y=x axis is the p-value
    std::unique_ptr<TH2D> lnLDrawHist;
    /// The 2D lnLHist, showing (draw vs data) and (draw vs draw fluct), anything above y=x axis is the p-value
    std::unique_ptr<TH2D> lnLFlucHist;

    /// The 2D lnLhist, showing (draw vs data) and (draw vs fluct), using rate, anything above y=x axis is the p-value
    std::unique_ptr<TH2D> lnLDrawHistRate;
    /// The 2D lnLHist but for ProjectionX histogram (pmu), showing (draw vs data) and (draw vs draw fluct), anything above y=x axis is the p-value
    std::unique_ptr<TH2D> lnLFlucHist_ProjectX;

    /// The histogram containing the lnL (draw vs data) for each throw for each sample
    std::vector<TH1D*> lnLHist_Sample_DrawData;
    /// The histogram containing the lnL (draw vs draw fluct) for each throw for each sample
    std::vector<TH1D*> lnLHist_Sample_DrawflucDraw;
    /// The histogram containing the lnL (draw vs pred fluct) for each throw for each sample
    std::vector<TH1D*> lnLHist_Sample_PredflucDraw;

    /// The LLH distribution in pmu cosmu for using the mean in each bin
    std::vector<TH2Poly*> lnLHist_Mean;
    /// The LLH distribution in pmu cosmu for using the mode in each bin
    std::vector<TH2Poly*> lnLHist_Mode;

    /// The LLH distribution in pmu using the mean in each bin
    std::vector<TH1D*> lnLHist_Mean_ProjectX;

    /// The posterior predictive distribution in pmu cosmu using the mean
    std::vector<TH2Poly*> MeanHist;
    /// The posterior predictive distribution in pmu cosmu using the mean after applying Barlow-Beeston Correction
    std::vector<TH2Poly*> MeanHistCorrected;
    /// The posterior predictive distribution in pmu cosmu using the mode
    std::vector<TH2Poly*> ModeHist;

    /// Holds the bin-by-bin LLH for the mean posterior predictive vs the data
    std::vector<TH1D*> lnLHist_Mean1D;
    /// Holds the bin-by-bin LLH for the mode posterior predictive vs the data
    std::vector<TH1D*> lnLHist_Mode1D;

    /// Holds the history of which entries have been drawn in the MCMC file
    std::unique_ptr<TH1D> RandomHist;

    /// Distribution of beta parameters in Barlow Beeston formalisms
    std::vector<std::vector<std::unique_ptr<TH1D>>> BetaHist;
    /// Are we making Beta Histograms
    bool DoBetaParam;

    /// Number of throws by user
    unsigned int nChainSteps;

    /// bool whether we have Prior or Posterior Predictive
    bool isPriorPredictive;

    /// bool whether to normalise each toy to have shape based p-value and pos pred distribution
    bool doShapeOnly;

    /// Number of throws
    unsigned int nThrows;

    /// Max Number of Bins per each sample
    std::vector<int> maxBins;

    /// Total LLH for the posterior predictive distribution
    double llh_total;

    /// Output filename
    std::string OutputName;
    /// Output filename
    TFile *Outputfile;
    /// Directory for each sample
    std::vector<TDirectory*> Dir;

    /// TTree which we save useful information to
    TTree *OutputTree;
    /// Data vs Draw
    std::vector<double> llh_data_draw;
    /// Fluctuated Draw vs Draw
    std::vector<double> llh_drawfluc_draw;
    /// Fluctuated Predictive vs Draw
    std::vector<double> llh_predfluc_draw;

    /// Data vs Draw using rate only
    std::vector<double> llh_rate_data_draw;
    /// Fluctuated Predictive vs Draw using rate only
    std::vector<double> llh_rate_predfluc_draw;

    /// Data vs Fluctuated Draw
    std::vector<double> llh_data_drawfluc;
    /// Data vs Fluctuated Predictive
    std::vector<double> llh_data_predfluc;
    /// Draw vs Predictive
    std::vector<double> llh_draw_pred;
    /// Fluctuated Draw vs Predictive
    std::vector<double> llh_drawfluc_pred;

    /// Fluctuated Predictive vs Predictive
    std::vector<double> llh_predfluc_pred;
    /// Fluctuated Draw vs Fluctuated Predictive
    std::vector<double> llh_drawfluc_predfluc;
    /// Fluctuated Data vs Draw
    std::vector<double> llh_datafluc_draw;

    /// Projection X (most likely muon momentum) of LLH
    std::vector<double> llh_data_draw_ProjectX;
    std::vector<double> llh_drawfluc_draw_ProjectX;

    /// LLH penalty for each throw
    double llh_penalty;

    /// Data vs Draw
    double total_llh_data_draw;
    /// Fluctuated Draw vs Draw
    double total_llh_drawfluc_draw;
    /// Fluctuated Predictive vs Draw
    double total_llh_predfluc_draw;

    ///  Rate Data vs Draw
    double total_llh_rate_data_draw;
    /// Fluctuated Predictive vs Draw using Rate
    double total_llh_rate_predfluc_draw;

    /// Data vs Fluctuated Predictive
    double total_llh_data_predfluc;
    /// Data vs Fluctuated Draw
    double total_llh_data_drawfluc;
    /// Draw vs Predictive
    double total_llh_draw_pred;
    /// Fluctuated Draw vs Predictive
    double total_llh_drawfluc_pred;
    /// Fluctuated Draw vs Fluctuated Predictive
    double total_llh_drawfluc_predfluc;
    /// Fluctuated Data vs Draw
    double total_llh_datafluc_draw;
    /// Fluctuated Predictive vs Predictive
    double total_llh_predfluc_pred;

    /// Data vs Draw for projection X (most likely muon momentum)
    double total_llh_data_draw_ProjectX;
    /// Fluctuated Draw vs Draw for projection X (most likely muon momentum)
    double total_llh_drawfluc_draw_ProjectX;

    /// By mode variables
    bool DoByModePlots;
    /// The posterior predictive distribution in pmu cosmu using the mean
    std::vector<std::vector<TH2Poly*>> MeanHist_ByMode;
    /// Histogram which corresponds to each bin in the sample's th2poly
    TH1D ****PosteriorHist_ByMode;
    
    /// Pointer to samplePDF object, mostly used to get sample names, binning etc.
    samplePDFBase* SamplePDF;

    /// MaCh3 Modes
    MaCh3Modes* Modes;

    /// Type of likelihood for example Poisson, Barlow-Beeston or Ice Cube
    TestStatistic likelihood;

    /// Number of parameters
    int nModelParams;

    /// Tells Debug level to save additional histograms
    int Debug;
};
