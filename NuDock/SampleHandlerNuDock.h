#include "Samples/SampleHandlerFD.h"
#include "NuDockFactory.h"
_MaCh3_Safe_Include_Start_ //{
#include "nudock.hpp"
_MaCh3_Safe_Include_End_ //}

class SampleHandlerNuDock : public SampleHandlerFD {
public:
  SampleHandlerNuDock(std::string configFile, ParameterHandlerGeneric* xsec_cov, const std::shared_ptr<OscillationHandler>& Oscillator = nullptr);
  virtual ~SampleHandlerNuDock();
  void Reweight() override;
  double GetLikelihood() const override;

  // HH: Some virtual functions
  void AddAdditionalWeightPointers() override {};
  void SetupSplines() override {};
  int SetupExperimentMC() override {return -1;};
  void SetupFDMC() override {};
  void RegisterFunctionalParameters() override {};
  double ReturnKinematicParameter(std::string KinematicParameter, int iEvent) override {(void) KinematicParameter; (void) iEvent; return -1.;};
  double ReturnKinematicParameter(int KinematicVariable, int iEvent) override {(void) KinematicVariable; (void) iEvent; return -1.;};
  const double* GetPointerToKinematicParameter(std::string KinematicParamter, int iEvent) override {(void) KinematicParamter; (void) iEvent; return nullptr;};
  const double* GetPointerToKinematicParameter(double KinematicVariable, int iEvent) override {(void) KinematicVariable; (void) iEvent; return nullptr;};

protected:
  void CleanMemoryBeforeFit() override {};
  void Init() override;
private:
  // Pointer to NuDock client
  std::unique_ptr<NuDock> nudock_ptr;
  std::unique_ptr<manager> SampleManager;
  bool verbose;
  std::vector<int> nudockParamInds;
};
