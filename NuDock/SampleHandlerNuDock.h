#include "Samples/SampleHandlerFD.h"
#include "nudock.hpp"

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

  std::unordered_map<std::string, std::string> NuDockOscNameMap = {
    {"Theta12", "sin2th_12"},
    {"Theta13", "sin2th_13"},
    {"Theta23", "sin2th_23"},
    {"DeltaCP", "delta_cp"},
    {"Deltam2_21", "delm2_12"},
    {"Deltam2_32", "delm2_23"},
  };

  std::unordered_map<std::string, std::string> NuDockOscNameMap_r = {
    {"sin2th_12", "Theta12"},
    {"sin2th_13", "Theta13"},
    {"sin2th_23", "Theta23"},
    {"delta_cp", "DeltaCP"},
    {"delm2_12", "Deltam2_21"},
    {"delm2_23", "Deltam2_32"},
  };
protected:
  void CleanMemoryBeforeFit() override {};
  void Init() override;
private:
  // Pointer to NuDock client
  std::unique_ptr<NuDock> nudock_ptr;
  std::unique_ptr<manager> SampleManager;
  bool verbose;
};
