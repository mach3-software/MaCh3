// MaCh3 includes
#include "mcmc/MaCh3Factory.h"

// ********************************************
std::unique_ptr<FitterBase> MaCh3FitterFactory(manager *fitMan) {
// ********************************************
  std::unique_ptr<FitterBase> MaCh3Fitter = nullptr;

  auto Algorithm = GetFromManager<std::string>(fitMan->raw()["General"]["FittingAlgorithm"], "MCMC");

  if(Algorithm == "MCMC") {
    MaCh3Fitter = std::make_unique<mcmc>(fitMan);
  } else if (Algorithm == "PSO") {
    MaCh3Fitter = std::make_unique<PSO>(fitMan);
  } else if (Algorithm == "Minuit2") {
    #ifdef MaCh3_MINUIT2
    MaCh3Fitter = std::make_unique<MinuitFit>(fitMan);
    #else
    MACH3LOG_ERROR("Trying to use Minuit2 however MaCh3 was compiled without Minuit2 support");
    throw MaCh3Exception(__FILE__ , __LINE__ );
    #endif
  } else {
    MACH3LOG_ERROR("You want to use algorithm {}, I don't recognize it, sry", Algorithm);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return MaCh3Fitter;
}

// ********************************************
covarianceXsec* MaCh3CovarianceFactory(manager *FitManager, const std::string& PreFix) {
// ********************************************
  return MaCh3CovarianceFactory<covarianceXsec>(FitManager, PreFix);
}

// ********************************************
std::unique_ptr<manager> MaCh3ManagerFactory(int argc, char **argv) {
// ********************************************
  if (argc < 2) {
    MACH3LOG_ERROR("Wrong usage of MaCh3 executable!");
    MACH3LOG_ERROR("Syntax is $: {} config.yaml", argv[0]);
    MACH3LOG_ERROR("Where config.yaml is a valid config file, compatible with the manager class (manager/manager.cpp/h)");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // Initialise manger responsible for config handling
  auto FitManager = std::make_unique<manager>(argv[1]);
  
  //KS: Lambda to make sure we are not overwriting setting which should be committed
  auto SanityOverwrite = [](const std::string& Name) {
    if (Name.find("Systematics") != std::string::npos ||
        Name.find("Samples") != std::string::npos)
    {
      MACH3LOG_CRITICAL("You are overwriting settings ({}) that are highly likely intended to be committed.", Name);
      throw MaCh3Exception(__FILE__ , __LINE__ );
    }
  };
  
  for (int i = 2; i < argc; ++i)
  {
    const std::string arg = argv[i];
    const size_t colonCount = std::count(arg.begin(), arg.end(), ':');

    /// @todo KS: May need some recursive magic to reduce amount of hardcoding
    if (colonCount == 1) {
      const size_t firstColon = arg.find(':');
      const std::string section = arg.substr(0, firstColon);
      const std::string value = arg.substr(firstColon + 1);

      MACH3LOG_INFO("Overriding setting: Section={}, Value={}", section, value);
      SanityOverwrite(section);
      FitManager->OverrideSettings(section, value);
    } else if (colonCount == 2) {
      const size_t firstColon = arg.find(':');
      const size_t secondColon = arg.find(':', firstColon + 1);

      const std::string section = arg.substr(0, firstColon);
      const std::string key = arg.substr(firstColon + 1, secondColon - firstColon - 1);
      const std::string value = arg.substr(secondColon + 1);

      MACH3LOG_INFO("Overriding setting: Section={}, Key={}, Value={}", section, key, value);
      SanityOverwrite(section);
      SanityOverwrite(key);
      FitManager->OverrideSettings(section, key, value);
    } else if (colonCount == 3) {
      const size_t firstColon = arg.find(':');
      const size_t secondColon = arg.find(':', firstColon + 1);
      const size_t thridColon = arg.find(':', secondColon + 1);

      const std::string section = arg.substr(0, firstColon);
      const std::string key = arg.substr(firstColon + 1, secondColon - firstColon - 1);
      const std::string key2 = arg.substr(secondColon + 1, thridColon - secondColon - 1);
      const std::string value = arg.substr(thridColon + 1);

      MACH3LOG_INFO("Overriding setting: Section={}, Key={}, Key={}, Value={}", section, key, key2, value);
      SanityOverwrite(section);
      SanityOverwrite(key);
      SanityOverwrite(key2);
      FitManager->OverrideSettings(section, key, key2, value);
    } else {
      MACH3LOG_ERROR("Invalid override argument format: {}", arg);
      MACH3LOG_ERROR("Expected format:Section:Key:Key:Valu, Section:Key:Value or Section:Value");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
  return FitManager;
}
