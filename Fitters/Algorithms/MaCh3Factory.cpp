// MaCh3 includes
#include "Fitters/MaCh3Factory.h"

// ********************************************
std::unique_ptr<FitterBase> MaCh3FitterFactory(Manager *fitMan) {
// ********************************************
  std::unique_ptr<FitterBase> MaCh3Fitter = nullptr;

  auto Algorithm = GetFromManager<std::string>(fitMan->raw()["General"]["FittingAlgorithm"], "MCMC", __FILE__ , __LINE__);

  if(Algorithm == "MCMC" || Algorithm == "MR2T2") 
  {
    MaCh3Fitter = std::make_unique<MR2T2>(fitMan);
  }
  else if (Algorithm == "DelayedMCMC"  || Algorithm == "DelayedMR2T2")
  {
    MaCh3Fitter = std::make_unique<DelayedMR2T2>(fitMan);
  }
  else if (Algorithm == "PSO")
  {
    MaCh3Fitter = std::make_unique<PSO>(fitMan);
  }
  else if (Algorithm == "Minuit2")
  {
    #ifdef MaCh3_MINUIT2
    MaCh3Fitter = std::make_unique<MinuitFit>(fitMan);
    #else
    MACH3LOG_ERROR("Trying to use Minuit2 however MaCh3 was compiled without Minuit2 support");
    throw MaCh3Exception(__FILE__ , __LINE__ );
    #endif
  }
  else
  {
    MACH3LOG_ERROR("You want to use algorithm {}, I don't recognize it, sry", Algorithm);
    throw MaCh3Exception(__FILE__ , __LINE__ );
  }
  return MaCh3Fitter;
}

// ********************************************
std::unique_ptr<Manager> MaCh3ManagerFactory(int argc, char **argv) {
// ********************************************
  if (argc < 2) {
    MACH3LOG_ERROR("Wrong usage of MaCh3 executable!");
    MACH3LOG_ERROR("Syntax is $: {} config.yaml", argv[0]);
    MACH3LOG_ERROR("Where config.yaml is a valid config file, compatible with the Manager class (Manager/Manager.cpp/h)");
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  std::string arg2;
  if (argc > 2) {
    arg2 = argv[2];
  }
  // Check if we are using --override mode, or for some CLI whether we simply passed .yaml
  if ((argc == 4 && arg2 == "--override") ||
      (argc == 3 && arg2.size() >= 5 && arg2.compare(arg2.size() - 5, 5, ".yaml") == 0)) {
    std::string overrideFile;

    if (arg2 == "--override") {
      overrideFile = argv[3];
    } else {
      overrideFile = arg2;
    }
    MACH3LOG_INFO("Merging configuration files: base config '{}', override config '{}'. "
                  "Options in '{}' will take precedence over '{}'.",
                  argv[1], overrideFile, overrideFile, argv[1]);

    // Load the two YAML files
    YAML::Node config1 = M3OpenConfig(argv[1]);
    YAML::Node config2 = M3OpenConfig(overrideFile);

    // Merge them
    YAML::Node merged = MergeNodes(config1, config2);
    auto FitManager = std::make_unique<Manager>(merged);

    return FitManager;
  }

  // Initialise manger responsible for config handling
  auto FitManager = std::make_unique<Manager>(argv[1]);
  
  //KS: Lambda to make sure we are not overwriting setting which should be committed
  auto SanityOverwrite = [](const std::string& Name) {
    if (Name.find("Systematics") != std::string::npos ||
        Name.find("Samples") != std::string::npos)
    {
      MACH3LOG_CRITICAL("You are overwriting settings ({}) that are highly likely intended to be committed.", Name);
      /// @todo DL: Should probably replace this with something that doesn't require modifying core code
      MACH3LOG_CRITICAL("If you're sure you want to do this, e.g. for testing or step size tuning, you can remove the throw that lives here:");
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
    } else if (colonCount == 4) {
      const size_t firstColon = arg.find(':');
      const size_t secondColon = arg.find(':', firstColon + 1);
      const size_t thirdColon = arg.find(':', secondColon + 1);
      const size_t fourthColon = arg.find(':', thirdColon + 1);

      const std::string section = arg.substr(0, firstColon);
      const std::string key = arg.substr(firstColon + 1, secondColon - firstColon - 1);
      const std::string key2 = arg.substr(secondColon + 1, thirdColon - secondColon - 1);
      const std::string key3 = arg.substr(thirdColon + 1, fourthColon - thirdColon - 1);
      const std::string value = arg.substr(fourthColon + 1);

      MACH3LOG_INFO(
        "Overriding setting: Section={}, Key={}, Key={}, Key={}, Value={}",
        section, key, key2, key3, value);

      SanityOverwrite(section);
      SanityOverwrite(key);
      SanityOverwrite(key2);
      SanityOverwrite(key3);

      FitManager->OverrideSettings(section, key, key2, key3, value);
    } else {
      MACH3LOG_ERROR("Invalid override argument format: {}", arg);
      MACH3LOG_ERROR("Expected format:Section:Key:Key:Value, Section:Key:Value or Section:Value");
      throw MaCh3Exception(__FILE__, __LINE__);
    }
  }
  return FitManager;
}
