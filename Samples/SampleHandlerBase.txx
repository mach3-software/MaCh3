/// @brief HH - a helper function for RegisterFunctionalParameter
template <typename EventType, typename SFType>
void SampleHandlerBase::RegisterIndividualFunctionalParameter(
    std::vector<EventType> &ExptEvents,
    std::vector<std::string> const &par_names, SFType shift_func) {

  static_assert(
      std::is_same_v<
          std::function<void(std::vector<double> const &, EventType &)>,
          decltype(std::function(shift_func))>,
      "Function call signature for single parameter Functional shift must be "
      "void(std::vector<double> const &, EventType &).");

  auto sample_func_pars =
      ParHandler->GetFunctionalParametersFromSampleName(SampleHandlerName);

  std::stringstream ss_pars, ss_miss;
  std::vector<FunctionalParameter const *> matched_pars;
  for (auto const &par_name : par_names) {
    ss_pars << par_name << " ";
    bool found = false;
    for (auto const &fp : sample_func_pars) {
      if (fp.name == par_name) {
        matched_pars.push_back(&fp);
        found = true;
      }
    }
    if (!found) {
      ss_miss << par_name;
    }
  }

  // allows experiments to effectively disable functional parameters by not
  // supplying the YAML defining the parameters
  if (!matched_pars.size()) {
    MACH3LOG_INFO("Functional shift consuming parameters: [ {}], found no "
                  "matching parameters for "
                  "sample handler: {}",
                  ss_pars.str(), SampleHandlerName);
    return;
  } else if (matched_pars.size() !=
             par_names.size()) { // not well defined how to procede with
                                 // partially defined parameter sets, so don't
    MACH3LOG_ERROR(
        "Functional shift consuming parameters: [ {}], only partially "
        "applys to sample handler: {}, missed parameters: [ {}]",
        ss_pars.str(), SampleHandlerName, ss_miss.str());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (!functional.event_shifts.size()) {
    functional.event_shifts.resize(ExptEvents.size());
  } else if (functional.event_shifts.size() != ExptEvents.size()) {
    MACH3LOG_ERROR("When registering functional shift consuming parameters: "
                   "[ {}], SampleHandler: {} has an allocated event map of "
                   "size: {}, but passed a vector of experiment events size: "
                   "{}. SampleHandler must have a unique set of event "
                   "indices so this indicates something has gone wrong.",
                   ss_pars.str(), SampleHandlerName,
                   functional.event_shifts.size(), ExptEvents.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  if (ExptEvents.size() != MCEvents.size()) {
    MACH3LOG_ERROR("When registering functional shift consuming parameters: "
                   "[ {}], SampleHandler: {} knows about {} MCEvents, but "
                   "passed a vector of experiment events size: "
                   "{}. SampleHandler must have a unique set of event "
                   "indices so this indicates something has gone wrong.",
                   ss_pars.str(), SampleHandlerName, MCEvents.size(),
                   ExptEvents.size());
    throw MaCh3Exception(__FILE__, __LINE__);
  }

  // If we've made it here, then this functional shift is correctly
  //   configured for this SampleHandler

  int iShift = int(functional.shifts.size());
  M3::detail::Functional::Shift shift;

  std::vector<int> par_indices;
  for (auto const &fp : matched_pars) {
    par_indices.push_back(fp->index);
  }

  shift.par_vals.resize(par_indices.size());
  shift.update_vals = [=](std::vector<double> &par_vals) {
    for (size_t i = 0; i < par_indices.size(); ++i) {
      auto const &fpi = par_indices[i];
      par_vals[i] = ParHandler->GetParPropVec()[fpi];
    }
  };

  shift.apply = [shift_func, &ExptEvents](std::vector<double> const &par_vals,
                                          int iEvent) {
    shift_func(par_vals, ExptEvents[iEvent]);
  };

  functional.shifts.push_back(std::move(shift));

  // For each event, make a vector of pointers to the functional parameters
  int NEvents = GetNEvents();
  int neventsaffected = 0;
  for (int iEvent = 0; iEvent < NEvents; ++iEvent) {
    int nmatch = 0;
    for (auto const &par : matched_pars) {
      const int Mode = static_cast<int>(
          std::round(ReturnKinematicParameter("Mode", iEvent)));
      if (!MatchCondition(par->modes, Mode)) {
        MACH3LOG_TRACE("Event {}, missed Mode check ({}) for dial {}", iEvent,
                       Mode, par->name);
        break;
      }
      if (!PassesSelection((*par), iEvent)) {
        MACH3LOG_TRACE("Event {}, missed Kinematic var check for dial {}",
                       iEvent, par->name);
        break;
      }
      nmatch++;
    }
    if (!nmatch) {
      continue;
    }
    if (nmatch != matched_pars.size()) {
      MACH3LOG_ERROR(
          "When determining wether functional shift consuming parameters: "
          "[ {}], in SampleHandler: {} should apply to event: {}, only {}/{} "
          "parameters matched. Partially applied shifts are ill-defined.",
          ss_pars.str(), SampleHandlerName, iEvent, nmatch,
          matched_pars.size());
      throw MaCh3Exception(__FILE__, __LINE__);
    }
    functional.event_shifts[iEvent].push_back(iShift);
    neventsaffected++;
  }
  MACH3LOG_INFO(
      "Registered functional shift consuming parameters: "
      "[ {}] for SampleHandler: {}, with {} par vals affecting {} events.",
      ss_pars.str(), SampleHandlerName, shift.par_vals.size(), neventsaffected);
}

template <typename EventType, typename SFType>
void SampleHandlerBase::RegisterIndividualFunctionalParameter(
    std::vector<EventType> &ExptEvents, std::string const &par_name,
    SFType shift_func) {

  static_assert(
      std::is_same_v<std::function<void(double const &, EventType &)>,
                     decltype(std::function(shift_func))>,
      "Function call signature for single parameter Functional shift must be "
      "void(double const &, EventType&).");

  RegisterIndividualFunctionalParameter(
      ExptEvents,
      std::vector<std::string>{
          par_name,
      },
      [=](std::vector<double> const &par_vals, EventType &ev) {
        shift_func(par_vals[0], ev);
      });
}

// ***************************************************************************
template <typename ParT>
bool SampleHandlerBase::PassesSelection(const ParT& Par, std::size_t iEvent) {
// ***************************************************************************
  bool IsSelected = true;
  if (Par.hasKinBounds) {
    const auto& kinVars = Par.KinematicVarStr;
    const auto& selection = Par.Selection;

    for (std::size_t iKinPar = 0; iKinPar < kinVars.size(); ++iKinPar) {
      const double kinVal = ReturnKinematicParameter(kinVars[iKinPar], static_cast<int>(iEvent));

      bool passedAnyBound = false;
      const auto& boundsList = selection[iKinPar];

      for (const auto& bounds : boundsList) {
        if (kinVal > bounds[0] && kinVal <= bounds[1]) {
          passedAnyBound = true;
          break;
        }
      }

      if (!passedAnyBound) {
        MACH3LOG_TRACE("Event {}, missed kinematic check ({}) for dial {}",
                       iEvent, kinVars[iKinPar], Par.name);
        IsSelected = false;
        break;
      }
    }
  }
  return IsSelected;
}

