#include <iostream>
#include <cmath>
#include <boost/math/special_functions/lambert_w.hpp>
#include <sstream>
#include <fstream>

#include "SimTracker/SiPhase2Digitizer/plugins/PixelDigitizerAlgorithm.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationOfflineSimService.h"

// Geometry
#include "CondFormats/SiPixelObjects/interface/GlobalPixel.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelLorentzAngle.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "CondFormats/SiPixelObjects/interface/SiPixelQuality.h"
#include "CondFormats/SiPixelObjects/interface/PixelROC.h"
#include "CondFormats/SiPixelObjects/interface/LocalPixel.h"
#include "CondFormats/SiPixelObjects/interface/CablingPathToDetUnit.h"

using namespace edm;
using namespace sipixelobjects;

void PixelDigitizerAlgorithm::init(const edm::EventSetup& es) {
  if (use_ineff_from_db_)  // load gain calibration service fromdb...
    theSiPixelGainCalibrationService_->setESObjects(es);

  if (use_deadmodule_DB_)
    siPixelBadModule_ = &es.getData(siPixelBadModuleToken_);

  if (use_LorentzAngle_DB_)  // Get Lorentz angle from DB record
    siPixelLorentzAngle_ = &es.getData(siPixelLorentzAngleToken_);

  // gets the map and geometry from the DB (to kill ROCs)
  fedCablingMap_ = &es.getData(fedCablingMapToken_);
  geom_ = &es.getData(geomToken_);
  if (useChargeReweighting_) {
    theSiPixelChargeReweightingAlgorithm_->init(es);
  }
}

PixelDigitizerAlgorithm::PixelDigitizerAlgorithm(const edm::ParameterSet& conf, edm::ConsumesCollector iC)
    : Phase2TrackerDigitizerAlgorithm(conf.getParameter<ParameterSet>("AlgorithmCommon"),
                                      conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm"),
                                      iC),
      odd_row_interchannelCoupling_next_row_(conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm")
                                                 .getParameter<double>("Odd_row_interchannelCoupling_next_row")),
      even_row_interchannelCoupling_next_row_(conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm")
                                                  .getParameter<double>("Even_row_interchannelCoupling_next_row")),
      odd_column_interchannelCoupling_next_column_(
          conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm")
              .getParameter<double>("Odd_column_interchannelCoupling_next_column")),
      even_column_interchannelCoupling_next_column_(
          conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm")
              .getParameter<double>("Even_column_interchannelCoupling_next_column")),
      apply_timewalk_(conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<bool>("ApplyTimewalk")),
      timewalk_model_(
          conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<edm::ParameterSet>("TimewalkModel")),
      fedCablingMapToken_(iC.esConsumes()),
      geomToken_(iC.esConsumes()) {
      // Waveform configuration
      waveformModelEnabled_ = conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").exists("waveformModelEnabled") ? conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<bool>("waveformModelEnabled") : false;
      Krummenacher_ = conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").exists("Krummenacher") ? conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<double>("Krummenacher") : 0.0;
      riseTimeSignal_ = conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").exists("riseTimeSignal") ? conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<double>("riseTimeSignal") : 0.0;
      phase_ = conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").exists("phase") ? conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<double>("phase") : 0.0;
      asynchronous_ = conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").exists("asynchronous") ? conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<bool>("asynchronous") : false;
      timewindow_ = conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").exists("timewindow") ? conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<double>("timewindow") : 42;
      ToT80_ = conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").exists("ToT80") ? conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<bool>("ToT80") : false;
      ChosenBX_ = conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").exists("ChosenBX") ? conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<int>("ChosenBX") : 0;
  if (use_deadmodule_DB_)
    siPixelBadModuleToken_ = iC.esConsumes();
  if (use_LorentzAngle_DB_)
    siPixelLorentzAngleToken_ = iC.esConsumes();
  pixelFlag_ = true;
  LogDebug("PixelDigitizerAlgorithm") << "Algorithm constructed "
                                      << "Configuration parameters:"
                                      << "Threshold/Gain = "
                                      << "threshold in electron Endcap = " << theThresholdInE_Endcap_
                                      << "threshold in electron Barrel = " << theThresholdInE_Barrel_ << " "
                                      << theElectronPerADC_ << " " << theAdcFullScale_
                                      << " The delta cut-off is set to " << tMax_ << " pix-inefficiency "
                                      << addPixelInefficiency_;
}

PixelDigitizerAlgorithm::~PixelDigitizerAlgorithm() { LogDebug("PixelDigitizerAlgorithm") << "Algorithm deleted"; }

//
// -- Select the Hit for Digitization
//
bool PixelDigitizerAlgorithm::select_hit(const PSimHit& hit, double tCorr, double& sigScale) const {
  // in case of signal-shape emulation do not apply [TofLower,TofUpper] selection
  double toa = hit.tof() - tCorr;
  return apply_timewalk_ || (toa >= theTofLowerCut_ && toa < theTofUpperCut_);
}

// ======================================================================
//
//  Add  Cross-talk contribution
//
// ======================================================================
void PixelDigitizerAlgorithm::add_cross_talk(const Phase2TrackerGeomDetUnit* pixdet) {
  if (!pixelFlag_)
    return;

  const Phase2TrackerTopology* topol = &pixdet->specificTopology();

  // cross-talk calculation valid for the case of 25x100 pixels
  const float pitch_first = 0.0025;
  const float pitch_second = 0.0100;

  // 0.5 um tolerance when comparing the pitch to accommodate the small changes in different TK geometrie (temporary fix)
  const double pitch_tolerance(0.0005);

  if (std::abs(topol->pitch().first - pitch_first) > pitch_tolerance ||
      std::abs(topol->pitch().second - pitch_second) > pitch_tolerance)
    return;

  uint32_t detID = pixdet->geographicalId().rawId();
  signal_map_type& theSignal = _signal[detID];
  signal_map_type signalNew;

  int numRows = topol->nrows();
  int numColumns = topol->ncolumns();

  for (auto& s : theSignal) {
    float signalInElectrons = s.second.ampl();  // signal in electrons

    auto hitChan = PixelDigi::channelToPixel(s.first);

    float signalInElectrons_odd_row_Xtalk_next_row = signalInElectrons * odd_row_interchannelCoupling_next_row_;
    float signalInElectrons_even_row_Xtalk_next_row = signalInElectrons * even_row_interchannelCoupling_next_row_;
    float signalInElectrons_odd_column_Xtalk_next_column =
        signalInElectrons * odd_column_interchannelCoupling_next_column_;
    float signalInElectrons_even_column_Xtalk_next_column =
        signalInElectrons * even_column_interchannelCoupling_next_column_;

    // subtract the charge which will be shared
    s.second.set(signalInElectrons - signalInElectrons_odd_row_Xtalk_next_row -
                 signalInElectrons_even_row_Xtalk_next_row - signalInElectrons_odd_column_Xtalk_next_column -
                 signalInElectrons_even_column_Xtalk_next_column);

    if (hitChan.first != 0) {
      auto XtalkPrev = std::make_pair(hitChan.first - 1, hitChan.second);
      int chanXtalkPrev = pixelFlag_ ? PixelDigi::pixelToChannel(XtalkPrev.first, XtalkPrev.second)
                                     : Phase2TrackerDigi::pixelToChannel(XtalkPrev.first, XtalkPrev.second);
      if (hitChan.first % 2 == 1)
        signalNew.emplace(chanXtalkPrev,
                          digitizerUtility::Ph2Amplitude(signalInElectrons_even_row_Xtalk_next_row, nullptr, -1.0));
      else
        signalNew.emplace(chanXtalkPrev,
                          digitizerUtility::Ph2Amplitude(signalInElectrons_odd_row_Xtalk_next_row, nullptr, -1.0));
    }
    if (hitChan.first < numRows - 1) {
      auto XtalkNext = std::make_pair(hitChan.first + 1, hitChan.second);
      int chanXtalkNext = pixelFlag_ ? PixelDigi::pixelToChannel(XtalkNext.first, XtalkNext.second)
                                     : Phase2TrackerDigi::pixelToChannel(XtalkNext.first, XtalkNext.second);
      if (hitChan.first % 2 == 1)
        signalNew.emplace(chanXtalkNext,
                          digitizerUtility::Ph2Amplitude(signalInElectrons_odd_row_Xtalk_next_row, nullptr, -1.0));
      else
        signalNew.emplace(chanXtalkNext,
                          digitizerUtility::Ph2Amplitude(signalInElectrons_even_row_Xtalk_next_row, nullptr, -1.0));
    }

    if (hitChan.second != 0) {
      auto XtalkPrev = std::make_pair(hitChan.first, hitChan.second - 1);
      int chanXtalkPrev = pixelFlag_ ? PixelDigi::pixelToChannel(XtalkPrev.first, XtalkPrev.second)
                                     : Phase2TrackerDigi::pixelToChannel(XtalkPrev.first, XtalkPrev.second);
      if (hitChan.second % 2 == 1)
        signalNew.emplace(
            chanXtalkPrev,
            digitizerUtility::Ph2Amplitude(signalInElectrons_even_column_Xtalk_next_column, nullptr, -1.0));
      else
        signalNew.emplace(
            chanXtalkPrev,
            digitizerUtility::Ph2Amplitude(signalInElectrons_odd_column_Xtalk_next_column, nullptr, -1.0));
    }
    if (hitChan.second < numColumns - 1) {
      auto XtalkNext = std::make_pair(hitChan.first, hitChan.second + 1);
      int chanXtalkNext = pixelFlag_ ? PixelDigi::pixelToChannel(XtalkNext.first, XtalkNext.second)
                                     : Phase2TrackerDigi::pixelToChannel(XtalkNext.first, XtalkNext.second);
      if (hitChan.second % 2 == 1)
        signalNew.emplace(
            chanXtalkNext,
            digitizerUtility::Ph2Amplitude(signalInElectrons_odd_column_Xtalk_next_column, nullptr, -1.0));
      else
        signalNew.emplace(
            chanXtalkNext,
            digitizerUtility::Ph2Amplitude(signalInElectrons_even_column_Xtalk_next_column, nullptr, -1.0));
    }
  }
  for (auto const& l : signalNew) {
    int chan = l.first;
    auto iter = theSignal.find(chan);
    if (iter != theSignal.end()) {
      iter->second += l.second.ampl();
    } else {
      theSignal.emplace(chan, digitizerUtility::Ph2Amplitude(l.second.ampl(), nullptr, -1.0));
    }
  }
}

PixelDigitizerAlgorithm::TimewalkCurve::TimewalkCurve(const edm::ParameterSet& pset)
    : x_(pset.getParameter<std::vector<double>>("charge")), y_(pset.getParameter<std::vector<double>>("delay")) {
  if (x_.size() != y_.size())
    throw cms::Exception("Configuration")
        << "Timewalk model error: the number of charge values does not match the number of delay values!";
}

double PixelDigitizerAlgorithm::TimewalkCurve::operator()(double x) const {
  auto it = std::lower_bound(x_.begin(), x_.end(), x);
  if (it == x_.begin())
    return y_.front();
  if (it == x_.end())
    return y_.back();
  int index = std::distance(x_.begin(), it);
  double x_high = *it;
  double x_low = *(--it);
  double p = (x - x_low) / (x_high - x_low);
  return p * y_[index] + (1 - p) * y_[index - 1];
}

PixelDigitizerAlgorithm::TimewalkModel::TimewalkModel(const edm::ParameterSet& pset) {
  threshold_values = pset.getParameter<std::vector<double>>("ThresholdValues");
  const auto& curve_psetvec = pset.getParameter<std::vector<edm::ParameterSet>>("Curves");
  if (threshold_values.size() != curve_psetvec.size())
    throw cms::Exception("Configuration")
        << "Timewalk model error: the number of threshold values does not match the number of curves.";
  for (const auto& curve_pset : curve_psetvec)
    curves.emplace_back(curve_pset);
}

double PixelDigitizerAlgorithm::TimewalkModel::operator()(double q_in, double q_threshold) const {
  auto index = find_closest_index(threshold_values, q_threshold);
  return curves[index](q_in);
}

std::size_t PixelDigitizerAlgorithm::TimewalkModel::find_closest_index(const std::vector<double>& vec,
                                                                       double value) const {
  auto it = std::lower_bound(vec.begin(), vec.end(), value);
  if (it == vec.begin())
    return 0;
  else if (it == vec.end())
    return vec.size() - 1;
  else {
    auto it_upper = it;
    auto it_lower = --it;
    auto closest = (value - *it_lower > *it_upper - value) ? it_upper : it_lower;
    return std::distance(vec.begin(), closest);
  }
}
//
// -- Compare Signal with Threshold
//
bool PixelDigitizerAlgorithm::isAboveThreshold(const digitizerUtility::SimHitInfo* hitInfo,
                                               float charge,
                                               float thr) const {
  if (charge < thr)
    return false;
  if (apply_timewalk_ && hitInfo) {
    float corrected_time = hitInfo->time();
    double time = corrected_time + timewalk_model_(charge, thr);
    return (time >= theTofLowerCut_ && time < theTofUpperCut_);
  } else
    return true;
}
//
// -- Read Bad Channels from the Condidion DB and kill channels/module accordingly
//
void PixelDigitizerAlgorithm::module_killing_DB(const Phase2TrackerGeomDetUnit* pixdet) {
  bool isbad = false;
  uint32_t detID = pixdet->geographicalId().rawId();
  int ncol = pixdet->specificTopology().ncolumns();
  if (ncol < 0)
    return;
  std::vector<SiPixelQuality::disabledModuleType> disabledModules = siPixelBadModule_->getBadComponentList();

  SiPixelQuality::disabledModuleType badmodule;
  for (const auto& mod : disabledModules) {
    if (detID == mod.DetID) {
      isbad = true;
      badmodule = mod;
      break;
    }
  }

  if (!isbad)
    return;

  signal_map_type& theSignal = _signal[detID];  // check validity
  if (badmodule.errorType == 0) {               // this is a whole dead module.
    for (auto& s : theSignal)
      s.second.set(0.);  // reset amplitude
  } else {               // all other module types: half-modules and single ROCs.
    // Get Bad ROC position:
    // follow the example of getBadRocPositions in CondFormats/SiPixelObjects/src/SiPixelQuality.cc
    std::vector<GlobalPixel> badrocpositions;
    for (size_t j = 0; j < static_cast<size_t>(ncol); j++) {
      if (siPixelBadModule_->IsRocBad(detID, j)) {
        std::vector<CablingPathToDetUnit> path = fedCablingMap_->pathToDetUnit(detID);
        for (auto const& p : path) {
          const PixelROC* myroc = fedCablingMap_->findItem(p);
          if (myroc->idInDetUnit() == j) {
            LocalPixel::RocRowCol local = {39, 25};  //corresponding to center of ROC row, col
            GlobalPixel global = myroc->toGlobal(LocalPixel(local));
            badrocpositions.push_back(global);
            break;
          }
        }
      }
    }

    for (auto& s : theSignal) {
      std::pair<int, int> ip;
      if (pixelFlag_)
        ip = PixelDigi::channelToPixel(s.first);
      else
        ip = Phase2TrackerDigi::channelToPixel(s.first);

      for (auto const& p : badrocpositions) {
        for (auto& k : badPixels_) {
          if (p.row == k.getParameter<int>("row") && ip.first == k.getParameter<int>("row") &&
              std::abs(ip.second - p.col) < k.getParameter<int>("col")) {
            s.second.set(0.);
          }
        }
      }
    }
  }
}

void PixelDigitizerAlgorithm::digitize(const Phase2TrackerGeomDetUnit* pixdet,
                                      std::map<int, digitizerUtility::DigiSimInfo>& digi_map,
                                      const TrackerTopology* tTopo) {
  //bool waveformModelEnabled_ = false;
  //throw cms::Exception("LogicError") << "I made this!!!! \n";
  if (!waveformModelEnabled_) {
     return Phase2TrackerDigitizerAlgorithm::digitize(pixdet, digi_map, tTopo);
  }
  uint32_t detID = pixdet->geographicalId().rawId();
  auto it = _signal.find(detID);
  if (it == _signal.end())
    return;

  const signal_map_type& theSignal = _signal[detID];

  uint32_t Sub_detid = DetId(detID).subdetId();

  float theThresholdInE = 0.;
  float theHIPThresholdInE = 0.;
  // Define Threshold
  /* if (Sub_detid == PixelSubdetector::PixelBarrel || Sub_detid == StripSubdetector::TOB) {  // Barrel modules
    theThresholdInE = addThresholdSmearing_ ? smearedThreshold_Barrel_->fire()             // gaussian smearing
                                            : theThresholdInE_Barrel_;                     // no smearing
    theHIPThresholdInE = theHIPThresholdInE_Barrel_;
  } else {                                                                      // Forward disks modules
    theThresholdInE = addThresholdSmearing_ ? smearedThreshold_Endcap_->fire()  // gaussian smearing
                                            : theThresholdInE_Endcap_;          // no smearing
    theHIPThresholdInE = theHIPThresholdInE_Endcap_;
  } */

  if (Sub_detid == PixelSubdetector::PixelBarrel ) { // || Sub_detid == StripSubdetector::TOB) {
    theThresholdInE = theThresholdInE_Barrel_;
  } else {
    theThresholdInE = theThresholdInE_Endcap_;
  }
 
  //  if (addNoise) add_noise(pixdet, theThresholdInE/theNoiseInElectrons_);  // generate noise
  if (addNoise_)
    add_noise(pixdet);  // generate noise
  if (addXtalk_)
    add_cross_talk(pixdet);
  if (addNoisyPixels_) {
    float thresholdInNoiseUnits = 99.9;
    if (theNoiseInElectrons_)
      thresholdInNoiseUnits = theThresholdInE / theNoiseInElectrons_;
    add_noisy_cells(pixdet, thresholdInNoiseUnits);
  }

  // Do only if needed
  if (addPixelInefficiency_ && !theSignal.empty()) {
    if (use_ineff_from_db_)
      pixel_inefficiency_db(detID);
    else
      pixel_inefficiency(subdetEfficiencies_, pixdet, tTopo);
  }
  if (use_module_killing_) {
    if (use_deadmodule_DB_)  // remove dead modules using DB
      module_killing_DB(pixdet);
    else  // remove dead modules using the list in cfg file
      module_killing_conf(detID);
  }

  /* int counter = 0;
  std::ostringstream outputString; */
  
  std::pair<float, float> SignalPeak = CalculateSignalPeak(theThresholdInE);
  float waveModelThreshold = SignalPeak.first;
  float t_peak = SignalPeak.second;

  /* outputString << "     Configuration parameters:\n" << "waveformModelEnabled_: " << waveformModelEnabled_
               << "\nKrummenacher_: " << Krummenacher_ << "\nriseTimeSignal_: " << riseTimeSignal_ << "\nphase_ :" << phase_
               << "\nasynchronous_: " << asynchronous_ << "\ntimewindow_: " << timewindow_ << "\nthePhase2ReadoutMode_: " <<thePhase2ReadoutMode_; */
 
  //std::ostringstream outputStringTiming;

  // Digitize if the signal is greater than threshold
  for (auto const& s : theSignal) {
    const digitizerUtility::Ph2Amplitude& sig_data = s.second;
    float signalInElectrons = sig_data.ampl();

    const auto& info_list = sig_data.simInfoList();
    const digitizerUtility::SimHitInfo* hitInfo = nullptr;
    if (!info_list.empty())
      hitInfo = std::max_element(info_list.begin(), info_list.end())->second.get();
    
    if (!CoarseFiltering(signalInElectrons, t_peak, hitInfo -> time(), ChosenBX_, theThresholdInE)) {
      continue;
    }
    //auto start = std::chrono::high_resolution_clock::now();
    std::pair<float, float> times = crossThresholdTimes(signalInElectrons, waveModelThreshold);
    //auto time1 = std::chrono::high_resolution_clock::now();
    int AssignedBX = CalculateAssignedBX(times.first, times.second, hitInfo -> time(), signalInElectrons);
    //auto time2 = std::chrono::high_resolution_clock::now();
    //outputString << "  Signal: " << signalInElectrons << " t1: " << times.first << " t2: " << times.second << " hitInfo -> time(): " << hitInfo -> time() << " Assigned BX: " << AssignedBX;
    if (AssignedBX == ChosenBX_) {
      int ToT = convertSignalToADCWaveform(times.first, times.second, hitInfo -> time(), AssignedBX);
      digitizerUtility::DigiSimInfo info;
      info.sig_tot = ToT;  // adc
      info.ot_bit = signalInElectrons > theHIPThresholdInE ? true : false;
      if (makeDigiSimLinks_) {
        for (auto const& l : sig_data.simInfoList()) {
          float charge_frac = l.first / signalInElectrons;
          if (l.first > -5.0)
            info.simInfoList.push_back({charge_frac, l.second.get()});
        }
      }
      //outputString << ", ToT: " << ToT << "\n";
      digi_map.insert({s.first, info});
    } //else {outputString << "\n";}
    //throw cms::Exception("LogicError") << "This is the new version";
    /*auto time3 = std::chrono::high_resolution_clock::now();

    auto duration1 = std::chrono::duration_cast<std::chrono::nanoseconds>(time1 - start);
    auto duration2 = std::chrono::duration_cast<std::chrono::nanoseconds>(time2 - time1);
    auto duration3 = std::chrono::duration_cast<std::chrono::nanoseconds>(time3 - time2);
    
    outputStringTiming << "crossThresholdTimes: " << duration1.count() << " CalculateAssignedBX: " << duration2.count() << " ToT converter: " << duration3.count() << std::endl;*/
    /* ++counter;
    if (counter >= 32) {
      throw cms::Exception("LogicError") << "I made this!!!! \n" << outputString.str();
    } */
  }
  /* std::ofstream outFile("DebugHits.txt", std::ios::app);
  outFile << outputString.str();
  outFile.close(); */
}


// =========================== Waveform model =========================================
// Assumes the signal for a CROC can be approximated to this signal: 
//            S(t, Q) = Q(1 - exp{-t/tau}) - (I_k / 2) t, 
// Where tau is the rise time and I_k is the krummenacher current.
// By finding the times the signal crosses the threshold, the signal is converted into ToT
// To activate this model "waveformModelEnabled" must be set to True in the python configuration
// Important parameters:
//    Krummenacher Current : The current used for discharging the signal
//    Rise time            : The rise time of the signal
//    Phase                : The relation between corrTime and sampling edge
//    Async or sync        : All signals above threshold is seen in async mode, just the signals that are above threshold on a sampling edge is seen in Sync mode
//    Time window          : The lenght of time window in ns
//    ToT 80               : If enabled, counts ToT with twice the precision on the falling edge
//    Chosen BX            : Decides which assigned BX that should be digitized and saved to file.


// Signal peak, used for finding the waveform threshold and maximum timewalk
std::pair<float, float> PixelDigitizerAlgorithm::CalculateSignalPeak(float charge) {
    float t_peak = riseTimeSignal_ * log(2 * charge / (Krummenacher_ * riseTimeSignal_));
    float S_peak = charge * (1.0 - exp(- t_peak / riseTimeSignal_)) - 0.5 * Krummenacher_ * t_peak;
    return {S_peak, t_peak};
}

// Calculate the times the signal crosses the threshold
std::pair<float, float> PixelDigitizerAlgorithm::crossThresholdTimes(float charge, float thr) {
    float t1 = std::numeric_limits<float>::quiet_NaN();
    float t2 = std::numeric_limits<float>::quiet_NaN();

    long double B = 0.5 * Krummenacher_;
    long double A = charge - thr;

    long double z = - (charge / (B * riseTimeSignal_)) * std::exp(-A / (B * riseTimeSignal_));

    // Using Lambert W function to solve for times
    // W_0 branch, falling edge
    long double w0 = boost::math::lambert_w0(z);
    if (std::isfinite(w0) && A / B + riseTimeSignal_ * w0 > 0) {
        t2 = A / B + riseTimeSignal_ * w0;
    }

    // W_-1 branch, rising edge
    if (charge > B * riseTimeSignal_ * 11000) {  // If too high charge, we use an approximation (11000 is found be testing)
      t1 = - riseTimeSignal_ * log(1 - thr / charge);
    } else {
      long double w1 = boost::math::lambert_wm1(z);
      if (std::isfinite(w1) && A / B + riseTimeSignal_ * w1 > 0) {
          t1 = A / B + riseTimeSignal_ * w1;
      }
    }
    return {t1, t2};
}

// Calculate the assigned BX based on threshold crossing times
int PixelDigitizerAlgorithm::CalculateAssignedBX(float t1, float t2, float corrTime, float charge) {
    int AssignedBX = std::floor((t1 + corrTime - phase_) / timewindow_);
    if (asynchronous_) {
        return AssignedBX;
    } else {
        return ((t2 + corrTime) > ((AssignedBX + 1) * timewindow_ + phase_)) ? AssignedBX : - 1000;
    }
}

// Convert charge to ToT using waveform model
int PixelDigitizerAlgorithm::convertSignalToADCWaveform(float t1, float t2, float corrTime, int AssignedBX) {
    int signal_in_adc;
    int temp_signal;
    const int max_limit = 10;

    if (thePhase2ReadoutMode_ == 0) {
      temp_signal = theAdcFullScale_;
    } else {
      float denom = ToT80_ ? (timewindow_ / 2.0) : timewindow_;
      float samplingTime = (AssignedBX + 1) * timewindow_ + phase_;
      temp_signal = std::floor((t2 + corrTime - samplingTime) / denom);
      //throw cms::Exception("LogicError") << "I made this!!!! \n" << "Sampling Edge: " << samplingTime << " t2: " << t2 << " corrTime: " << corrTime << " denom: " << denom;
      if (thePhase2ReadoutMode_ != - 1) {
        // calculate the kink point and the slope
        int dualslope_param = std::min(std::abs(thePhase2ReadoutMode_), max_limit);
        int kink_point = static_cast<int>(theAdcFullScale_ / 2) + 1;
        // C-ROC: first valid ToT code above threshold is 0b0000
        if (temp_signal > kink_point) {
          temp_signal = std::floor((temp_signal - kink_point) / (pow(2, dualslope_param - 1))) + kink_point;
        }
      }
    }
    signal_in_adc = std::min(temp_signal, theAdcFullScale_);
    return signal_in_adc;
}   


bool PixelDigitizerAlgorithm::CoarseFiltering(float signalInElectrons, float t_peak, float corrTime, int ChosenBX, float thresholdINElectrons) {
    // A coarse filtering to avoid the computation of crossing times.
    // Check if the signal is inside the chosen BX with maximum timewalk and with no timewalk
    // Check if the signal is above the threshold
    return (corrTime + t_peak > ChosenBX * timewindow_ + phase_) && (corrTime < (ChosenBX + 1) * timewindow_ + phase_) && signalInElectrons >= thresholdINElectrons;
} 
