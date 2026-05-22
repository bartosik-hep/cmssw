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

#include "CLHEP/Random/RandGaussQ.h"

using namespace edm;
using namespace sipixelobjects;

void PixelDigitizerAlgorithm::init(const edm::EventSetup& es) {
  if (use_ineff_from_db_)  // load gain calibration service fromdb...
    theSiPixelGainCalibrationService_->setESObjects(es);

  if (use_deadmodule_DB_)
    siPixelBadModule_ = &es.getData(siPixelBadModuleToken_);

  if (use_LorentzAngle_DB_)  // Get Lorentz angle from DB record
    siPixelLorentzAngle_ = &es.getData(siPixelLorentzAngleToken_);

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
      waveformModelEnabled_(conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<bool>("waveformModelEnabled")),
      WaveformModel_(conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<edm::ParameterSet>("WaveFormModel")),

      apply_timewalk_(conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<bool>("ApplyTimewalk")),
      timewalk_model_(
          conf.getParameter<ParameterSet>("PixelDigitizerAlgorithm").getParameter<edm::ParameterSet>("TimewalkModel")),
      geomToken_(iC.esConsumes()) {
      
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

  if (waveformModelEnabled_) {
    uint32_t Sub_detid = DetId(hit.detUnitId()).subdetId(); 
    float LowerBound;
    if (Sub_detid == PixelSubdetector::PixelBarrel) { // Barrel modules  // || Sub_detid == StripSubdetector::TOB
      LowerBound = addThresholdSmearing_ ? theThresholdInE_Barrel_ - 2 * theThresholdSmearing_Barrel_           // 2 Sigma lower
                                          : theThresholdInE_Barrel_; // no smearing
    } else { // Forward disks modules
      LowerBound = addThresholdSmearing_ ? theThresholdInE_Endcap_ - 2 * theThresholdSmearing_Endcap_           // 2 Sigma lower
                                          : theThresholdInE_Endcap_; // no smearing
    } 

    WaveformModel::EffectiveParams p = WaveformModel_.computeEffectiveParams(LowerBound, WaveformModel_.nominalKrummenacher_);
    WaveformModel::SignalPeak signalPeak = WaveformModel_.calculateSignalPeak(LowerBound, p);

    return WaveformModel_.coarseTimeFiltering(signalPeak.time, float(toa));
    
  } else if (apply_timewalk_) {
    return true;
  } else {
    return (toa >= theTofLowerCut_ && toa < theTofUpperCut_);
  }
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
  throw cms::Exception("PixelDigitizerAlgorithm") << "Trying to kill modules from the pixel digitizer."
                                                  << " This method is not yet implemented!";
}

void PixelDigitizerAlgorithm::digitize(const Phase2TrackerGeomDetUnit* pixdet,
                                      std::map<int, digitizerUtility::DigiSimInfo>& digi_map,
                                      const TrackerTopology* tTopo) {
  if (!waveformModelEnabled_) {
     Phase2TrackerDigitizerAlgorithm::digitize(pixdet, digi_map, tTopo);
     return;
  }
  //throw cms::Exception("LogicError") << " I made this!!!!";
  uint32_t detID = pixdet->geographicalId().rawId();
  auto it = _signal.find(detID);
  if (it == _signal.end())
    return;

  const signal_map_type& theSignal = _signal[detID];

  uint32_t Sub_detid = DetId(detID).subdetId();

  float theThresholdInE = 0.;
  float theHIPThresholdInE = 0.;
  // Define Threshold
  if (Sub_detid == PixelSubdetector::PixelBarrel) {  // Barrel modules  // || Sub_detid == StripSubdetector::TOB
    theThresholdInE = addThresholdSmearing_ ? smearedThreshold_Barrel_->fire()             // gaussian smearing
                                            : theThresholdInE_Barrel_;                     // no smearing
    theHIPThresholdInE = theHIPThresholdInE_Barrel_;
  } else {                                                                      // Forward disks modules
    theThresholdInE = addThresholdSmearing_ ? smearedThreshold_Endcap_->fire()  // gaussian smearing
                                            : theThresholdInE_Endcap_;          // no smearing
    theHIPThresholdInE = theHIPThresholdInE_Endcap_;
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

  double sampledKrummenacher;
  if (WaveformModel_.addKrummSmearing_) {
    sampledKrummenacher = CLHEP::RandGaussQ::shoot(rengine_, WaveformModel_.nominalKrummenacher_, (WaveformModel_.krummSmearing_ / 100) * WaveformModel_.nominalKrummenacher_);
  } else {
    sampledKrummenacher = WaveformModel_.nominalKrummenacher_;
  }
  
  WaveformModel::EffectiveParams pThr = WaveformModel_.computeEffectiveParams(theThresholdInE, sampledKrummenacher);

  WaveformModel::SignalPeak signalPeak = WaveformModel_.calculateSignalPeak(theThresholdInE, pThr);
  float waveModelThreshold = signalPeak.amplitude;
  float t_peak = signalPeak.time;

  for (auto const& s : theSignal) {
    const digitizerUtility::Ph2Amplitude& sig_data = s.second;
    float signalInElectrons = sig_data.ampl();

    const auto& info_list = sig_data.simInfoList();
    const digitizerUtility::SimHitInfo* hitInfo = nullptr;
    int ToT;

    /* std::ofstream outFile("LogicTesting/TestWXtalk2.txt", std::ios::app);
    outFile << info_list.empty() << "\n";
    outFile.close(); */

    if (signalInElectrons < theThresholdInE) {
      continue;
    }

    if (info_list.empty()) { // No timing information, using baseline digitizer
      ToT = Phase2TrackerDigitizerAlgorithm::convertSignalToAdc(detID, signalInElectrons, theThresholdInE);
    } else { // Waveform model
      hitInfo = std::max_element(info_list.begin(), info_list.end())->second.get();
      float corrTime = hitInfo->time();

      if (!WaveformModel_.coarseTimeFiltering(t_peak, corrTime)) {
        continue;
      }

      WaveformModel::EffectiveParams p = WaveformModel_.computeEffectiveParams(signalInElectrons, sampledKrummenacher);

        std::ofstream outFile("LogicTesting/Krumm.txt", std::ios::app);
        outFile << "sampledKrummenacher: " << sampledKrummenacher << ", charge: " << signalInElectrons << ", param.krumm: " << p.krumm << "\n";
        outFile.close();


      std::pair<double, double> crossingTimes = WaveformModel_.crossThresholdTimes(signalInElectrons, waveModelThreshold, p);
      int assignedBX = WaveformModel_.calculateAssignedBX(crossingTimes.first, crossingTimes.second, corrTime, signalInElectrons);

      if (assignedBX != WaveformModel_.chosenBX_) {
        continue;
      }
      ToT = WaveformModel_.convertSignalToADCWaveform(crossingTimes.first, crossingTimes.second, corrTime, assignedBX, theAdcFullScale_, thePhase2ReadoutMode_);
    }

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
    digi_map.insert({s.first, info});
     
   
    // Debugging 
    //WaveformModel::SignalPeak signalPeak = WaveformModel_.calculateSignalPeak(signalInElectrons, p);
    std::ofstream outFile("LogicTesting/Test23rd.txt", std::ios::app);  // << ", Signal Peak: " << signalPeak.amplitude
    outFile << "Charge: " << signalInElectrons << ", Threshold: " << theThresholdInE << ", effective krumm: " << p.krumm << ", WaveformThr: " <<  waveModelThreshold <<
    //", Times: " << crossingTimes.first << ", " << crossingTimes.second << "\n";
    ", ToT: " << ToT << "\n";
    outFile.close();
  }
}


// =========================== Waveform model =========================================
// Assumes the signal for a CROC can be approximated to this signal: 
//            S(t, Q) = Q (1 - exp{-t/tau}) - I_k t, 
// Where tau is the rise time and I_k is the krummenacher current.
// By finding the times the signal crosses the threshold, the signal is converted into ToT
// To activate this model, "waveformModelEnabled" must be set to True in the python configuration
// Important parameters:
//    Krummenacher Current : The current used for discharging the signal
//    Rise time            : The rise time of the signal
//    Phase                : The relation between corrTime and sampling edge
//    Async or sync        : All signals above threshold is seen in async mode, just the signals that are above threshold on a sampling edge is seen in Sync mode
//    Time window          : The length of time window in ns
//    ToT 80               : If enabled, counts ToT with twice the precision on the falling edge
//    Chosen BX            : Decides which assigned BX that should be digitized and saved to file.

PixelDigitizerAlgorithm::WaveformModel::WaveformModel(const edm::ParameterSet& pset) {
  // Waveform parameters:
  //  Krummm:
  nominalKrummenacher_ = pset.getParameter<double>("Krummenacher");
  krummSatCharge_ = pset.getParameter<double>("KrummSatCharge");
  krummChargeOffset_ = pset.getParameter<double>("KrummChargeOffset");

  // Risetime:
  riseTimeSignal_ = pset.getParameter<double>("riseTimeSignal");
  riseTimeOffset_ = pset.getParameter<double>("riseTimeOffset");
  riseTimeSlope_ = pset.getParameter<double>("riseTimeSlope");

  // Fake t0, ground and gain
  timeOffset_ = pset.getParameter<double>("TimeOffset");
  groundOffset_ = pset.getParameter<double>("GroundOffset");
  signalOffset_ = pset.getParameter<double>("SignalOffset");

  // Krummenacher distribution
  addKrummSmearing_ = pset.getParameter<bool>("AddKrummSmearing");
  krummSmearing_ = pset.getParameter<double>("KrummSmearing");
  
  // Sampling parameters:
  phase_ = pset.getParameter<double>("phase");
  asynchronous_ = pset.getParameter<bool>("asynchronous");
  timewindow_ = pset.getParameter<double>("timewindow");
  tot80_ = pset.getParameter<bool>("ToT80");
  chosenBX_ = pset.getParameter<int>("ChosenBX");
}

PixelDigitizerAlgorithm::WaveformModel::EffectiveParams PixelDigitizerAlgorithm::WaveformModel::computeEffectiveParams(float charge, double sampledKrummenacher) const {
  EffectiveParams param;

  param.krumm = sampledKrummenacher * (1 - std::exp(- (charge -  krummChargeOffset_) / krummSatCharge_));

  param.risetime = charge <= riseTimeOffset_ ? riseTimeSignal_ : riseTimeSlope_ * (charge - riseTimeOffset_) + riseTimeSignal_;

  return param;
}

// Signal peak, used for finding the waveform threshold and maximum timewalk
PixelDigitizerAlgorithm::WaveformModel::SignalPeak PixelDigitizerAlgorithm::WaveformModel::calculateSignalPeak(float charge, EffectiveParams p) const {
    SignalPeak signalPeakResult;

    signalPeakResult.time = p.risetime * log(charge / (p.krumm * p.risetime));
    signalPeakResult.amplitude = charge * (1.0 - exp(- signalPeakResult.time / p.risetime)) - p.krumm * signalPeakResult.time;

    return signalPeakResult;
}

// Calculate the times the signal crosses the threshold
std::pair<double, double> PixelDigitizerAlgorithm::WaveformModel::crossThresholdTimes(float charge, float thr, EffectiveParams p) const {
    double effectiveCharge = charge + signalOffset_;
    double effectiveThr = thr + groundOffset_;

    double t1 = std::numeric_limits<double>::quiet_NaN();
    double t2 = std::numeric_limits<double>::quiet_NaN();

    double A = effectiveCharge - effectiveThr;

    double z = - (effectiveCharge / (p.krumm * p.risetime)) * std::exp(- A / (p.krumm * p.risetime));
    double edge = - 1 / std::exp(1);

    if (z <= edge) {
      std::ofstream outFile("LogicTesting/ShouldNeverHappen.txt", std::ios::app);
      outFile << z << "\n";
      outFile.close();
      return {t1, t2};
    }

    // Using Lambert W function to solve for times
    // W_0 branch, falling edge
    double w0 = boost::math::lambert_w0(z);
    if (std::isfinite(w0) && A / p.krumm + p.risetime * w0 > 0) {
        t2 = A / p.krumm + p.risetime * w0 - timeOffset_;
    }

    // W_-1 branch, rising edge
    if (effectiveCharge > p.krumm * p.risetime * 800) {  // If too high charge, we use an approximation (11000 is found be testing)
      t1 = - p.risetime * log(1 - effectiveThr / effectiveCharge) - timeOffset_;
    } else {
      double w1 = boost::math::lambert_wm1(z);
      if (std::isfinite(w1) && A / p.krumm + p.risetime * w1 > 0) {
          t1 = A / p.krumm + p.risetime * w1 - timeOffset_;
      }
    }
    return {t1, t2};
}

// Calculate the assigned BX based on threshold crossing times
int PixelDigitizerAlgorithm::WaveformModel::calculateAssignedBX(float t1, float t2, float corrTime, float charge) const {
    int AssignedBX = std::floor((t1 + corrTime - phase_) / timewindow_);
    if (asynchronous_) {
        return AssignedBX;
    } else {
        return ((t2 + corrTime) > ((AssignedBX + 1) * timewindow_ + phase_)) ? AssignedBX : - 1000;
    }
}

// Convert charge to ToT using waveform model
int PixelDigitizerAlgorithm::WaveformModel::convertSignalToADCWaveform(float t1, float t2, float corrTime, int AssignedBX, int theAdcFullScale, int thePhase2ReadoutMode) const {
    int signal_in_adc;
    int temp_signal;
    const int max_limit = 10;

    if (thePhase2ReadoutMode == 0) {
      temp_signal = theAdcFullScale;
    } else {
      float denom = tot80_ ? (timewindow_ / 2.0) : timewindow_;
      float samplingTime = (AssignedBX + 1) * timewindow_ + phase_;
      temp_signal = std::floor((t2 + corrTime - samplingTime) / denom);
      if (thePhase2ReadoutMode != - 1) {
        // calculate the kink point and the slope
        int dualslope_param = std::min(std::abs(thePhase2ReadoutMode), max_limit);
        int kink_point = static_cast<int>(theAdcFullScale / 2) + 1;
        // C-ROC: first valid ToT code above threshold is 0b0000
        if (temp_signal > kink_point) {
          temp_signal = std::floor((temp_signal - kink_point) / (pow(2, dualslope_param - 1))) + kink_point;
        }
      }
    }
    signal_in_adc = std::min(temp_signal, theAdcFullScale);
    return signal_in_adc;
}   


bool PixelDigitizerAlgorithm::WaveformModel::coarseTimeFiltering(float t_peak, float corrTime) const {
    // A coarse filtering to avoid the computation of crossing times.
    // Check if the signal is inside the chosen BX with maximum timewalk and with no timewalk
    return (corrTime + t_peak > chosenBX_ * timewindow_ + phase_) && (corrTime < (chosenBX_ + 1) * timewindow_ + phase_);
} 
