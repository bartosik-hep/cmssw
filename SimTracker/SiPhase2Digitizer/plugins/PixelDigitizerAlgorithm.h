#ifndef _SimTracker_SiPhase2Digitizer_PixelDigitizerAlgorithm_h
#define _SimTracker_SiPhase2Digitizer_PixelDigitizerAlgorithm_h

#include "CondFormats/SiPixelObjects/interface/GlobalPixel.h"
#include "CondFormats/DataRecord/interface/SiPixelQualityRcd.h"
#include "CondFormats/DataRecord/interface/SiPixelLorentzAngleSimRcd.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "SimTracker/SiPhase2Digitizer/plugins/Phase2TrackerDigitizerAlgorithm.h"

class PixelDigitizerAlgorithm : public Phase2TrackerDigitizerAlgorithm {
private:
  // A list of 2d points
  class TimewalkCurve {
  public:
    // pset must contain "charge" and "delay" of type vdouble
    TimewalkCurve(const edm::ParameterSet& pset);

    // linear interpolation
    double operator()(double x) const;

  private:
    std::vector<double> x_;
    std::vector<double> y_;
  };

  // Holds the timewalk model data
  class TimewalkModel {
  public:
    TimewalkModel(const edm::ParameterSet& pset);

    // returns the delay for given input charge and threshold
    double operator()(double q_in, double q_threshold) const;

  private:
    std::size_t find_closest_index(const std::vector<double>& vec, double value) const;

    std::vector<double> threshold_values;
    std::vector<TimewalkCurve> curves;
  };

  // Holds the waveform parameters
  class WaveformModel {
  public:
    WaveformModel(const edm::ParameterSet& pset);

    struct SignalPeak {
        float amplitude;
        float time;
    };

    struct EffectiveParams {
        double krumm;
        double risetime;
      };

    EffectiveParams computeEffectiveParams(float charge, double sampledKrummenacher) const;

    SignalPeak calculateSignalPeak(float charge, EffectiveParams p) const;
    std::pair<double, double> crossThresholdTimes(float charge, float thr, EffectiveParams p) const;
    int calculateAssignedBX(float t1, float t2, float corrTime, float charge) const;
    int convertSignalToADCWaveform(float t1, float t2, float corrTime, int AssignedBX, int theAdcFullScale, int thePhase2ReadoutMode) const;
    bool coarseTimeFiltering(float t_peak, float corrTime) const;

    int chosenBX_;

    double nominalKrummenacher_;
    bool addKrummSmearing_;
    double krummSmearing_;


  private:
    double krummSatCharge_;
    double krummChargeOffset_;

    double riseTimeSignal_;
    double riseTimeOffset_;
    double riseTimeSlope_;

    double timeOffset_;
    double groundOffset_;
    double signalOffset_;

    double phase_; 
    bool asynchronous_;
    double timewindow_;
    bool tot80_;
  };

public:
  PixelDigitizerAlgorithm(const edm::ParameterSet& conf, edm::ConsumesCollector iC);
  ~PixelDigitizerAlgorithm() override;

  // initialization that cannot be done in the constructor
  void init(const edm::EventSetup& es) override;

  bool select_hit(const PSimHit& hit, double tCorr, double& sigScale) const override;
  bool isAboveThreshold(const digitizerUtility::SimHitInfo* hitInfo, float charge, float thr) const override;
  void add_cross_talk(const Phase2TrackerGeomDetUnit* pixdet) override;
  void module_killing_DB(const Phase2TrackerGeomDetUnit* pixdet) override;
  void digitize(const Phase2TrackerGeomDetUnit* pixdet, std::map<int, digitizerUtility::DigiSimInfo>& digi_map, const TrackerTopology* tTopo) override;

  // Addition four xtalk-related parameters to PixelDigitizerAlgorithm specific parameters initialized in Phase2TrackerDigitizerAlgorithm
  double odd_row_interchannelCoupling_next_row_;
  double even_row_interchannelCoupling_next_row_;
  double odd_column_interchannelCoupling_next_column_;
  double even_column_interchannelCoupling_next_column_;

  // Waveform parameteres
  bool waveformModelEnabled_;
  const WaveformModel WaveformModel_;

  // Timewalk parameters
  bool apply_timewalk_;
  const TimewalkModel timewalk_model_;

  edm::ESGetToken<SiPixelQuality, SiPixelQualityRcd> siPixelBadModuleToken_;
  edm::ESGetToken<SiPixelLorentzAngle, SiPixelLorentzAngleSimRcd> siPixelLorentzAngleToken_;
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> geomToken_;
};
#endif
