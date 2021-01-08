#ifndef RecoLocalCalo_EcalRecProducers_EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2_hh
#define RecoLocalCalo_EcalRecProducers_EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2_hh

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "RecoLocalCalo/EcalRecProducers/interface/EcalUncalibRecHitWorkerRunOneDigiBase.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitFixedAlphaBetaAlgo.h"

#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"
namespace edm {
  class Event;
  class EventSetup;
  class ParameterSet;
  class ParameterSetDescription;
}  // namespace edm

class EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2 final : public EcalUncalibRecHitWorkerRunOneDigiBase {
public:
  EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2(const edm::ParameterSet& ps, edm::ConsumesCollector&);
  EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2(){};
  ~EcalUncalibRecHitWorkerFixedAlphaBetaFitPh2(){};

  void set(const edm::EventSetup& es) override;
  bool run(const edm::Event& evt, const EcalDigiCollectionPh2::const_iterator& digi, EcalUncalibratedRecHitCollection& result);

  edm::ParameterSetDescription getAlgoDescription();

private:
  double AmplThrEB_;
  double AmplThrEE_;

  EcalUncalibRecHitFixedAlphaBetaAlgo<EBDataFrame> algoEB_;
  EcalUncalibRecHitFixedAlphaBetaAlgo<EEDataFrame> algoEE_;

  double alphaEB_;
  double betaEB_;
  double alphaEE_;
  double betaEE_;
  std::vector<std::vector<std::pair<double, double> > > alphaBetaValues_;  // List of alpha and Beta values [SM#][CRY#](alpha, beta)
  bool useAlphaBetaArray_;
  std::string alphabetaFilename_;

  bool setAlphaBeta();  // Sets the alphaBetaValues_ vectors by the values provided in alphabetaFilename_

  edm::ESHandle<EcalGainRatios> pRatio;
  edm::ESHandle<EcalPedestals> pedHandle;
  edm::ESGetToken<EcalGainRatios, EcalGainRatiosRcd> ratiosToken_;
  edm::ESGetToken<EcalPedestals, EcalPedestalsRcd> pedestalsToken_;
};
#endif
