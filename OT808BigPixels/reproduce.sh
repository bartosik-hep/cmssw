#!/bin/bash
# Workflow 38810 (FourMuExtendedPt1_200+Run4D129), same cmsDriver commands as runTheMatrix.py.
#   bash reproduce.sh           steps 1-2: aborts at event 60
#   STEP3=1 bash reproduce.sh   steps 1-3 (DQMIO)
set -e

NEVENTS=${NEVENTS:-100}
NTHREADS=${NTHREADS:-1}
SEED=${SEED:-6}
STEP3=${STEP3:-0}
WORKDIR=${WORKDIR:-$CMSSW_BASE/src/OT808BigPixels/ot808-bigpixels-seed$SEED}

SEEDCMD="from IOMC.RandomEngine.RandomServiceHelper import RandomNumberServiceHelper; randHelper=RandomNumberServiceHelper(process.RandomNumberGeneratorService); randHelper.resetSeeds($SEED)"

mkdir -p "$WORKDIR"
cd "$WORKDIR"

cmsDriver.py FourMuExtendedPt_1_200_pythia8_cfi -s GEN,SIM -n $NEVENTS \
  --conditions auto:phase2_realistic_T40_13TeV --beamspot NoSmear --datatier GEN-SIM --eventcontent FEVTDEBUG \
  --geometry ExtendedRun4D129 --era Phase2C17I13M9 --relval 10000,100 --fileout file:step1.root \
  --nThreads $NTHREADS --customise_commands "$SEEDCMD" > step1_FourMuExtendedPt1_200+Run4D129.log 2>&1

cmsDriver.py step2 -s DIGI:pdigi_valid,L1TrackTrigger,L1,L1P2GT,DIGI2RAW,HLT:@relvalRun4 \
  --conditions auto:phase2_realistic_T40 --datatier GEN-SIM-DIGI-RAW -n $NEVENTS --eventcontent FEVTDEBUGHLT \
  --geometry ExtendedRun4D129 --era Phase2C17I13M9 --filein file:step1.root --fileout file:step2.root \
  --nThreads $NTHREADS --customise_commands "$SEEDCMD" > step2_FourMuExtendedPt1_200+Run4D129.log 2>&1

[ "$STEP3" = 1 ] || exit 0

cmsDriver.py step3 -s RAW2DIGI,RECO,RECOSIM,PAT,VALIDATION:@phase2Validation+@miniAODValidation,DQM:@phase2+@miniAODDQM \
  --conditions auto:phase2_realistic_T40 --datatier DQMIO -n $NEVENTS --eventcontent DQM \
  --geometry ExtendedRun4D129 --era Phase2C17I13M9 --filein file:step2.root --fileout file:step3.root \
  --nThreads $NTHREADS --customise_commands "$SEEDCMD" > step3_FourMuExtendedPt1_200+Run4D129.log 2>&1
