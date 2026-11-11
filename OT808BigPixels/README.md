# Setting up the environment

```bash
# pwd = $CMSSW_BASE/src  (CMSSW_20_1_0_pre4)

git remote add gabe https://github.com/gabrielribcesario/cmssw.git

git fetch gabe

git switch -c local-ot808-bigpixels-reproducer gabe/ot808-bigpixels-reproducer

git cms-addpkg Configuration/AlCa Configuration/Geometry Configuration/PyReleaseValidation Configuration/StandardSequences \
               Geometry/TrackerCommonData Geometry/TrackerRecoData Geometry/CMSCommonData Geometry/TrackerSimData \
               L1Trigger/TrackFindingTracklet \
               tkLa-dev/DesignAlignment \
               OT808BigPixels

scram b -j 8
```

D129 is OT v8.0.8 + IT v7.5.0. OT v8.0.8 is like OT v8.0.7 but with PS module type `ptPSlarger_320_bigPixels`: the PS-p sensors declare `BigPixels_X = 2`, `BigPixels_Pitch_X = 0.02`, `BigPixels_TwoEdges_X = true`.

# MatchProcessor.cc

TL;DR: one stub on a tilted L2 PS module (row 476) has $|\Delta\varphi| = 1.30$ mrad, over the 1 mrad tolerance (`MatchProcessor.cc:519`).

Build without the "Add MatchProcessor debugging info" commit (e.g. `git revert --no-commit` it), then:

```bash
bash OT808BigPixels/reproduce.sh

grep Assertion OT808BigPixels/ot808-bigpixels-seed6/step2_FourMuExtendedPt1_200+Run4D129.log
```

Workflow 38810 (`FourMuExtendedPt1_200+Run4D129`) steps 1–2, 100 events, 1 thread, seed 6. Step 2 should abort at event 60, where the debug build below logs the stub:

```
MatchProcessor.cc:519: ... Assertion `std::abs(dphi) < 0.001' failed.
```

Hypothesis, from reading the code (not verified):

1. `RectangularPixelPhase2Topology` puts all 16 big pixels of the 8-ROC PS-p sensor in one central block (rows 470–485), instead of 2 at each edge of every MPA.
2. The float $\varphi$ (`tt::Setup::stubPos`, `topol->localPosition(mp)`) follows that topology.
3. The DTC $\varphi$ (`TrackerDTC/src/Stub.cc:28,41`) is linear in `(row - numRows/2) * pitch`.
4. Inside the block the two differ by `(row - 470) * 0.01 - 0.02` cm: 0.045 cm at row 476.5, i.e. 1.22 mrad at r = 37 cm against 1.30 mrad measured.

Debugging info: build with the "Add MatchProcessor debugging info" commit, which replaces the asserts with a `LogError`, rerun, then:

```bash
grep -A1 MatchProcessorDphi OT808BigPixels/ot808-bigpixels-seed6/step2_FourMuExtendedPt1_200+Run4D129.log
```

```
MP_L2PHIA dphi=-0.00129818 detId=438578185 layerdisk=1 tilted=1 tiltedRing=2 endcapRing=0 PS=1 flipped=1 strip=476 bend=-0.5 DTClink=negPS10G_3_A region=1 r=36.9811 z=-32.8999 phi=0.359787 rapprox=36.9434 zapprox=-32.9297 phiapprox=0.361085 phimin=0.18425 iphi=2819 nbitsphi=14 stubword=F9B9C5819
```

# TTStubAlgorithm_official.cc

TL;DR: compared with D121, ~99% of PS stubs are missing and PS-p clusters are halved. 2S stubs and PS-s clusters are unchanged.

With the debug build (or `SEED=1` without it):

```bash
STEP3=1 bash OT808BigPixels/reproduce.sh
```

Hypothesis for the stub loss, from reading the code (not verified). The PS-p cluster loss is unexplained.

1. `TTStubAlgorithm_official.cc:87` correlates clusters by raw row: `dispD = 2 * (mp1.x() - mp0.x()) * (pitch0.first / pitch1.first)`.
2. PS-p has 944 rows (16 big), PS-s has 960, both at 0.0100 cm.
3. For a straight track the strip − pixel row offset is 12 half-strips below the block and 44 above it.
4. Only tracks whose real bend cancels the offset form a stub.

Step 3 DQM vs. D121 (same workflow, 1000 events, 16 threads):

| | D121 | D129 |
|---|---|---|
| `OuterTracker/L1Stub_isPS` mean | 0.637 | **0.014** |
| Barrel L1 / L2 / L3 stubs / event | 3.05 / 2.31 / 1.79 | **0.00 / 0.03 / 0.02** |
| Barrel L4 / L5 / L6 stubs / event (2S) | 1.41 / 1.18 / 1.06 | 1.37 / 1.37 / 1.09 |
| PS-p clusters (HLT), barrel L1 / L2 / L3 / event | 3.58 / 2.58 / 2.02 | **1.65 / 1.41 / 1.04** |
| PS-s clusters (HLT), barrel L1 / L2 / L3 / event | 3.56 / 2.55 / 1.99 | 3.56 / 2.53 / 1.85 |

Stubs: `OuterTracker/Barrel/Layer{1..6}/Num_L1Stubs_Per_Event`. Clusters: `HLT/TrackerPhase2OTCluster/Barrel/Layer{1..3}/Num_Clusters_Per_Event_{P,S}`.
