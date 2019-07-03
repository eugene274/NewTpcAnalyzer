Scripts to run TPC calibration software in SHINE

https://twiki.cern.ch/twiki/bin/viewauth/NA61/TpcDvCalibrationInstructions


### Building NewAnalyzer

To build analyzer, go to source location and type

```


SetShineEnv <...>
mkdir build; cd build
cmake ../

make
```

New analyzer is compatible (and tested) with
**v1r13p0/x86_64-centos7-gcc8-opt** ShineOffline environment.


### Phylosophy

The script is based on Stage2 algorithm described here

https://twiki.cern.ch/twiki/pub/NA61/TpcDvCalibrationInstructions/20149424_Calibration_chain_with_shine_offline.pdf

For the calibration we use track matching details between TPCs
(TOF-TPC):

time - Unix time

master_Y - track c-te at the upstream TPC

slave_Y - extrapolated track c-te of the downstream TPC to the plane of
upstream TPC


Non-zero slope '_S_' of dependency `dY = (slave_Y - master_Y)_ vs
_master_Y` indicates difference between true drift velocity and vD in
database.

Corrected vD is calculated using formula

``` 
vD(true) = vD(DB) * (1 + s)^-1
```

Since drift velocity is time-dependent quantity, measurement of slope is
performed in time bins. Time slice is formed with configurable fixed
number of tracks (2000 by default).

#### Offset/slope propagation

For `dY = (slave_Y - master_Y)` we believe upstream TPC represented by
_master_Y_ is calibrated, so _master_Y_ is true Y in TPC. To satisfy
this condition we estimate true 'Y' in upstream TPC.

```
Y_prim = Y - Y*slope - offset
```

where _slope_ and _offset_ are taken from previous step. This step
produces instability increasing from step to step: fluctuation of
slope/offset in MTPCL-TOFL will affect VTPC2-MTPCL and so on.

#### Smoothing (new)

To get rid of instability, time-smoothing procedure is introduced.
Before going to the downstream TPC, offset and slope vs time is smoothed
with Lowess algorithm with the specified time window (5 min for slope
and 1 hour for offset).


### Configuration (TODO)


### Running Stage2

To run calibration, type

` ./NewAnalyzer -i "<input files location>/trackMatchDump-*" -o
"vDCalibOutput.root"`

`vDCalibOutput.root` contains QA histograms dY vs Y in time slices,
TGraphs with offsets and slope:

MTPCL/grSlope - slope vs Time

MTPCL/grSlopeLowess - slope vs Time smoothed with Lowess algorithm

MTPCL/grOffset - offset vs Time

MTPCL/grOffsetBottom - offset vs Time in the bottom point of TPC

**MTPCL/grOffsetLowess** - offset vs Time in the bottom point of TPC
smoothed with Lowess algorithm

MTPCL/grRecVDrift - Drift velocity from the DB

**MTPCL/grCalibVDrift** - Drift Velocity after calibration

#### Notes








