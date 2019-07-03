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


### Algorithm




### Configuration


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








