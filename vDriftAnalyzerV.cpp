
#include <cassert>

#include "vDriftAnalyzerV.h"
#include <algorithm>
#include <iomanip>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TChain.h>
#include <TGraphSmooth.h>
#include <TSystem.h>
#include <TError.h>
#include <TFileCollection.h>

struct CalibTask_t {
  std::string tgtName{};
  ushort tgt{kDetPairs};
  ushort upstream{kDetPairs};

  /*** options ***/
  bool swap_{false};
  bool splitEvents_{false};
  int nEntries_{1000};

  int ny_{25};
  int ndy_{20};

  double slopeSmoothingSpan_{600};
  double offsetSmoothingSpan_{3600};

  double yBottom_{0};
  double yTop_{0};

  double xLo{-999};
  double xHi{999};

  /*** slice data ***/
  vDriftTreeStructure trackMatchData_{};

  TH2D *hdYvsY{nullptr};
  TProfile *pdYvsY{nullptr};

  TChain *sourceChain_{nullptr};

  /*** results ***/
  TGraph *recVDriftGraph{nullptr};
  TGraphErrors *calibVDriftGraph_{nullptr};

  TGraphErrors *offsetGraph_{nullptr};
  TGraphErrors *offsetBottomGraph_{nullptr};
  TGraphErrors *offsetBottomGraphSmooth_{nullptr};
  TGraphErrors *offsetTopGraph_{nullptr};
  TGraphErrors *offsetTopGraphSmooth_{nullptr};

  TGraphErrors *offsetXGraph_{nullptr};
  TGraphErrors *offsetXGraphSmooth_{nullptr};

  TGraphErrors *slopeGraph_{nullptr};
  TGraphErrors *slopeGraphSmooth_{nullptr};

  bool isDone_{false};

  /*** c-tors ***/
  CalibTask_t() = default;
  CalibTask_t(const std::string &tgtName, ushort tgt, ushort upstream = kDetPairs)
      : tgtName(tgtName), tgt(tgt), upstream(upstream) {}

  /*** builder ***/
  CalibTask_t &swap(bool swap = true) {
    this->swap_ = swap;
    return *this;
  }

  CalibTask_t &ne(int ne) {
    this->nEntries_ = ne;
    return *this;
  }

  CalibTask_t &ndy(int ndy) {
    this->ndy_ = ndy;
    return *this;
  }

  CalibTask_t &xLim(double lo, double hi) {
    this->xLo = lo;
    this->xHi = hi;
    return *this;
  }

  CalibTask_t &bottom(double yb) {
    this->yBottom_ = yb;
    return *this;
  }

  CalibTask_t &top(double yb) {
    this->yTop_ = yb;
    return *this;
  }

  CalibTask_t &smoSpanSlope(double span) {
    this->slopeSmoothingSpan_ = span;
    return *this;
  }

  /*** initialization ***/
  void initHistograms() {
    delete hdYvsY;
    delete pdYvsY;
    delete slopeGraph_;
    delete offsetGraph_;
    delete offsetXGraph_;

    float yLo = gSTATIC_INFO.at(tgt).yLimits.first;
    float yHi = gSTATIC_INFO.at(tgt).yLimits.second;
    hdYvsY = new TH2D("hdYvsY", ";Y (cm); dY (cm)", ny_, yLo, yHi, ndy_, -2., 2.);
    hdYvsY->SetDirectory(nullptr);
    pdYvsY = new TProfile("pdYvsY", ";Y (cm); dY (cm)", ny_, yLo, yHi);
    pdYvsY->SetDirectory(nullptr);

    recVDriftGraph = new TGraph;

    offsetGraph_ = new TGraphErrors;
    offsetBottomGraph_ = new TGraphErrors;
    offsetTopGraph_ = new TGraphErrors;

    offsetXGraph_ = new TGraphErrors;

    slopeGraph_ = new TGraphErrors;
  }

  void initInput(const std::string &fileName) {
    delete sourceChain_;

    sourceChain_ = new TChain(gSTATIC_INFO.at(tgt).treeName);

    if (endsWith(fileName, ".list")) {
      TFileCollection fc("fc", "", fileName.c_str());
      sourceChain_->AddFileInfoList(reinterpret_cast<TCollection *>(fc.GetList()));
    } else {
      sourceChain_->Add(fileName.c_str());
    }

    sourceChain_->ls();
    ReadBranchesFromTree(sourceChain_, trackMatchData_, swap_ ? "swap" : "");
  }

  /*** getters ***/
  bool hasUpstream() const { return upstream < kDetPairs; };

  static bool endsWith(const std::string &str, const std::string &suffix) {
    if (suffix.size() > str.size()) {
      return false;
    }

    if (suffix.empty()) {
      throw std::logic_error("suffix is empty");
    }

    std::string strEnd(str.end() - suffix.size(), str.end());
    return strEnd == suffix;
  }

};

int main(int argc, char **argv) {
  std::cout << "Starting analyser." << std::endl;

  std::string iFileName;
  std::string outputFileName{};
  TFile *outputFile;
  bool is_no_multiply;
  bool is_no_TOF_factor;
  // long int prev_time;
  std::string s_cut;
  is_no_multiply = false;
  is_no_TOF_factor = false;

  //-------------------------------USAGE-----------------------------------------//
  if (argc < 3) {
    Usage(argv);
    return error_code[0];
  }

  for (int i = 1; i < argc; i++) {
    if (
        std::string(argv[i]) != "-i" &&
            std::string(argv[i]) != "-o" &&
            std::string(argv[i]) != "--no-factor-multiplication" &&
            std::string(argv[i]) != "--no-TOF-MTPC-factor"
        ) {
      std::cerr << "\n[ERROR]: Unknown parameter " << i << ": " << argv[i] << std::endl;
      Usage(argv);
      std::cout << "Analyzer finished with exit code " << error_code[1] << std::endl;
      return error_code[1];
    } else {
      // input file
      if (std::string(argv[i]) == "-i" && i != argc - 1) {
        iFileName = argv[++i];
        continue;
      }
      if (std::string(argv[i]) == "-i" && i == argc - 1) {
        std::cerr << "\n[ERROR]: File name was not specified " << std::endl;
        Usage(argv);
        std::cout << "Analyzer finished with exit code " << error_code[2] << std::endl;
        return error_code[2];
      }

      // output file
      if (std::string(argv[i]) == "-o" && i != argc - 1) {
        if (i != argc - 1) {
          outputFileName = std::string(argv[++i]);
        } else {
          std::cerr << "\n[ERROR]: Output file name was not specified " << std::endl;
          Usage(argv);
          std::cout << "Analyzer finished with exit code " << error_code[2] << std::endl;
          return error_code[2];
        }
        continue;
      }

      if (std::string(argv[i]) == "--no-factor-multiplication") {
        is_no_multiply = true;
        continue;
      }
      if (std::string(argv[i]) == "--no-TOF-MTPC-factor") {
        is_no_TOF_factor = true;
        continue;
      }
    }
  }
  //-------------------------------USAGE-----------(end)-------------------------//

  if (outputFileName.empty()) {
    outputFileName = std::string(gSystem->DirName(iFileName.c_str())) + "/vdCalibOutput.root";
  }

  std::cout << "\nStarting options:" << std::endl;
  std::cout << "\nInput file        : " << iFileName << std::endl;
  std::cout << "\nOutput file        : " << outputFileName << std::endl;
  if (is_no_multiply) std::cout << "Multiplication factor : OFF - !!! Make sure this is intentional !!!" << std::endl;
  if (is_no_TOF_factor) std::cout << "MTPC-TOF factor       : OFF - !!! Make sure this is intentional !!!" << std::endl;
  std::cout << std::endl;

  outputFile = new TFile(outputFileName.c_str(), "recreate");

  bool usePropagation = true;

  std::vector<CalibTask_t> calibTasks{
      CalibTask_t("MTPCL", kMTPCLvsTOFL).ne(20000).smoSpanSlope(1800)
          .bottom(-60).top(60),
      CalibTask_t("VTPC2", kVTPC2vsMTPCL, kMTPCLvsTOFL).ne(2000).ndy(40).xLim(-110, 60)
          .bottom(-55).top(55),
      CalibTask_t("VTPC1", kVTPC1vsVTPC2, kVTPC2vsMTPCL).ne(2000).ndy(40)
          .bottom(-35).top(35),
      CalibTask_t("MTPCR", kVTPC2vsMTPCR, kVTPC2vsMTPCL).ne(2000).swap(),
      CalibTask_t("MTPCRfromTOFR", kMTPCRvsTOFR).ne(20000).smoSpanSlope(1800),
  };

  // TODO Order calibTasks
  std::vector<CalibTask_t> calibTasksOrdered{calibTasks};

  std::cout << "Tasks: ";
  for (const auto &task : calibTasksOrdered) {
    std::cout << gSTATIC_INFO.at(task.tgt).treeName << (task.swap_ ? " (swap)" : "") << " -> ";
  }
  std::cout << "End" << std::endl;

  for (auto &calibTask : calibTasksOrdered) {
    calibTask.initHistograms();
    calibTask.initInput(iFileName);
    auto &trackMatchData = calibTask.trackMatchData_;

    bool hasUpstream = calibTask.hasUpstream();
    CalibTask_t *upstreamTask = nullptr;
    for (auto &tt : calibTasksOrdered) {
      if (tt.tgt == calibTask.upstream) {
        upstreamTask = &tt;
        break;
      }
    }

    TDirectory *outDir = outputFile->mkdir(calibTask.tgtName.c_str(), "");
    TDirectory *qaDir = outDir->mkdir("qa");

    Info(__func__, "Processing %s", gSTATIC_INFO.at(calibTask.tgt).treeName);

    long nEntries = calibTask.sourceChain_->GetEntries();
    long iEntry = 0;

    int iSlice = 0;
    int iSliceEntries = 0;

    long currentEventNumber = -1;
    long currentUnixTime = -1;

    long sliceStartUnixTime = -1;

    int iRecVD = 0;
    while (iEntry < nEntries) {
      calibTask.sourceChain_->GetEntry(iEntry);

      if (iSliceEntries == 0) {
        calibTask.pdYvsY->Reset();
        calibTask.hdYvsY->Reset();

        sliceStartUnixTime = trackMatchData.eventUnixTime;
      }

      long unixTime = trackMatchData.eventUnixTime;
      float dy = trackMatchData.slave_Y - trackMatchData.master_Y;
      float y = trackMatchData.master_Y;
      float recVDrift = trackMatchData.slave_recVDrift;

      if (iRecVD == 0 || unixTime - calibTask.recVDriftGraph->GetX()[iRecVD - 1] > 5.) {
        calibTask.recVDriftGraph->SetPoint(iRecVD, unixTime, 1000 * recVDrift);
        ++iRecVD;
      }

      if (usePropagation && upstreamTask) {
        assert(upstreamTask->isDone_);
        double tofFactor = gSTATIC_INFO.at(upstreamTask->tgt).tofFactor;
        double upstreamSlope = tofFactor * upstreamTask->slopeGraphSmooth_->Eval(unixTime);
        double upstreamOffset = tofFactor * upstreamTask->offsetTopGraphSmooth_->Eval(unixTime);

        dy += upstreamSlope * y + upstreamOffset;
      }

      if (abs(dy) < 2. && (calibTask.xLo < trackMatchData.master_X && trackMatchData.master_X < calibTask.xHi)) {
        calibTask.hdYvsY->Fill(y, dy);
        calibTask.pdYvsY->Fill(y, dy);

        ++iSliceEntries;
      }

      if (iSliceEntries > calibTask.nEntries_) {
        if (calibTask.splitEvents_ || (currentEventNumber != trackMatchData.eventNumber)) {

          long sliceEndUnixTime = trackMatchData.eventUnixTime;
          long dT = sliceEndUnixTime - sliceStartUnixTime;
          double sliceT = (sliceStartUnixTime + sliceEndUnixTime) / 2.;

          calibTask.pdYvsY->Fit("pol1", "Q", "");

          auto *fitFun = calibTask.pdYvsY->GetFunction("pol1");

          double offset = fitFun->GetParameter(0);
          double offsetError = fitFun->GetParError(0);

          double offsetBottom = fitFun->Eval(calibTask.yBottom_);
          double offsetBottomError = offsetError; // TODO
          double offsetTop = fitFun->Eval(calibTask.yTop_);
          double offsetTopError = offsetError; // TODO

          double slope = fitFun->GetParameter(1);
          double slopeError = fitFun->GetParError(1);

          calibTask.offsetGraph_->SetPoint(iSlice, sliceT, offset);
          calibTask.offsetGraph_->SetPointError(iSlice, 0., offsetError);

          calibTask.offsetBottomGraph_->SetPoint(iSlice, sliceT, offsetBottom);
          calibTask.offsetBottomGraph_->SetPointError(iSlice, 0., offsetBottomError);

          calibTask.offsetTopGraph_->SetPoint(iSlice, sliceT, offsetTop);
          calibTask.offsetTopGraph_->SetPointError(iSlice, 0., offsetTopError);

          calibTask.slopeGraph_->SetPoint(iSlice, sliceT, slope);
          calibTask.slopeGraph_->SetPointError(iSlice, 0., slopeError);

          qaDir->WriteObject(calibTask.hdYvsY, Form("hdYvsY_%d", iSlice));
          qaDir->WriteObject(calibTask.pdYvsY, Form("pdYvsY_%d", iSlice));

          Info(__func__, "Slice %d: %d entries", iSlice, iSliceEntries);
          Info(__func__, "T = %ld; dT = %ld sec", long(sliceT), dT);
          Info(__func__, "offset = %e +- %e (cm)", offset, offsetError);
          Info(__func__, "slope = %e +- %e (cm)", slope, slopeError);

          ++iSlice;
          iSliceEntries = 0;
        }
      }

      currentEventNumber = trackMatchData.eventNumber;
      currentUnixTime = trackMatchData.eventUnixTime;

      ++iEntry;
    }

    double tStart = calibTask.slopeGraph_->GetX()[0];
    double tEnd = calibTask.slopeGraph_->GetX()[calibTask.slopeGraph_->GetN() - 1];
    double slopeSpan = calibTask.slopeSmoothingSpan_ / (tEnd - tStart);
    double offsetSpan = calibTask.offsetSmoothingSpan_ / (tEnd - tStart);

    {
      TGraphSmooth graphSmoothManager;
      calibTask.slopeGraphSmooth_ =
          new TGraphErrors(*dynamic_cast<TGraphErrors *>(graphSmoothManager.SmoothLowess(calibTask.slopeGraph_,
                                                                                         "",
                                                                                         slopeSpan, 5)));
      outDir->WriteObject(calibTask.slopeGraph_, "grSlope");
      outDir->WriteObject(calibTask.slopeGraphSmooth_, "grSlopeLowess");

      calibTask.calibVDriftGraph_ = new TGraphErrors(calibTask.recVDriftGraph->GetN());
      for (int ip = 0; ip < calibTask.calibVDriftGraph_->GetN(); ++ip) {
        double unixTime{-1};
        double recVDrift{-1};
        calibTask.recVDriftGraph->GetPoint(ip, unixTime, recVDrift);
        double slope = calibTask.slopeGraphSmooth_->Eval(unixTime);
        double tofFactor = gSTATIC_INFO.at(calibTask.tgt).tofFactor;

        double calibVDrift = recVDrift * 1.0 / (1 + tofFactor * slope);
        calibTask.calibVDriftGraph_->SetPoint(ip, unixTime, calibVDrift);
      }
      outDir->WriteObject(calibTask.calibVDriftGraph_, "grCalibVDrift");
      outDir->WriteObject(calibTask.recVDriftGraph, "grRecVDrift");

    }

    {
      TGraphSmooth graphSmoothManager;
      calibTask.offsetBottomGraphSmooth_ =
          new TGraphErrors(*dynamic_cast<TGraphErrors *>(graphSmoothManager.SmoothLowess(calibTask.offsetBottomGraph_,
                                                                                         "",
                                                                                         offsetSpan, 5)));
      outDir->WriteObject(calibTask.offsetGraph_, "grOffset");
      outDir->WriteObject(calibTask.offsetBottomGraph_, "grOffsetBottom");
      outDir->WriteObject(calibTask.offsetBottomGraphSmooth_, "grOffsetBottomLowess");
    }


    {
      TGraphSmooth graphSmoothManager;
      calibTask.offsetTopGraphSmooth_ =
          new TGraphErrors(*dynamic_cast<TGraphErrors *>(graphSmoothManager.SmoothLowess(calibTask.offsetTopGraph_,
                                                                                         "",
                                                                                         offsetSpan, 5)));
      outDir->WriteObject(calibTask.offsetTopGraph_, "grOffsetTop");
      outDir->WriteObject(calibTask.offsetTopGraphSmooth_, "grOffsetTopLowess");
    }


    calibTask.isDone_ = true;
  }

  outputFile->Close();

  Info(__func__, "File %s is written succesfully", outputFile->GetName());
  Info(__func__, "Analyzer finished successfully.");
  return 0;
}
