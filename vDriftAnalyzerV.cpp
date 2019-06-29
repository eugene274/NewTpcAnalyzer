#include "vDriftAnalyzerV.h"
#include <algorithm>
#include <iomanip>
#include <TGraphErrors.h>
#include <TChain.h>
#include <TGraphSmooth.h>

int main(int argc, char **argv) {
  std::cout << "Starting analyser." << std::endl;

  std::string iFileName;
  TFile *inputFile;
  TFile *f_out;
  bool is_verbose;
  bool is_root_output;
  bool is_no_multiply;
  bool is_no_TOF_factor;
  // long int prev_time;
  std::string s_cut;
  long int prev_time = 0;
  is_verbose = false;
  is_root_output = false;
  is_no_multiply = false;
  is_no_TOF_factor = false;

  //-------------------------------USAGE-----------------------------------------//
  if (argc < 3) {
    Usage(argv);
    return error_code[0];
  }

  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) != "-i" &&
        std::string(argv[i]) != "-v" &&
        std::string(argv[i]) != "--root-output" &&
        std::string(argv[i]) != "--no-factor-multiplication" &&
        std::string(argv[i]) != "--no-TOF-MTPC-factor" &&
        std::string(argv[i]) != "--check-algorithm") {
      std::cerr << "\n[ERROR]: Unknown parameter " << i << ": " << argv[i] << std::endl;
      Usage(argv);
      std::cout << "Analyzer finished with exit code " << error_code[1] << std::endl;
      return error_code[1];
    } else {
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
      if (std::string(argv[i]) == "-v") {
        is_verbose = true;
        continue;
      }
      if (std::string(argv[i]) == "--root-output") {
        is_root_output = true;
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
      if (std::string(argv[i]) == "--check-algorithm") {
        is_no_multiply = true;
        is_no_TOF_factor = true;
        continue;
      }
    }
  }
  //-------------------------------USAGE-----------(end)-------------------------//

  std::cout << "\nStarting options:" << std::endl;
  std::cout << "\nFile        : " << iFileName.c_str() << std::endl;
  if (is_verbose) std::cout << "Verbose mode          : ON" << std::endl;
  if (!is_verbose) std::cout << "Verbose mode          : OFF" << std::endl;
  if (is_root_output) std::cout << "ROOT output           : ON" << std::endl;
  if (!is_root_output) std::cout << "ROOT output           : OFF" << std::endl;
  if (is_no_multiply) std::cout << "Multiplication factor : OFF - !!! Make sure this is intentional !!!" << std::endl;
  if (is_no_TOF_factor) std::cout << "MTPC-TOF factor       : OFF - !!! Make sure this is intentional !!!" << std::endl;
  std::cout << std::endl;

  if (check_file_exists(iFileName.c_str())) {
    inputFile = new TFile(iFileName.c_str(), "ro");
  } else {
    std::cerr << "\n[ERROR]: No file " << iFileName.c_str() << " was found!" << std::endl;
    std::cout << "Analyzer finished with exit code " << error_code[3] << std::endl;
    return error_code[3];
  }
  //-----------------------------------------------------------------------------//

  const double Ylimit = 2.;
  const int Nbins = 20;
  std::map<std::string, double> Xlimit;
  Xlimit["MTPCLvsTOFL"] = 55;
  Xlimit["MTPCRvsTOFR"] = 55;
  Xlimit["VTPC2vsMTPCL"] = 40;
  Xlimit["VTPC2vsMTPCR"] = 40;
  Xlimit["VTPC1vsVTPC2"] = 20;
  Xlimit["GTPCvsVTPC2"] = 15; //added by Scott; not sure what the best value should be      

  std::cout << "Getting TTree* from file." << std::endl;

  std::map<std::string, TTree *> dataTrees;

  GetCalibTreesFromFile(inputFile, dataTrees);

  std::map<std::string, vDriftTreeStructure> inData;

  if (is_root_output) f_out = new TFile((iFileName + "_output.root").c_str(), "recreate");

  std::vector<std::string> v_order;
  std::vector<double> v_slope;
  std::vector<double> v_slope_prev;
  std::map<std::string, TGraph *> m_vDslope;
  std::map<std::string, TGraph *> m_vDtime;
  std::map<std::string, TGraph *> m_vDnom;
  std::map<std::string, std::string> filename;
  std::map<std::string, int> TofFactor;

  double v_drift = 0.;
  double v_drift_average;
  long int time_average = 0;
  double v_drift_prev;

  //-----------------------TPC-order-of-calibration------------------------------//
  v_order.push_back("MTPCLvsTOFL");
  v_order.push_back("VTPC2vsMTPCL");
  v_order.push_back("VTPC1vsVTPC2");
  v_order.push_back("MTPCRvsTOFR");
  v_order.push_back("VTPC2vsMTPCR");
//  v_order.push_back("GTPCvsVTPC2");

  std::vector<float> v_tpc_ymax = {55, 55, 35, 55, 35, 15};

  //-----------------------TPC-order-of-calibration-(end)------------------------//

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

    double smoothingSpan_{600};

    double yBottom_{0};

    double xLo{-999};
    double xHi{999};

    /*** slice data ***/
    vDriftTreeStructure trackMatchData_{};

    TH2D *hdYvsY{nullptr};
    TProfile *pdYvsY{nullptr};

    TTree *sourceTree_{nullptr};

    /*** results ***/
    TGraph *recVDriftGraph{nullptr};
    TGraphErrors *calibVDriftGraph_{nullptr};

    TGraphErrors *offsetGraph_{nullptr};
    TGraphErrors *offsetBottomGraph_{nullptr};
    TGraphErrors *offsetBottomGraphSmooth_{nullptr};

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

    CalibTask_t &smoSpan(double span) {
      this->smoothingSpan_ = span;
      return *this;
    }

    /*** initialization ***/
    void initHistograms() {
      delete hdYvsY;
      delete pdYvsY;
      delete slopeGraph_;
      delete offsetGraph_;

      float yLo = gSTATIC_INFO.at(tgt).yLimits.first;
      float yHi = gSTATIC_INFO.at(tgt).yLimits.second;
      hdYvsY = new TH2D("hdYvsY", ";Y (cm); dY (cm)", ny_, yLo, yHi, ndy_, -2., 2.);
      hdYvsY->SetDirectory(nullptr);
      pdYvsY = new TProfile("pdYvsY", ";Y (cm); dY (cm)", ny_, yLo, yHi);
      pdYvsY->SetDirectory(nullptr);

      recVDriftGraph = new TGraph;

      offsetGraph_ = new TGraphErrors;
      offsetBottomGraph_ = new TGraphErrors;

      slopeGraph_ = new TGraphErrors;
    }

    void initInput(TFile &fIn) {
      delete sourceTree_;
      fIn.GetObject(gSTATIC_INFO.at(tgt).treeName, sourceTree_);
      ReadBranchesFromTree(sourceTree_, trackMatchData_, swap_ ? "swap" : "");
    }

    /*** getters ***/
    bool hasUpstream() const { return upstream < kDetPairs; };
  };

  bool usePropagation = true;

  std::vector<CalibTask_t> calibTasks{
      CalibTask_t("MTPCL", kMTPCLvsTOFL).ne(20000).smoSpan(1800).bottom(-60),
      CalibTask_t("VTPC2", kVTPC2vsMTPCL, kMTPCLvsTOFL).ne(2000).ndy(40).xLim(-110, 60).bottom(-55),
      CalibTask_t("VTPC1", kVTPC1vsVTPC2, kVTPC2vsMTPCL).ne(2000).ndy(40).bottom(-35),
      CalibTask_t("MTPCR", kVTPC2vsMTPCR, kVTPC2vsMTPCL).ne(2000).swap(),
      CalibTask_t("MTPCRfromTOFR", kMTPCRvsTOFR).ne(20000).smoSpan(1800),
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
    calibTask.initInput(*inputFile);
    auto &trackMatchData = calibTask.trackMatchData_;

    bool hasUpstream = calibTask.hasUpstream();
    CalibTask_t *upstreamTask = nullptr;
    for (auto &tt : calibTasksOrdered) {
      if (tt.tgt == calibTask.upstream) {
        upstreamTask = &tt;
        break;
      }
    }

    TDirectory *outDir = f_out->mkdir(calibTask.tgtName.c_str(), "");

    Info(__func__, "Processing %s", gSTATIC_INFO.at(calibTask.tgt).treeName);


    long nEntries = calibTask.sourceTree_->GetEntries();
    long iEntry = 0;

    int iSlice = 0;
    int iSliceEntries = 0;

    long currentEventNumber = -1;
    long currentUnixTime = -1;

    long sliceStartUnixTime = -1;

    int iRecVD = 0;
    while (iEntry < nEntries) {
      calibTask.sourceTree_->GetEntry(iEntry);

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
        calibTask.recVDriftGraph->SetPoint(iRecVD, unixTime, 1000*recVDrift);
        ++iRecVD;
      }

      if (usePropagation && upstreamTask) {
        assert(upstreamTask->isDone_);
        double tofFactor = gSTATIC_INFO.at(upstreamTask->tgt).tofFactor;
        double upstreamSlope = tofFactor * upstreamTask->slopeGraphSmooth_->Eval(unixTime);
        double upstreamOffset = tofFactor * upstreamTask->offsetBottomGraphSmooth_->Eval(unixTime);

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

          double slope = fitFun->GetParameter(1);
          double slopeError = fitFun->GetParError(1);

          calibTask.offsetGraph_->SetPoint(iSlice, sliceT, offset);
          calibTask.offsetGraph_->SetPointError(iSlice, 0., offsetError);

          calibTask.offsetBottomGraph_->SetPoint(iSlice, sliceT, offsetBottom);
          calibTask.offsetBottomGraph_->SetPointError(iSlice, 0., offsetBottomError);

          calibTask.slopeGraph_->SetPoint(iSlice, sliceT, slope);
          calibTask.slopeGraph_->SetPointError(iSlice, 0., slopeError);

          outDir->WriteObject(calibTask.hdYvsY, Form("hdYvsY_%d", iSlice));
          outDir->WriteObject(calibTask.pdYvsY, Form("pdYvsY_%d", iSlice));

          Info(__func__, "Slice %d: %d entries", iSlice, iSliceEntries);
          Info(__func__, "dT = %ld sec", dT);
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
    double span = calibTask.smoothingSpan_ / (tEnd - tStart);

    {
      TGraphSmooth graphSmoothManager;
      calibTask.slopeGraphSmooth_ =
          new TGraphErrors(*dynamic_cast<TGraphErrors *>(graphSmoothManager.SmoothLowess(calibTask.slopeGraph_,
                                                                                         "",
                                                                                         span)));
      outDir->WriteObject(calibTask.slopeGraph_, "grSlope");
      outDir->WriteObject(calibTask.slopeGraphSmooth_, "grSlopeLowess");
      calibTask.slopeGraph_->SetBit(TGraph::kIsSortedX);
      calibTask.slopeGraphSmooth_->SetBit(TGraph::kIsSortedX);

      calibTask.calibVDriftGraph_ = new TGraphErrors(calibTask.recVDriftGraph->GetN());
      for (int ip = 0; ip < calibTask.calibVDriftGraph_->GetN(); ++ip) {
        double unixTime{-1};
        double recVDrift{-1};
        calibTask.recVDriftGraph->GetPoint(ip, unixTime, recVDrift);
        double slope = calibTask.slopeGraphSmooth_->Eval(unixTime);
        double tofFactor = gSTATIC_INFO.at(calibTask.tgt).tofFactor;

        double calibVDrift = recVDrift*1.0/(1 + tofFactor*slope);
        calibTask.calibVDriftGraph_->SetPoint(ip, unixTime, calibVDrift);
      }
      outDir->WriteObject(calibTask.calibVDriftGraph_, "grCalibVDrift");

    }

    {
      TGraphSmooth graphSmoothManager;
      calibTask.offsetBottomGraphSmooth_ =
          new TGraphErrors(*dynamic_cast<TGraphErrors *>(graphSmoothManager.SmoothLowess(calibTask.offsetBottomGraph_,
                                                                                         "",
                                                                                         span)));
      outDir->WriteObject(calibTask.offsetGraph_, "grOffset");
      outDir->WriteObject(calibTask.offsetBottomGraph_, "grOffsetBottom");
      outDir->WriteObject(calibTask.offsetBottomGraphSmooth_, "grOffsetBottomLowess");
      outDir->WriteObject(calibTask.recVDriftGraph, "grRecVDrift");

      calibTask.offsetGraph_->SetBit(TGraph::kIsSortedX);
      calibTask.offsetBottomGraph_->SetBit(TGraph::kIsSortedX);
      calibTask.offsetBottomGraphSmooth_->SetBit(TGraph::kIsSortedX);
    }
    calibTask.isDone_ = true;
  }


  /*** =================================================================== ***/

  filename["MTPCLvsTOFL"] = "MTPCL";
  filename["MTPCRvsTOFR"] = "MTPCRfromTOFR";
  filename["VTPC2vsMTPCL"] = "VTPC2";
  filename["VTPC2vsMTPCR"] = "MTPCR";
  filename["VTPC1vsVTPC2"] = "VTPC1";
  filename["GTPCvsVTPC2"] = "GTPC";

  if (!is_no_TOF_factor) TofFactor["MTPCLvsTOFL"] = 2;
  if (is_no_TOF_factor) TofFactor["MTPCLvsTOFL"] = 1;
  if (!is_no_TOF_factor) TofFactor["MTPCRvsTOFR"] = 2;
  if (is_no_TOF_factor) TofFactor["MTPCRvsTOFR"] = 1;
  TofFactor["VTPC2vsMTPCL"] = 1;
  TofFactor["VTPC2vsMTPCR"] = 1;
  TofFactor["VTPC1vsVTPC2"] = 1;
  TofFactor["GTPCvsVTPC2"] = 1;

  if (is_verbose) {
    std::cout << "\nCurrent TPC calibration order:" << std::endl;
    for (unsigned int i = 0; i < v_order.size(); i++) {
      std::cout << v_order.at(i) << std::endl;
    }
    std::cout << std::endl;
  }

  // // // // // // // // // // // // // //

  std::map<std::string, int> nEntriesPerBin;

  //Scott: for low multiplicity, I used 1000 for drift velocity (stage2) and 10000 for y0 (stage4)
  nEntriesPerBin["MTPCLvsTOFL"] = 10000;
  nEntriesPerBin["MTPCRvsTOFR"] = 10000;
  nEntriesPerBin["VTPC2vsMTPCL"] = 10000;
  nEntriesPerBin["VTPC2vsMTPCR"] = 10000;
  nEntriesPerBin["VTPC1vsVTPC2"] = 10000;
  nEntriesPerBin["GTPCvsVTPC2"] = 10000;

  std::vector<TProfile2D *> v_dY_Y;
  std::vector<TProfile *> v_vD;
  std::vector<TH1D *> v_vD_corr;
  std::vector<TH1D *> v_hist_slope;
  std::vector<TH1D *> v_hist_offset;

  double max_time = dataTrees["VTPC1vsVTPC2"]->GetMaximum("eventUnixTime");
  double min_time = dataTrees["VTPC1vsVTPC2"]->GetMinimum("eventUnixTime");

  // // // // // // // // // // // // // // 


  for (unsigned int i = 0; i < v_order.size(); i++) {

    const long int nEntries = dataTrees[v_order.at(i)]->GetEntriesFast();

    if (is_verbose) std::cout << "\n" << v_order.at(i).c_str() << ": " << nEntries << " entries" << std::endl;
    if (!is_verbose) std::cout << "\n" << (v_order.at(i)).c_str() << std::endl;
    if (is_root_output) {
      f_out->cd();
      f_out->mkdir(v_order.at(i).c_str());
      f_out->cd(v_order.at(i).c_str());
    }
    if (nEntries < 2000) {
      std::cerr << "\n[WARNING] No data in this tree. Skipping..." << std::endl;
      continue;
    }

    const int nBins = nEntries / nEntriesPerBin[v_order.at(i)] + 1;

    TString name = "dY_" + v_order.at(i);
    v_dY_Y.push_back(new TProfile2D(name,
                                    name,
                                    nBins,
                                    min_time,
                                    max_time,
                                    Nbins,
                                    -Xlimit[v_order.at(i)] - v_tpc_ymax.at(i),
                                    Xlimit[v_order.at(i)] - v_tpc_ymax.at(i),
                                    -Ylimit,
                                    Ylimit));

    name = "vD_" + v_order.at(i);
    v_vD.push_back(new TProfile(name, name, nBins, min_time, max_time));

    std::cout << "Collecting unbinned data." << std::endl;

    if (v_order.at(i) != "VTPC2vsMTPCR") ReadBranchesFromTree(dataTrees[v_order.at(i)], inData[v_order.at(i)]);
    if (v_order.at(i) == "VTPC2vsMTPCR") ReadBranchesFromTree(dataTrees[v_order.at(i)], inData[v_order.at(i)], "swap");

    std::map<int, std::vector<double> > vY;
    std::map<int, std::vector<double> > vDY;

    for (long int i_tree = 0; i_tree < nEntries; i_tree++) {
      dataTrees[v_order.at(i)]->GetEntry(i_tree);

      const float
          y = inData[v_order.at(i)].master_Y - v_tpc_ymax.at(i);  // shift to have 0 at ~ y=y_max for a given TPC
      float dy = inData[v_order.at(i)].slave_Y - inData[v_order.at(i)].master_Y;
      const long int time = inData[v_order.at(i)].eventUnixTime;

      if (v_order.at(i) == "VTPC1vsVTPC2" || v_order.at(i) == "VTPC2vsMTPCR" || v_order.at(i) == "VTPC2vsMTPCL"
          || v_order.at(i) == "GTPCvsVTPC2") {
        int iPrevTPC = i - 1;
        if (v_order.at(i) == "VTPC2vsMTPCR") {
          auto it = find(v_order.begin(), v_order.end(), "VTPC2vsMTPCL");
          if (it == v_order.end()) {
            std::cout << "Error" << std::endl;
            exit;
          }
          iPrevTPC = it - v_order.begin();
          //                     std::cout << "iPrevTPC = " << iPrevTPC << std::endl;
        } else if (v_order.at(i) == "GTPCvsVTPC2") {
          auto it = find(v_order.begin(), v_order.end(), "VTPC2vsMTPCL");
          if (it == v_order.end()) {
            std::cout << "Error" << std::endl;
            exit;
          }
          iPrevTPC = it - v_order.begin();
        }
        const int bin = v_hist_slope.at(iPrevTPC)->GetXaxis()->FindBin(time);
        if (bin > 0) {
          const float prev_slope = TofFactor[v_order.at(iPrevTPC)] * v_hist_slope.at(iPrevTPC)->GetBinContent(bin);
          const float prev_offset = TofFactor[v_order.at(iPrevTPC)] * v_hist_offset.at(iPrevTPC)->GetBinContent(bin);
          dy += y * prev_slope + prev_offset;
        }
      }
      if (TMath::Abs(dy) <= Ylimit) {
        v_dY_Y.back()->Fill(time, y, dy);
        v_vD.back()->Fill(time, inData[v_order.at(i)].slave_recVDrift);
      }
    }

    name = "slope_" + v_order.at(i);
    v_hist_slope.push_back(new TH1D(name, name, nBins, min_time, max_time));
    name = "offset_" + v_order.at(i);
    v_hist_offset.push_back(new TH1D(name, name, nBins, min_time, max_time));

    name = "vD_calib_" + v_order.at(i);
    v_vD_corr.push_back(new TH1D(name, name, nBins, min_time, max_time));

    for (int iBin = 0; iBin < v_dY_Y.back()->GetXaxis()->GetNbins(); ++iBin) {
      TH1D *h = v_dY_Y.back()->ProjectionY(Form("y_%d", iBin), iBin + 1, iBin + 1);

      //if (h->GetEntries() == 0) continue;  //Scott: this did not work for me - the variable is not set correctly
      //h->Fit("pol1", "Q");
      if (h->Fit("pol1", "Q") != 0) continue; //Scott: this prevents bad fits from crashing the routine
      const float slope = h->GetFunction("pol1")->GetParameter(1);
      const float slope_err = h->GetFunction("pol1")->GetParError(1);

      const float offset = h->GetFunction("pol1")->GetParameter(0);
      const float offset_err = h->GetFunction("pol1")->GetParError(0);

      std::cout << slope << " +- " << slope_err << std::endl;

      if (is_root_output)
        h->Write();

      v_hist_slope.back()->SetBinContent(iBin + 1, slope);
      v_hist_slope.back()->SetBinError(iBin + 1, slope_err);

      v_hist_offset.back()->SetBinContent(iBin + 1, offset);
      v_hist_offset.back()->SetBinError(iBin + 1, offset_err);

      const float vD_dcs = v_vD.back()->GetBinContent(iBin + 1);
      const float vD_calib = vD_dcs / (1 + TofFactor[v_order.at(i)] * slope);
      const float vD_calib_err = vD_dcs * slope_err;

      v_vD_corr.back()->SetBinContent(iBin + 1, vD_calib);
      v_vD_corr.back()->SetBinError(iBin + 1, vD_calib_err);
    }
    if (is_root_output) {
      v_hist_slope.back()->Write();
      v_hist_offset.back()->Write();
      v_dY_Y.back()->Write();
      v_vD.back()->Write();
      v_vD_corr.back()->Write();
    }
    //         std::ofstream myfile;
    //         std::cout << "Writing txt output file: " << filename[v_order.at(i)] << ".txt" << std::endl;
    // 		myfile.open ((filename[v_order.at(i)]+".txt").c_str());
    //         
    //         for (int iTime=0; iTime<v_vD_corr.back()->GetNbinsX(); ++iTime)
    //         {
    //             long int time = v_vD_corr.back()->GetBinCenter(iTime+1);
    //             if (iTime == 0) time -= v_vD_corr.back()->GetBinWidth(0)/2;
    //             if (iTime == v_vD_corr.back()->GetNbinsX() - 1 ) time += v_vD_corr.back()->GetBinWidth(0)/2;
    //             
    //             myfile << std::setprecision(10)  << time << " " << 1e3 * v_vD_corr.back()->GetBinContent(iTime+1) << std::endl;
    //         }
    //         myfile.close();
  }

  if (is_root_output) f_out->Close();

  std::cout << "\nAnalyzer finished successfully." << std::endl;
  return 0;
}
