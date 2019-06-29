//
// Created by eugene on 6/29/19.
//

#include <iostream>
#include <vector>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TFile.h>
#include <fstream>

using namespace std;

struct CalibResults_t {
  std::string det{};

  TMultiGraph mgRecVDrift{};
  TMultiGraph mgCalibVDrift{};

  bool recreateASCII{true};

  CalibResults_t () = default;
  CalibResults_t (const CalibResults_t& other) : det(other.det) {};
  explicit CalibResults_t(const string &det) : det(det) {}
};

void exportASCII(TGraph *graph, const std::string& filename, bool recreate = false) {
  if (!graph) {
    return;
  }

  auto mode = recreate? ios_base::trunc : ios_base::app;
  ofstream outputASCII(filename, mode);

  for (int ip = 0; ip < graph->GetN(); ++ip) {
    double unixTime{-1};
    double vDrift{-1};
    graph->GetPoint(ip, unixTime, vDrift);
    outputASCII << long(unixTime) << " " << vDrift << endl;
  }
}

int main(int argc, char **argv) {


  vector<CalibResults_t> calibResults{
      CalibResults_t("MTPCL"),
      CalibResults_t("VTPC2"),
      CalibResults_t("VTPC1"),
      CalibResults_t("MTPCR"),
  };

  const vector<string> files{&argv[1], &argv[argc]};

  cout << "Files: " << endl;
  for (const auto &file : files) {
    cout << file << endl;
  }
  cout << endl;

  for (const auto &fileName : files) {
    TFile ff(fileName.c_str(), "read");
    if (!ff.IsOpen()) {
      Warning(__func__, "Skipping file %s", fileName.c_str());
      continue;
    }

    for (auto &cr : calibResults) {
      TGraph *recVDriftGraph = nullptr;
      ff.GetObject((cr.det + "/grRecVDrift").c_str(), recVDriftGraph);

      recVDriftGraph->SetLineColor(kRed);
      recVDriftGraph->SetMarkerColor(kRed);
      recVDriftGraph->SetLineWidth(2);
      cr.mgRecVDrift.Add(recVDriftGraph, "p");

      TGraph *calibVDriftGraph = nullptr;
      ff.GetObject((cr.det + "/grCalibVDrift").c_str(), calibVDriftGraph);

      calibVDriftGraph->SetLineColor(kBlue);
      calibVDriftGraph->SetMarkerColor(kBlue);

      calibVDriftGraph->SetLineWidth(2);
      cr.mgCalibVDrift.Add(calibVDriftGraph, "p");

      exportASCII(calibVDriftGraph, cr.det + ".txt", cr.recreateASCII);
      cr.recreateASCII = false;
    }
  }

  TFile outputFile("output.root", "recreate");
  for (auto &cr : calibResults) {
    outputFile.WriteObject(&cr.mgRecVDrift, Form("mgRecVDrift_%s", cr.det.c_str()));
    outputFile.WriteObject(&cr.mgCalibVDrift, Form("mgCalibVDrift_%s", cr.det.c_str()));
  }

  return 0;
}