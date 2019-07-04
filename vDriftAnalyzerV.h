#ifndef VD_HEADER

#define VD_HEADER
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TCanvas.h"

#define _N_AGREGATE 300

//Patch namespace - used to fix known bugs & issues of the software
namespace patch {
template<typename T>
std::string to_string(const T &n)
//Bug in the old g++ compilator that crushes when using std::to_string function
{
  std::ostringstream stm;
  stm << n;
  return stm.str();
}
}

///////////////////////////////////////
const int error_code[] = {10, 11, 12, 13, 14};
//ERROR CODES:
////    Usage & input files:
////    10 - uncorrect usage of the file (see usage)
////    11 - uncorrect parameter (see usage)
////    12 - no file name after-i option (see usage)
////    13 - file not found
////    14 - uncorrect parameter in GetSlope("string",std::map<...>)


#pragma pack(push, 1)
struct vDriftTreeStructure {
  double slave_Y;
  double slave_X;
  double master_Y;
  double master_X;
  double slave_recVDrift;
  int sectionID;
  unsigned int runNumber;
  unsigned int eventUnixTime;
  unsigned int eventNumber;
};
#pragma pack(pop)

enum EDetectorPair : ushort { kMTPCLvsTOFL = 0, kMTPCRvsTOFR, kVTPC2vsMTPCL, kVTPC2vsMTPCR, kVTPC1vsVTPC2, kDetPairs };

struct DetPairStaticInfo {
  ushort id{kDetPairs};
  const char *treeName{};
  const char *masterYLeafName{};
  const char *slaveYLeafName{};
  float tofFactor{1.};
  std::pair<float,float> yLimits{};
};

std::array<DetPairStaticInfo, std::size_t(kDetPairs)> gSTATIC_INFO{
    {
        {kMTPCLvsTOFL, "MTPCLvsTOFL", "tofY", "slave_Y", 2.0, {-60, 60}},
        {kMTPCRvsTOFR, "MTPCRvsTOFR", "tofY", "slave_Y", 2.0, {-60, 60}},
        {kVTPC2vsMTPCL, "VTPC2vsMTPCL", "master_Y", "slave_Y", 1.0, {-40, 40}},
        {kVTPC2vsMTPCR, "VTPC2vsMTPCR", "master_Y", "slave_Y", 1.0, {-40, 40}},
        {kVTPC1vsVTPC2, "VTPC1vsVTPC2", "master_Y", "slave_Y", 1.0, {-20, 20}},
    }
};

void Usage(char **arg_string) {
  std::cerr << "\nUsage: " << arg_string[0] << " [options] [argument]" << std::endl;
  std::cerr << "  Options:" << std::endl;
  std::cerr << "    -i                         : set input file or filelist (with .list suffix) " << std::endl;
  std::cerr << "    -o                         : set output file " << std::endl;
//  std::cerr << "    --no-factor-multiplication : disable multiplication rule for calibration factor " << std::endl;
//  std::cerr << "    --no-TOF-MTPC-factor       : disable TOF-MTPC matching correction " << std::endl;
  std::cerr << "  Examples: " << std::endl;
  std::cerr << "    ./vDriftAnalyzer -i input_file.root -o test.root" << std::endl;
  std::cerr << "    ./vDriftAnalyzer -i \"*.root\"" << std::endl;
  std::cerr << "    ./vDriftAnalyzer -i \"filelist.list\"\n" << std::endl;
}

void ReadBranchesFromTree(TTree *tree, vDriftTreeStructure &data, std::string optSwap = "") {
  std::string s_name;
  s_name = (std::string) tree->GetName();
  if (optSwap == "swap") std::cout << "ReadBranchesFromTree(): Option 'swap' is ON." << std::endl;

  if (optSwap == "") {
    if (s_name.find("TOF") != std::string::npos) {
      tree->SetBranchAddress("tofY", &data.master_Y);
      tree->SetBranchAddress("tofX", &data.master_X);
    } else {
      tree->SetBranchAddress("master_Y", &data.master_Y);
      tree->SetBranchAddress("master_X", &data.master_X);
    }
    tree->SetBranchAddress("slave_Y", &data.slave_Y);
    tree->SetBranchAddress("slave_X", &data.slave_X);
    tree->SetBranchAddress("slave_recVDrift", &data.slave_recVDrift);
  }
  if (optSwap == "swap") {
    if (s_name.find("TOF") != std::string::npos) {
      tree->SetBranchAddress("tofY", &data.slave_Y);
      tree->SetBranchAddress("tofX", &data.slave_X);
      std::cerr << "[WARNING] ReadBranchesFromTree(): Attempt to read drift velocity from TOF!" << std::endl;
    } else {
      tree->SetBranchAddress("master_Y", &data.slave_Y);
      tree->SetBranchAddress("master_X", &data.slave_X);
      tree->SetBranchAddress("master_recVDrift", &data.slave_recVDrift);
    }
    tree->SetBranchAddress("slave_Y", &data.master_Y);
    tree->SetBranchAddress("slave_X", &data.master_X);
  }
  tree->SetBranchAddress("sectionId", &data.sectionID);
  tree->SetBranchAddress("runNumber", &data.runNumber);
  tree->SetBranchAddress("eventNumber", &data.eventNumber);
  tree->SetBranchAddress("eventUnixTime", &data.eventUnixTime);
}

#endif
