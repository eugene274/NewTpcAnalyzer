#include <fstream>
#include <vector>
#include <iomanip> 
#include <iostream>

#include <TString.h>
#include <TFile.h>
#include <TH1.h>
#include <TProfile.h>
#include <TGraphErrors.h>

void merge_output(const TString runlist="runs.txt", const TString indir="/eos/user/v/viklochk/public/vD//no_commit/20180605_163002/")
{
  const int nChambers = 6;

  std::ifstream infile(runlist);

  TFile *fIn {nullptr};

  std::vector <TString> tpcs = {"MTPCLvsTOFL", "MTPCRvsTOFR",	"VTPC2vsMTPCL", "VTPC2vsMTPCR",	"VTPC1vsVTPC2", "GTPCvsVTPC2"};
  std::vector <TString> tpc = {"MTPCL", "MTPCRvsTOFR",	"VTPC2", "MTPCR",	"VTPC1", "GTPC"};

  std::vector <double> v_vD_raw[nChambers];
  std::vector <double> v_vD_calib[nChambers];
  std::vector <double> v_slope[nChambers];
  std::vector <double> v_offset[nChambers];

  std::vector <double> v_vD_calib_err[nChambers];
  std::vector <double> v_slope_err[nChambers];
  std::vector <double> v_offset_err[nChambers];

  std::vector <double> v_time[nChambers];

  int run {-1};
  while (infile >> run)
  {
    std::cout << "Run # " << run << std::endl;

    fIn = TFile::Open( Form("%s/run_0%d_calib.root", indir.Data(), run), "read");
    if ( !fIn || fIn == nullptr )
      continue;


    for (int i=0; i<nChambers; ++i)
    {
      TH1D *h_offset = (TH1D*) fIn->Get( tpcs.at(i) + "/offset_" + tpcs.at(i) );

      //             std::cout << tpcs.at(i) << std::endl;

      if ( !h_offset || h_offset == nullptr )
        continue;

      TH1D *h_slope = (TH1D*) fIn->Get( tpcs.at(i) + "/slope_" + tpcs.at(i) );
      TH1D *h_vD_calib = (TH1D*) fIn->Get( tpcs.at(i) + "/vD_calib_" + tpcs.at(i) );
      TProfile *h_vD_raw = (TProfile*) fIn->Get( tpcs.at(i) + "/vD_" + tpcs.at(i) );

      for (int iTime=0; iTime<h_vD_calib->GetNbinsX(); ++iTime)
      {
        long int time = h_vD_calib->GetBinCenter(iTime+1);
        if (iTime == 0) time -= h_vD_calib->GetBinWidth(0)/2;
        if (iTime == h_vD_calib->GetNbinsX() - 1 ) time += h_vD_calib->GetBinWidth(0)/2;

        //                 std::cout << time << std::endl;

        v_time[i].push_back(time);
        v_vD_calib[i].push_back(h_vD_calib->GetBinContent(iTime+1));
        v_vD_raw[i].push_back(h_vD_raw->GetBinContent(iTime+1));
        v_slope[i].push_back(h_slope->GetBinContent(iTime+1));
        v_offset[i].push_back(h_offset->GetBinContent(iTime+1));

        v_vD_calib_err[i].push_back(h_vD_calib->GetBinError(iTime+1));
        v_slope_err[i].push_back(h_slope->GetBinError(iTime+1));
        v_offset_err[i].push_back(h_offset->GetBinError(iTime+1));

        //                 std::cout << std::setprecision(10)  << time << " " << 1e3 * h_vD_calib->GetBinContent(iTime+1) << std::endl;
      }
    }

  }

  TFile *fOut = TFile::Open("out.root", "recreate");

  for (int i=0; i<nChambers; ++i)
  {
    fOut->cd();
    fOut->mkdir(tpcs.at(i));
    fOut->cd(tpcs.at(i));

    TGraphErrors *g_vD_calib = new TGraphErrors( v_time[i].size(), &(v_time[i].at(0)), &(v_vD_calib[i].at(0)), 0, &(v_vD_calib_err[i].at(0)) );
    TGraphErrors *g_vD_raw = new TGraphErrors( v_time[i].size(), &(v_time[i].at(0)), &(v_vD_raw[i].at(0)), 0, 0 );
    TGraphErrors *g_slope = new TGraphErrors( v_time[i].size(), &(v_time[i].at(0)), &(v_slope[i].at(0)), 0, &(v_slope_err[i].at(0)) );
    TGraphErrors *g_offset = new TGraphErrors( v_time[i].size(), &(v_time[i].at(0)), &(v_offset[i].at(0)), 0, &(v_offset_err[i].at(0)) );

    g_vD_calib->Write("vD_calib");
    g_vD_raw->Write("vD_raw");
    g_slope->Write("vD_slope");
    g_offset->Write("vD_offset");

    std::ofstream myfile;
    std::cout << "Writing txt output file: " << tpc.at(i) << ".txt" << std::endl;
    myfile.open ( tpc.at(i)+".txt" );

    for (int iTime=0; iTime<g_vD_calib->GetN(); ++iTime)
    {
      const long int time = g_vD_calib->GetX()[iTime];
      const float vD = 1e3 * g_vD_calib->GetY()[iTime];
      if (vD > 0)
        myfile << std::setprecision(10)  << time << " " << vD << std::endl;
    }
    myfile.close();
  }

  fOut->Close();
}

int main(int argc, char** argv)
{
	if (argc<3){
		return -1;
	}
    
    merge_output( argv[1], argv[2] );
    return 0;
}
