#include "vDriftAnalyzerV.h"
#include <algorithm>
#include <iomanip>

int main(int argc, char** argv)
{
  std::cout << "Starting analyser." << std::endl;

  std::string iFileName;
  TFile* inputFile;
  TFile* f_out;
  bool is_verbose;
  bool is_root_output;
  bool is_no_multiply;
  bool is_no_TOF_factor;
  // long int prev_time;
  std::string s_cut;
  long int prev_time=0;
  is_verbose       = false;
  is_root_output   = false;
  is_no_multiply   = false;
  is_no_TOF_factor = false;

  //-------------------------------USAGE-----------------------------------------//
  if (argc<3){
    Usage(argv);
    return error_code[0];
  }

  for(int i=1; i<argc; i++){
    if(std::string(argv[i]) != "-i" &&
        std::string(argv[i]) != "-v" &&
        std::string(argv[i]) != "--root-output" &&
        std::string(argv[i]) != "--no-factor-multiplication" &&
        std::string(argv[i]) != "--no-TOF-MTPC-factor" &&
        std::string(argv[i]) != "--check-algorithm"){
      std::cerr << "\n[ERROR]: Unknown parameter " << i << ": " <<  argv[i] << std::endl;
      Usage(argv);
      std::cout << "Analyzer finished with exit code " << error_code[1] << std::endl;
      return error_code[1];
    } else {
      if(std::string(argv[i]) == "-i" && i!=argc-1) {
        iFileName = argv[++i];
        continue;
      }
      if(std::string(argv[i]) == "-i" && i==argc-1) {
        std::cerr << "\n[ERROR]: File name was not specified " << std::endl;
        Usage(argv);
        std::cout << "Analyzer finished with exit code " << error_code[2] << std::endl;
        return error_code[2];
      }
      if(std::string(argv[i]) == "-v") {
        is_verbose = true;
        continue;
      }
      if(std::string(argv[i]) == "--root-output") {
        is_root_output = true;
        continue;
      }
      if(std::string(argv[i]) == "--no-factor-multiplication") {
        is_no_multiply = true;
        continue;
      }
      if(std::string(argv[i]) == "--no-TOF-MTPC-factor") {
        is_no_TOF_factor = true;
        continue;
      }
      if(std::string(argv[i]) == "--check-algorithm") {
        is_no_multiply   = true;
        is_no_TOF_factor = true;
        continue;
      }
    }
  }
  //-------------------------------USAGE-----------(end)-------------------------//

  std::cout << "\nStarting options:" << std::endl;
  std::cout << "\nFile        : " << iFileName.c_str() << std::endl;
  if (is_verbose)  std::cout      << "Verbose mode          : ON" << std::endl;
  if (!is_verbose) std::cout      << "Verbose mode          : OFF" << std::endl;
  if (is_root_output) std::cout   << "ROOT output           : ON" << std::endl;
  if (!is_root_output) std::cout  << "ROOT output           : OFF" << std::endl;
  if (is_no_multiply) std::cout   << "Multiplication factor : OFF - !!! Make sure this is intentional !!!" << std::endl;
  if (is_no_TOF_factor) std::cout << "MTPC-TOF factor       : OFF - !!! Make sure this is intentional !!!" << std::endl;
  std::cout << std::endl;

  if (check_file_exists(iFileName.c_str())){
    inputFile = new TFile(iFileName.c_str(),"ro");
  } else {
    std::cerr << "\n[ERROR]: No file " << iFileName.c_str() << " was found!" << std::endl;
    std::cout << "Analyzer finished with exit code " << error_code[3] << std::endl;
    return error_code[3];
  }
  //-----------------------------------------------------------------------------//

  const double Ylimit = 2.;
  const int    Nbins  = 20;
  std::map<std::string,double> Xlimit;
  Xlimit["MTPCLvsTOFL"] = 55;
  Xlimit["MTPCRvsTOFR"] = 55;
  Xlimit["VTPC2vsMTPCL"] = 40;
  Xlimit["VTPC2vsMTPCR"] = 40;
  Xlimit["VTPC1vsVTPC2"] = 20;
  Xlimit["GTPCvsVTPC2"] = 15; //added by Scott; not sure what the best value should be      

  std::cout << "Getting TTree* from file." << std::endl;

  std::map<std::string, TTree*> dataTrees;

  GetCalibTreesFromFile(inputFile, dataTrees);

  std::map<std::string, vDriftTreeStructure> inData;

  if (is_root_output) f_out = new TFile((iFileName+"_output.root").c_str(),"recreate");

  std::vector<std::string>          v_order;
  std::vector<double>               v_slope;
  std::vector<double>               v_slope_prev;
  std::map<std::string,TGraph*>     m_vDslope;
  std::map<std::string,TGraph*>     m_vDtime;
  std::map<std::string,TGraph*>     m_vDnom;
  std::map<std::string,std::string> filename;
  std::map<std::string, int>        TofFactor;

  double v_drift      = 0.;
  double v_drift_average;
  long int time_average=0;
  double v_drift_prev;

  //-----------------------TPC-order-of-calibration------------------------------//
  v_order.push_back("MTPCLvsTOFL");
  v_order.push_back("VTPC2vsMTPCL");
  v_order.push_back("VTPC1vsVTPC2");
  v_order.push_back("MTPCRvsTOFR");
  v_order.push_back("VTPC2vsMTPCR");
//  v_order.push_back("GTPCvsVTPC2");

  std::vector <float> v_tpc_ymax = { 55, 55, 35, 55, 35, 15};

  //-----------------------TPC-order-of-calibration-(end)------------------------//

  filename["MTPCLvsTOFL"]   = "MTPCL";
  filename["MTPCRvsTOFR"]   = "MTPCRfromTOFR";
  filename["VTPC2vsMTPCL"]  = "VTPC2";
  filename["VTPC2vsMTPCR"]  = "MTPCR";
  filename["VTPC1vsVTPC2"]  = "VTPC1";
  filename["GTPCvsVTPC2"]  = "GTPC";

  if (!is_no_TOF_factor) TofFactor["MTPCLvsTOFL"] = 2;
  if (is_no_TOF_factor) TofFactor["MTPCLvsTOFL"] = 1;
  if (!is_no_TOF_factor) TofFactor["MTPCRvsTOFR"] = 2;
  if (is_no_TOF_factor) TofFactor["MTPCRvsTOFR"] = 1;
  TofFactor["VTPC2vsMTPCL"] = 1;
  TofFactor["VTPC2vsMTPCR"] = 1;
  TofFactor["VTPC1vsVTPC2"] = 1;
  TofFactor["GTPCvsVTPC2"] = 1;

  if (is_verbose){
    std::cout << "\nCurrent TPC calibration order:" << std::endl;
    for (unsigned int i=0; i<v_order.size();i++){
      std::cout << v_order.at(i) << std::endl;
    }
    std::cout << std::endl;
  }

  // // // // // // // // // // // // // //

  std::map<std::string, int> nEntriesPerBin;

  //Scott: for low multiplicity, I used 1000 for drift velocity (stage2) and 10000 for y0 (stage4)
  nEntriesPerBin["MTPCLvsTOFL"]  = 10000;
  nEntriesPerBin["MTPCRvsTOFR"]  = 10000;
  nEntriesPerBin["VTPC2vsMTPCL"] = 10000;
  nEntriesPerBin["VTPC2vsMTPCR"] = 10000;
  nEntriesPerBin["VTPC1vsVTPC2"] = 10000;
  nEntriesPerBin["GTPCvsVTPC2"] = 10000;

  std::vector <TProfile2D*> v_dY_Y;
  std::vector <TProfile*> v_vD;
  std::vector <TH1D*> v_vD_corr;
  std::vector <TH1D*> v_hist_slope;
  std::vector <TH1D*> v_hist_offset;

  double max_time = dataTrees["VTPC1vsVTPC2"]->GetMaximum("eventUnixTime");
  double min_time = dataTrees["VTPC1vsVTPC2"]->GetMinimum("eventUnixTime");

  // // // // // // // // // // // // // // 


  for (unsigned int i=0;i<v_order.size();i++){

    const long int nEntries = dataTrees[v_order.at(i)]->GetEntriesFast();

    if (is_verbose) std::cout << "\n" << v_order.at(i).c_str()<<  ": " << nEntries << " entries" << std::endl;
    if (!is_verbose) std::cout << "\n" << (v_order.at(i)).c_str() << std::endl;
    if (is_root_output)
    {
      f_out->cd();
      f_out->mkdir(v_order.at(i).c_str());
      f_out->cd(v_order.at(i).c_str());
    }
    if (nEntries < 2000){
      std::cerr << "\n[WARNING] No data in this tree. Skipping..." << std::endl;
      continue;
    }

    const int nBins = nEntries / nEntriesPerBin[v_order.at(i)] + 1;

    TString name = "dY_" + v_order.at(i);
    v_dY_Y.push_back( new TProfile2D( name, name, nBins, min_time, max_time, Nbins, -Xlimit[v_order.at(i)] - v_tpc_ymax.at(i), Xlimit[v_order.at(i)] - v_tpc_ymax.at(i), -Ylimit, Ylimit) );

    name = "vD_" + v_order.at(i);
    v_vD.push_back( new TProfile( name, name, nBins, min_time, max_time) );

    std::cout << "Collecting unbinned data." << std::endl;

    if (v_order.at(i) != "VTPC2vsMTPCR") ReadBranchesFromTree(dataTrees[v_order.at(i)], inData[v_order.at(i)]);
    if (v_order.at(i) == "VTPC2vsMTPCR") ReadBranchesFromTree(dataTrees[v_order.at(i)], inData[v_order.at(i)], "swap");

    std::map<int,std::vector<double> > vY;
    std::map<int,std::vector<double> > vDY;

    for (long int i_tree=0;i_tree<nEntries;i_tree++){
      dataTrees[v_order.at(i)]->GetEntry(i_tree);

      const float y = inData[v_order.at(i)].master_Y - v_tpc_ymax.at(i);  // shift to have 0 at ~ y=y_max for a given TPC
      float dy = inData[v_order.at(i)].slave_Y - inData[v_order.at(i)].master_Y;
      const long int time = inData[v_order.at(i)].eventUnixTime;

      if ( v_order.at(i) == "VTPC1vsVTPC2" || v_order.at(i) == "VTPC2vsMTPCR" || v_order.at(i) == "VTPC2vsMTPCL" || v_order.at(i) == "GTPCvsVTPC2")
      {
        int iPrevTPC = i-1;
        if (v_order.at(i) == "VTPC2vsMTPCR")
        {
          auto it = find (v_order.begin(), v_order.end(), "VTPC2vsMTPCL");
          if (it == v_order.end())
          {
            std::cout << "Error" << std::endl;
            exit;
          }
          iPrevTPC = it - v_order.begin();
          //                     std::cout << "iPrevTPC = " << iPrevTPC << std::endl;
        }
        else if (v_order.at(i) == "GTPCvsVTPC2")
        {
          auto it = find (v_order.begin(), v_order.end(), "VTPC2vsMTPCL");
          if (it == v_order.end())
          {
            std::cout << "Error" << std::endl;
            exit;
          }
          iPrevTPC = it - v_order.begin();
        }
        const int bin = v_hist_slope.at(iPrevTPC)->GetXaxis()->FindBin(time);
        if (bin > 0)
        {
          const float prev_slope = TofFactor[v_order.at(iPrevTPC)] * v_hist_slope.at(iPrevTPC)->GetBinContent(bin);
          const float prev_offset = TofFactor[v_order.at(iPrevTPC)] * v_hist_offset.at(iPrevTPC)->GetBinContent(bin);
          dy += y * prev_slope + prev_offset;
        }
      }
      if ( TMath::Abs(dy) <= Ylimit ){
        v_dY_Y.back()->Fill( time, y, dy );
        v_vD.back()->Fill( time, inData[v_order.at(i)].slave_recVDrift );
      }
    }

    name = "slope_" + v_order.at(i);
    v_hist_slope.push_back ( new TH1D( name, name, nBins, min_time, max_time ));
    name = "offset_" + v_order.at(i);
    v_hist_offset.push_back ( new TH1D( name, name, nBins, min_time, max_time ));

    name = "vD_calib_" + v_order.at(i);
    v_vD_corr.push_back ( new TH1D( name, name, nBins, min_time, max_time ));

    for (int iBin = 0; iBin < v_dY_Y.back()->GetXaxis()->GetNbins(); ++iBin)
    {
      TH1D *h = v_dY_Y.back()->ProjectionY( Form("y_%d", iBin) , iBin+1, iBin+1);

      //if (h->GetEntries() == 0) continue;  //Scott: this did not work for me - the variable is not set correctly
      //h->Fit("pol1", "Q");
      if(h->Fit("pol1", "Q")!=0) continue; //Scott: this prevents bad fits from crashing the routine
      const float slope = h->GetFunction("pol1")->GetParameter(1);
      const float slope_err = h->GetFunction("pol1")->GetParError(1);

      const float offset = h->GetFunction("pol1")->GetParameter(0);
      const float offset_err = h->GetFunction("pol1")->GetParError(0);

      std::cout << slope << " +- " << slope_err << std::endl;

      if (is_root_output)
        h->Write();

      v_hist_slope.back()->SetBinContent(iBin+1, slope);
      v_hist_slope.back()->SetBinError(iBin+1, slope_err);

      v_hist_offset.back()->SetBinContent(iBin+1, offset);
      v_hist_offset.back()->SetBinError(iBin+1, offset_err);

      const float vD_dcs = v_vD.back()->GetBinContent(iBin+1);
      const float vD_calib = vD_dcs / ( 1 + TofFactor[v_order.at(i)] * slope);
      const float vD_calib_err = vD_dcs * slope_err;

      v_vD_corr.back()->SetBinContent(iBin+1, vD_calib);
      v_vD_corr.back()->SetBinError(iBin+1, vD_calib_err);
    }
    if (is_root_output)
    {
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
