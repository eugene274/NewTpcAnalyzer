// // root -l -b -q 'plot_res.C("/eos/user/v/viklochk/public/vD/26455_stage2/20180611_112751/merged_26455.root", "res_26455.png")'
// // root -l -b -q 'plot_res.C("/eos/user/v/viklochk/public/vD/26603_stage2/20180611_112916/merged_26603.root", "res_26603.png")'
// // root -l -b -q 'plot_res.C("/eos/user/v/viklochk/public/vD/26803_stage2/20180611_113517/trackMatchDump_run-026803x003.root", "res_26803.png")'
// // root -l -b -q 'plot_res.C("/eos/user/v/viklochk/public/vD/26860_stage2/20180611_113131/merged_26860.root", "res_26860.png")'

const int nTpcs = 5;
const float fit_range[nTpcs][2] = { {-20,20}, {-30,30}, {-30,30}, {-50,50}, {-50,50} };

void plot_axis(const float xmin, const float xmax, const float ymin, const float ymax, 
               const TString title="", const TString xtitle="", const TString ytitle="")
{
    TH1F *h = new TH1F("h" + title, "", 1, xmin, xmax);
    h->SetTitle(title);

    h->GetYaxis()->SetTitle(ytitle);
    h->GetXaxis()->SetTitle(xtitle);

    h->GetYaxis()->SetRangeUser(ymin, ymax);
    h->GetXaxis()->SetRangeUser(xmin, xmax);
    
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleSize(0.06);
    
    h->GetXaxis()->SetTitleOffset(0.78);
    h->GetYaxis()->SetTitleOffset(0.98);
    
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetLabelSize(0.05);
    
    h->GetXaxis()->SetNdivisions(505);
    h->GetYaxis()->SetNdivisions(505);
    
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.12);
    
    h->Draw("");
}

void plot_res(const TString filename, const TString outfilename)
{
    TFile *f = TFile::Open(filename, "read");
    gStyle->SetOptStat(00000);
    
    const TString tpcs[nTpcs] = {"VTPC1vsVTPC2", "VTPC2vsMTPCL", "VTPC2vsMTPCR", "MTPCLvsTOFL", "MTPCRvsTOFR"};
    
    TCanvas *c = new TCanvas("c", "", 1200, 800);
    c->Divide(3,2);
    
    for (int i=0; i<nTpcs; ++i)
    {
        c->cd(i+1);
        TTree *t = (TTree*) f->Get(tpcs[i]);
        if (tpcs[i] != "MTPCLvsTOFL" && tpcs[i] != "MTPCRvsTOFR")
        {
            t->Draw( Form("slave_Y-master_Y : master_Y >> h%d(100,-50,50,100,-1,1)", i), "", "colz");
            t->Draw( Form("slave_Y-master_Y : master_Y >> p%d(20,-50,50)", i), "abs(slave_Y-master_Y)<1.", "profile");
        }
        else
        {
            t->Draw( Form("slave_Y-tofY : tofY >> h%d(100,-50,50,100,-3,3)", i), "", "colz");
            t->Draw( Form("slave_Y-tofY : tofY >> p%d(20,-50,50,-3,3)", i), "", "profile");
        }
        
        TH2D *h = gDirectory->Get( Form("h%d", i) );
        TProfile *p = gDirectory->Get( Form("p%d", i) );
        
        if (tpcs[i] != "MTPCLvsTOFL" && tpcs[i] != "MTPCRvsTOFR")
            h->SetMinimum(10);
        
        if (tpcs[i] != "MTPCLvsTOFL" && tpcs[i] != "MTPCRvsTOFR")
            plot_axis(-50, 50, -1, 1, tpcs[i], "y (cm)", "#Delta y (cm)");
        else
            plot_axis(-50, 50, -3, 3, tpcs[i], "y (cm)", "#Delta y (cm)");
            
        h->Draw("colz, same");
        p->Draw("same");
        
        p->Fit("pol1", "", "", fit_range[i][0], fit_range[i][1] );
        
        const float slope = p->GetFunction("pol1")->GetParameter(1);
        const float offset = p->GetFunction("pol1")->GetParameter(0);
        
        std::cout << slope << " " << offset << std::endl;
        
        TLatex latex;
        latex.SetTextSize(0.06);
        latex.SetTextColor(kBlack);  //align at top
        latex.SetTextAlign(13);  //align at top
        latex.DrawLatexNDC(.2,.88, Form("slope = %f", slope) );
        latex.DrawLatexNDC(.2,.8, Form("offset = %f", offset) );
    }
    
    c->SaveAs(outfilename);
}
