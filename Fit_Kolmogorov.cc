#include <iostream>
#include <fstream>
#include <sstream>
#include <fstream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TVector.h"
#include "TVector3.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TFractionFitter.h>
#include <TLegend.h>
using namespace std;

const int eres_n_bin = 3;
const float eres_bin_min = 7;
const float eres_bin_max = 8;
const float eres_bin_width = (eres_bin_max-eres_bin_min)/(eres_n_bin-1);

const int gain_n_bin = 150;
const float gain_bin_min = 1/6e-05;
const float gain_bin_max = 1/3e-05;
const float gain_bin_width = (gain_bin_max-gain_bin_min)/(gain_n_bin-1);

const int charge_n_bin = 1024;
const float charge_bin_min = 0;
const float charge_bin_max = 200000;
const float charge_bin_width = (charge_bin_max-charge_bin_min)/(charge_n_bin-1);


TH1D* spectre_charge(int om_number )
  {
    TFile *file = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_energie_435.root", "READ");
    gROOT->cd();
    TH2F* charge = (TH2F*)file->Get("histo_pm_charge");

    TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
    file->Close();
    return spectre_charge;
    }


TH3D* MC_Simu(string name){

  TFile *file = new TFile(Form("Histo_simu/Simu_%s.root", name.c_str()), "READ");


  std::vector<int> *om_num = new std::vector<int>;
  std::vector<double> *energy = new std::vector<double>;

  TRandom3 rando;

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_id",1);
  tree->SetBranchAddress("om_id", &om_num);
  tree->SetBranchStatus("energy",1);
  tree->SetBranchAddress("energy", &energy);

  TH3D* MC_Simu = new TH3D(Form("MC_Simu_%s", name.c_str()), Form("MC_Simu_%s", name.c_str()),
                                eres_n_bin, eres_bin_min - (eres_bin_width/2), eres_bin_max + (eres_bin_width/2),
                                gain_n_bin, gain_bin_min - (gain_bin_width/2), gain_bin_max + (gain_bin_width/2),
                                charge_n_bin, charge_bin_min, charge_bin_max);


  for (int i = 0; i < 1000000; i++) {
  // for (int i = 0; i < tree->GetEntries(); i++) {
    double E_kolmo =0;
    tree->GetEntry(i);




    for (int j = 0; j < energy->size(); j++) {
      for (float eres = eres_bin_min; eres<= eres_bin_max; eres+= eres_bin_width) {
        float Evis = rando.Gaus(energy->at(j), (eres/235.482)*sqrt(energy->at(j)));
        for(float gain = gain_bin_min; gain<= gain_bin_max; gain+= gain_bin_width) {
            E_kolmo = gain * Evis;

            MC_Simu->Fill(eres, gain, E_kolmo);

        }
      }
    }
  }

  //
  // for (int bin_eres = 0; bin_eres < eres_n_bin; bin_eres++) {
  //   for(int bin_gain = 0; bin_gain< gain_n_bin; bin_gain++) {
  //     double integrale =
   // MC_Simu->Integral(bin_eres + 1, bin_eres + 1, bin_gain + 1, bin_gain + 1, charge_bin_min, charge_bin_max);
  //     // std::cout << integrale << '\n';
  //     for(int bin_charge = 0; bin_charge < charge_n_bin; bin_charge++){
  //       MC_Simu->SetBinContent(bin_eres+1, bin_gain+1, bin_charge +1, (MC_Simu->GetBinContent(bin_eres+1, bin_gain+1, bin_charge +1))/integrale);
  //       // MC_Simu->SetBinError(bin_eres+1, bin_gain+1, bin_charge +1, TMath::Sqrt((MC_Simu->GetBinContent(bin_eres+1, bin_gain+1, bin_charge +1)))/integrale);
  //     }
  //   }
  // }
 return MC_Simu;

}



void kolmo()
{
  gStyle->SetOptFit(1);
  TH1::SetDefaultSumw2();

  TH1D* spectre_om = spectre_charge(399);

        for (int bin =1; bin < 120; bin++) {
        spectre_om->SetBinContent(bin, 0);
      }
  TH3D* MC_Tl_208 = MC_Simu("Tl_208");
  TH3D* MC_Bi_214 = MC_Simu("Bi_214");
  TH3D* MC_K_40 = MC_Simu("K_40");

  double Chi2NDF = 0;
  double param1 = 0;
  double param2 = 0;
  double param3 = 0;
  double_t error1 =0;
  double_t error2 =0;
  double_t error3 = 0;

  float gain =0;
  float eres = 0;

  TFile *newfile = new TFile("histo_kolmo/Simu_kolmo.root", "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("Chi2NDF", &Chi2NDF);
  Result_tree.Branch("param1", &param1);
  Result_tree.Branch("param2", &param2);
  Result_tree.Branch("param3", &param3);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("eres", &eres);

  for (int eres_count = 3; eres_count < 4; eres_count++) {
    for (int gain_count = 1; gain_count <150; gain_count++) {
      TH1D *mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", eres_count, eres_count, gain_count, gain_count);    // first MC histogram
      TH1D *mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
      TH1D *mc2 = MC_K_40->ProjectionZ("Charge_K_40", eres_count, eres_count, gain_count, gain_count);    // second MC histogram

      for (int bin =1; bin < 120; bin++) {
        mc0->SetBinContent(bin, 0);
        mc1->SetBinContent(bin, 0);
        mc2->SetBinContent(bin, 0);
      }
                                              // retrieve histograms
     TObjArray *mc = new TObjArray(3);        // MC histograms are put in this array
     mc->Add(mc0);
     mc->Add(mc1);
     mc->Add(mc2);

     TFractionFitter* fit = new TFractionFitter(spectre_om, mc); // initialise
     fit->Constrain(0, 0, 1);               // constrain fraction 1 to be between 0 and 1
     fit->Constrain(1, 0, 1);
     fit->Constrain(2, 0, 1);


     // fit->SetRangeX(fitr,512);                // use only the first 15 bins in the fit
     Int_t status = fit->Fit();               // perform the fit

     std::cout << "fit status: " << status << std::endl;

     if (status == 0) {                       // check on fit status

       fit->GetResult(0, param1, error1);
       fit->GetResult(1, param2, error2);
       fit->GetResult(2, param3, error3);

       Chi2NDF = (fit->GetChisquare())/(fit->GetNDF());

       TH1D* result = (TH1D*) fit->GetPlot();
       TH1D* result_0 = (TH1D*) fit->GetMCPrediction(0);
       TH1D* result_1 = (TH1D*) fit->GetMCPrediction(1);
       TH1D* result_2 = (TH1D*) fit->GetMCPrediction(2);

       gain = 1/(gain_bin_min + gain_bin_width*(gain_count-1));
       eres = eres_bin_min + eres_bin_width*(eres_count-1);

       TCanvas* canvas = new TCanvas;
       canvas->SetLogy();
       spectre_om->Draw("same");
       spectre_om->SetTitle("Fit simu Bi+Tl+K");
       spectre_om->GetXaxis()->SetRangeUser(0, 120000);
       spectre_om->GetXaxis()->SetTitle("Charge (adc)");
       result_0->Draw("same");
       result_0->Scale(param1/result_0->Integral()*spectre_om->Integral());
       result_0->SetLineColor(kGreen);
       result_1->Draw("same");
       result_1->Scale(param2/result_1->Integral()*spectre_om->Integral());
       result_1->SetLineColor(kOrange);
       result_2->Draw("same");
       result_2->Scale(param3/result_2->Integral()*spectre_om->Integral());
       result_2->SetLineColor(kBlack);
       result->Draw("same");
       result->SetLineColor(kRed);
       auto legend = new TLegend(0.1,0.4,0.2,0.1);
       TLatex *t = new TLatex(.15,.15,Form("#Chi^{2}/NDF = %.3f",Chi2NDF));
       t->SetTextSize(0.04);
       t->Draw("same");
       legend->AddEntry(spectre_om, "data");
       legend->AddEntry(result_0, "Tl_208");
       legend->AddEntry(result_1, "Bi_214");
       legend->AddEntry(result_2, "K_40");
       legend->AddEntry(result, "fit");
       legend->Draw();
       canvas->SaveAs(Form("Fit_kolmo/triple/fit_kolmo_eres_%d_gain_%d.png", eres_count, gain_count));
       // canvas->SaveAs(Form("Fit_kolmo/double/fit_kolmo_eres_%d_gain_%d.png", eres_count, gain_count));




       Result_tree.Fill();
       delete fit;
      }
     }
 }


 newfile->cd();

 Result_tree.Write();

 newfile->Close();

}


int main(int argc, char const *argv[])
{
  kolmo();
  return 0;
}
