#define _USE_MATH_DEFINES
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

const int eres_n_bin = 53;
const float eres_bin_min = 7;
const float eres_bin_max = 20;
const float eres_bin_width = (eres_bin_max-eres_bin_min)/(eres_n_bin-1);

const int gain_n_bin = 221;
const float gain_bin_min = 10000;
const float gain_bin_max = 65000;
const float gain_bin_width = (gain_bin_max-gain_bin_min)/(gain_n_bin-1);

const int charge_n_bin = 1024;
const float charge_bin_min = 0e-05;
const float charge_bin_max = 200000;
const float charge_bin_width = (charge_bin_max-charge_bin_min)/(charge_n_bin-1);

TH1D* spectre_charge_it(int om_number )
{
  TFile *file = new TFile("histo_brut/histo_charge_amplitude_422.root", "READ");
  gROOT->cd();
  TH2F* charge = (TH2F*)file->Get("charge");

  TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  file->Close();
  return spectre_charge;
}

TH1D* spectre_charge_XW(int om_number )
{
  TFile *file = new TFile("histo_brut/histo_charge_amplitude_energie_449.root", "READ");
  gROOT->cd();
  TH2F* charge = (TH2F*)file->Get("histo_pm_charge");

  TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  file->Close();
  return spectre_charge;
}

TH1D* spectre_charge_fr(int om_number )
{
  TFile *file = new TFile("histo_brut/histo_charge_amplitude_energie_435.root", "READ");
  gROOT->cd();
  TH2F* charge = (TH2F*)file->Get("histo_pm_charge");

  TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
  file->Close();
  return spectre_charge;
}

TH1D* spectre_chooser(int om)
{  
  TH1D* spectre_chooser = NULL;
  if (om <260){
    spectre_chooser = spectre_charge_it(om);
  }
    else if (om <520) {
      spectre_chooser = spectre_charge_fr(om);
    }
    else if (om >519) {
      spectre_chooser = spectre_charge_XW(om);
    }
  return spectre_chooser;
}

double* om_gain_fit(int om)
{
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);

  double mean = 0;
  double sigma = 0;
  double n_evt = 0;
  double nbg = 0;
  double chi = 0;
  double ndf = 0;
  double chin = 0;

  double* tab = new double[7];
  TCanvas* canvas = new TCanvas;
  canvas->SetLogy();

  TH1D* spectre_om = NULL;
  if (om <260){
    spectre_om = spectre_charge_it(om);
  }
  else if (om <520) {
    spectre_om = spectre_charge_fr(om);
  }
  else if (om >519) {
    spectre_om = spectre_charge_XW(om);
  }

  TF1* f_ComptonEdgePoly = new TF1 ("f_ComptonEdgePoly","[0]/2.0*(1+TMath::Erf(([1]-x)/(TMath::Sqrt(2)*[2])))+[3]*x", 40000, 90000);
  f_ComptonEdgePoly->SetParNames("N_evt","Mean","Sigma","Nbg" );

  if ((om % 13) == 12 )        //om multiple de (13)-1
    {
      f_ComptonEdgePoly->SetParameters(120, 72367, 6893, 3.91e-5);
      f_ComptonEdgePoly->SetRange(6000,100000);
      f_ComptonEdgePoly->Draw("same");
      spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
      f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
      spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
      f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
      spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    }
  else if ((om % 13) == 0)       //om multiple de 13
    {
      f_ComptonEdgePoly->SetParameters(112, 68168, 5604, 1.2e-05);
      f_ComptonEdgePoly->SetRange(50000,100000);
      f_ComptonEdgePoly->Draw("same");
      spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
      f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
      spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
      f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-1.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
      spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    }
  else         //om normaux (8pouces)
    {
      f_ComptonEdgePoly->SetParameters(111, 60978, 3787, 4.19e-05);
      f_ComptonEdgePoly->SetRange(55000,100000);
      f_ComptonEdgePoly->Draw("same");
      spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
      f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+7.5*f_ComptonEdgePoly->GetParameter(2));
      spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
      f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-3.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+7.5*f_ComptonEdgePoly->GetParameter(2));
      spectre_om->Fit(f_ComptonEdgePoly, "RQ0");
    }
  chi = (f_ComptonEdgePoly->GetChisquare());
  std::cout << chi << '\n';
  ndf = (f_ComptonEdgePoly->GetNDF());
  chin = (f_ComptonEdgePoly->GetChisquare()/f_ComptonEdgePoly->GetNDF());

  if (chin < 1.5) {
    n_evt = f_ComptonEdgePoly->GetParameter(0);
    mean = (f_ComptonEdgePoly->GetParameter(1))/2.6;
    sigma = (f_ComptonEdgePoly->GetParameter(2))/2.6;
    nbg = f_ComptonEdgePoly->GetParameter(3);
  }
  else{
    n_evt = 0;
    mean = 0;
    sigma = 0;
    nbg = 0;
  }
  delete f_ComptonEdgePoly;
  delete canvas;
  tab[0] = mean;
  tab[1] = sigma;
  tab[2] = chin;
  tab[3] = chi;
  tab[4] = ndf;
  tab[5] = n_evt;
  tab[6] = nbg;
  return tab;
}

void fractionfitter(int om, int draw, double gain, double eres, TH1D* spectre_om, TH1D* mc0, TH1D* mc1, TH1D* mc2, double *fftab)
{

  double param1, param2, param3, error1, error2, error3 = 0;
  double result_0_scale, result_1_scale, result_2_scale = 0;
  mc0->Draw();
 

  cout << "ok" << endl;

  int lim = 105;
  for (int bin =1; bin < lim; bin++)
    {
      mc0->SetBinContent(bin, 0);
      mc1->SetBinContent(bin, 0);
      mc2->SetBinContent(bin, 0);
    }
  // retrieve histograms
  TObjArray *mc = new TObjArray(3);        // MC histograms are put in this array
  mc->Add(mc0);
  mc->Add(mc1);
  mc->Add(mc2);
  
  TFractionFitter* fit = new TFractionFitter(spectre_om, mc, "Q");          // initialise
  fit->Constrain(0, 0, 1);               // constrain fraction 1 to be between 0 and 1
  fit->Constrain(1, 0, 1);
  fit->Constrain(2, 0, 1);

  // fit->SetRangeX(fitr,512);                // use only the first 15 bins in the fit
  Int_t status = fit->Fit();               // perform the fit

  if (status == 0) {                       // check on fit status

    fit->GetResult(0, param1, error1);
    fit->GetResult(1, param2, error2);
    fit->GetResult(2, param3, error3);
    double Chi2NDF = (fit->GetChisquare())/(fit->GetNDF());
    fftab[0] = Chi2NDF;
    double param_tab [3];
    param_tab[0] = param1;
    param_tab[1] = param2;
    param_tab[2] = param3;

    fftab[1] = fit->EvaluateFCN(param_tab);
    cout << "          " << endl;
    cout << "FCN = " << fit->EvaluateFCN(param_tab) << endl;
    cout << "          " << endl;

    fftab[2] = param1;
    fftab[3] = param2;
    fftab[4] = param3;
    fftab[5] = error1;
    fftab[6] = error2;
    fftab[7] = error3;
    
    if (draw == 1 && fftab[0] < 0.4){

      TH1D* result = (TH1D*) fit->GetPlot();
      TH1D* result_0 = (TH1D*) fit->GetMCPrediction(0);
      TH1D* result_1 = (TH1D*) fit->GetMCPrediction(1);
      TH1D* result_2 = (TH1D*) fit->GetMCPrediction(2);

      TCanvas* canvas = new TCanvas;
      canvas->SetLogy();
      spectre_om->Draw("same");
      spectre_om->SetTitle(Form("Fit_simu Eres = %.2f, Gain = %.0f", eres, 1/gain));
      spectre_om->GetXaxis()->SetRangeUser(0, 120000);
      spectre_om->GetXaxis()->SetTitle("Charge (adc)");
      result_0->Draw("same");
      result_0_scale=(param1/result_0->Integral())*spectre_om->Integral();
      result_0->Scale(result_0_scale);
      result_0->SetLineColor(kGreen);
      result_1->Draw("same");
      result_1_scale=(param2/result_1->Integral())*spectre_om->Integral();
      result_1->Scale(result_1_scale);
      result_1->SetLineColor(kOrange);
      result_2->Draw("same");
      result_2_scale=(param3/result_2->Integral())*spectre_om->Integral();
      result_2->Scale(result_2_scale);
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
      
      canvas->SaveAs(Form("Best_fit/best_fit_om_%d_eres_%f_gain_%f.png", om, eres, 1/gain));

      delete canvas;	
    }
  } 
}


void kolmo(int draw)
{
  std::ifstream  charge("gain_data/Resultat_mean_energie_charge_total.txt");
  float charge_valeur_fit [712];
  memset(charge_valeur_fit, -1, 712*sizeof(float));
  int charge_om_num;
  while (charge >> charge_om_num)
  {
    charge >> charge_valeur_fit[charge_om_num];
  }

  TFile *histo_file_Tl = new TFile("Histo_simu_new/MC_simu_Tl_208_recOC_new_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_Tl->cd();
  TH3D* MC_Tl_208 = (TH3D*)histo_file_Tl->Get("MC_simu_Tl_208_recOC_new_ubc");

  TFile *histo_file_Bi = new TFile("Histo_simu_new/MC_simu_Bi_214_recOC_new_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_Bi->cd();
  TH3D* MC_Bi_214 = (TH3D*)histo_file_Bi->Get("MC_simu_Bi_214_recOC_new_ubc");

  TFile *histo_file_K = new TFile("Histo_simu_new/MC_simu_K_40_recOC_new_eres_53_gain_221_OC_MW8.root", "READ");
  histo_file_K->cd();
  TH3D* MC_K_40 = (TH3D*)histo_file_K->Get("MC_simu_K_40_recOC_new_ubc");


  // TFile *histo_file_Tl = new TFile("Histo_simu_new/Energy_MC_simu_Tl_208_recOC_new_eres_53_gain_221_wubc.root", "READ");
  // histo_file_Tl->cd();
  // TH3D* MC_Tl_208 = (TH3D*)histo_file_Tl->Get("MC_simu_Tl_208_recOC_new_ubc");

  // TFile *histo_file_Bi = new TFile("Histo_simu_new/Energy_MC_simu_Bi_214_recOC_new_eres_53_gain_221_wubc.root", "READ");
  // histo_file_Bi->cd();
  // TH3D* MC_Bi_214 = (TH3D*)histo_file_Bi->Get("MC_simu_Bi_214_recOC_new_ubc");

  // TFile *histo_file_K = new TFile("Histo_simu_new/Energy_MC_simu_K_40_recOC_new_eres_53_gain_221_wubc.root", "READ");
  // histo_file_K->cd();
  // TH3D* MC_K_40 = (TH3D*)histo_file_K->Get("MC_simu_K_40_recOC_new_ubc");



  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();

  int lim = 0;
  int om = 0;
  double param1 = 0;
  double param2 = 0;
  double param3 = 0;
  double_t error1 =0;
  double_t error2 =0;
  double_t error3 = 0;
  double Chi2NDF = 0;
  double mean_erf = 0;
  double sigma_erf = 0;
  double like =0;
  double* fftab = new double[8];
  double* tab = new double[7];
  float gain = 0;
  float eres = 0;

  TFile *newfile = new TFile("histo_test2.root", "RECREATE");
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("Chi2NDF", &Chi2NDF);
  Result_tree.Branch("param1", &param1);
  Result_tree.Branch("param2", &param2);
  Result_tree.Branch("param3", &param3);
  Result_tree.Branch("error1", &error1);
  Result_tree.Branch("error2", &error2);
  Result_tree.Branch("error3", &error3);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("eres", &eres);
  Result_tree.Branch("om", &om);
  Result_tree.Branch("lim", &lim);
  Result_tree.Branch("mean_erf", &mean_erf);
  Result_tree.Branch("sigma_erf", &sigma_erf);
  Result_tree.Branch("likelihood", &like);

  for (om = 19; om < 20; om++) 
    {
      if (((om%13)==12) || ((om%13)==0))continue; 
      //if ((om%13)==12)continue;
      //if ((om%13)==0)continue;
    
      TH1D* spectre_om = NULL;
      spectre_om = spectre_chooser(om);

      for (int bin =1; bin < lim; bin++) {
	spectre_om->SetBinContent(bin, 0);
      }
      if ((spectre_om->GetEntries() < 100) || (spectre_om->GetMean(1) < 20000)){
	delete spectre_om;
	om++;
	spectre_om = spectre_chooser(om);
      }
      if (charge_valeur_fit[om] != -1) {
	if (spectre_om->GetEntries() > 100){
	  tab = om_gain_fit(om);
	  if (tab[0] == 0) {
	    mean_erf = 30000;
	    sigma_erf = 0;
	  }
	  else{
	    sigma_erf = 1/tab[1];
	    mean_erf = 1/tab[0];
	  }
	}
      }
 
      lim = 105; 
      for (int bin = 1; bin < lim; bin++){
	spectre_om->SetBinContent(bin,0);
      }
      if (spectre_om->GetEntries() > 100){
	for (int eres_count = 0; eres_count < 53; eres_count++) {
	  for (int gain_count = 0; gain_count <221; gain_count++) {
	    std::cout << (gain_bin_min + gain_bin_width*(gain_count-1))<< "   et    lim_inf = " << 1/charge_valeur_fit[om]*0.9 << "   sup  ="  << 1/charge_valeur_fit[om]*1.1 << '\n';
	    if ((1/(gain_bin_min + gain_bin_width*(gain_count-1)) > charge_valeur_fit[om]*0.9) && (1/(gain_bin_min + gain_bin_width*(gain_count-1))<charge_valeur_fit[om]*1.1))	
	      // std::cout << (gain_bin_min + gain_bin_width*(gain_count-1))<< "   et    lim_inf = " << mean_erf*0.8 << "   sup  ="  << mean_erf*1.2 << '\n';
	    // if (((gain_bin_min + gain_bin_width*(gain_count-1)) > (mean_erf*0.85)) && (1/(gain_bin_min + gain_bin_width*(gain_count-1)) < (mean_erf*1.15)))
	      {
		gain = 1/(gain_bin_min + gain_bin_width*(gain_count-1));
		eres = eres_bin_min + eres_bin_width*(eres_count-1);

		std::cout << "eres_count = " << eres_count << '\n';
		std::cout << "gain_count = " << gain_count << '\n';
		TH1D *mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", eres_count, eres_count, gain_count, gain_count);    // first MC histogram
		TH1D *mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
		TH1D *mc2 = MC_K_40->ProjectionZ("Charge_K_40", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
		
		fractionfitter(om, draw, gain, eres, spectre_om, mc0, mc1, mc2, fftab);
		
		cout << fftab[0] << endl;

		Chi2NDF = fftab[0];
		like = fftab[1];
		param1 = fftab[2];
		param2 = fftab[3];
		param3 = fftab[4];
		error1 = fftab[5];
		error2 = fftab[6];
		error3 = fftab[7];		
		Result_tree.Fill();
		
		delete mc0;
		delete mc1;
		delete mc2;
	      }
	  }
	}
      }
      delete spectre_om;
    }
  
  // std::cout << "****************************************************" << '\n';
  // std::cout << "The best X2 is : " <<  Result_tree.GetMinimum ("Chi2NDF") << '\n';
  // std::cout << "****************************************************" << '\n';
  
  newfile->cd();
  Result_tree.Write();
  newfile->Close();
  
}

void distrib(string name) {
  TFile *eff_file = new TFile(Form("histo_fit/histo_%s.root", name.c_str()), "READ");
  std::ofstream outFile("Best_Khi2.txt");
  int om =0;
  double Chi2NDF = 0;
  float gain = 0;
  float eres = 0;
  TTree* eff_tree = (TTree*)eff_file->Get("Result_tree");
  eff_tree->SetBranchStatus("*",0);
  eff_tree->SetBranchStatus("om",1);
  eff_tree->SetBranchAddress("om", &om);
  eff_tree->SetBranchStatus("Chi2NDF",1);
  eff_tree->SetBranchAddress("Chi2NDF", &Chi2NDF);
  eff_tree->SetBranchStatus("gain",1);
  eff_tree->SetBranchAddress("gain", &gain);
  eff_tree->SetBranchStatus("eres",1);
  eff_tree->SetBranchAddress("eres", &eres);
  double n_entry = 0;
  
  double test = 10;
  if (name.compare("MW8") == 0){
    TH1D* dist = new TH1D ("distribution Chi2","distribution Chi2 OM MW 8p",100, 0, 2);
    dist->GetXaxis()->SetTitle("Khi2");
    for (int j = 1; j < 520; j++) {
      if (((j%13)!=0) && ((j%13)!=12)){
	for (double i = 0; i < eff_tree->GetEntries(); i++) {
	  eff_tree->GetEntry(i);
	  if (om == j) {
	    if (Chi2NDF < test) {
	      test = Chi2NDF;
	      n_entry = i;
	    }
	  }
	}
      }
      dist->Fill(test);
      test =10;
      if (((j%13)!=0) && ((j%13)!=12)){
	eff_tree->GetEntry(n_entry);
	outFile << om << "\t" << gain << "\t" << eres << endl;  
      }
    }
    dist->Draw();
  }
  else if (name.compare("MW5") == 0){
    TH1D* dist = new TH1D ("distribution Chi2","distribution Chi2 OM 5p MW",100, 0, 2);
    dist->GetXaxis()->SetTitle("Khi2");  
    for (int j = 0; j < 520; j++) {
      if (((j%13)==0) || ((j%13)==12)){
	for (double i = 0; i < eff_tree->GetEntries(); i++) {
	  eff_tree->GetEntry(i);
	  if (om == j) {
	    if (Chi2NDF < test) {
	      test = Chi2NDF;
	      n_entry = i;
	    }
	  }
	}
      }
      dist->Fill(test);
      cout << j << endl;
      test = 10;
      if (((j%13)==0) || ((j%13)==12)){
	eff_tree->GetEntry(n_entry);
	outFile << om << "\t" << gain << "\t" << eres << endl;  
      }
    }
    dist->Draw();
  }
  else if (name.compare("XW") == 0){
    TH1D* dist = new TH1D ("distribution Chi2","distribution Chi2 OM XW",100, 0, 2);
    dist->GetXaxis()->SetTitle("Khi2");
    for (int j = 520; j < 648; j++) {
      for (double i = 0; i < eff_tree->GetEntries(); i++) {
	eff_tree->GetEntry(i);
	if (om == j) {
	  if (Chi2NDF < test) {
	    test = Chi2NDF;
	    n_entry = i;
	  }
	}
      }
      dist->Fill(test);
      test =10;
      eff_tree->GetEntry(n_entry);
      outFile << om << "\t" << gain << "\t" << eres << endl;  
    }
    dist->Draw();
  }
  else if (name.compare("GV") == 0){
    TH1D* dist = new TH1D ("distribution Chi2","distribution Chi2 OM GV",100, 0, 11);
    dist->GetXaxis()->SetTitle("Khi2");
    for (int j = 648; j < 712; j++) {
      for (double i = 0; i < eff_tree->GetEntries(); i++) {
	eff_tree->GetEntry(i);
	if (om == j) {
	  if (Chi2NDF < test) {
	    test = Chi2NDF;
	    n_entry = i;
	  }
	}
      } 
      dist->Fill(test);
      test =10;
      eff_tree->GetEntry(n_entry);
      outFile << om << "\t" << gain << "\t" << eres << endl;  
    }
    dist->Draw();
  }
}

int main(int argc, char const *argv[])
{
  kolmo(1);
  return 0;
}
