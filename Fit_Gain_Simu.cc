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
const float gain_bin_min = 20000;
const float gain_bin_max = 75000;
const float gain_bin_width = (gain_bin_max-gain_bin_min)/(gain_n_bin-1);

const int charge_n_bin = 1024;
const float charge_bin_min = 0e-05;
const float charge_bin_max = 200000;
const float charge_bin_width = (charge_bin_max-charge_bin_min)/(charge_n_bin-1);

void test()
{

  for(float gain = gain_bin_min; gain<= gain_bin_max; gain+= gain_bin_width) {
    std::cout << "gain = " << gain  << " bin = " << gain/gain_bin_width - 79 <<'\n';
  }
  // for (float eres = eres_bin_min; eres<= eres_bin_max; eres+= eres_bin_width) {
  //   std::cout << "eres = " << eres << " bin = " << eres/eres_bin_width - 27 <<'\n';
  // }
}

TH1D* spectre_charge_it(int om_number )
{
    TFile *file = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_422.root", "READ");
    gROOT->cd();
    TH2F* charge = (TH2F*)file->Get("charge");

    TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
    file->Close();
    return spectre_charge;
  }

TH1D* spectre_charge_XW(int om_number )
{
    TFile *file = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_energie_449.root", "READ");
    gROOT->cd();
    TH2F* charge = (TH2F*)file->Get("histo_pm_charge");

    TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
    file->Close();
    return spectre_charge;
  }

TH1D* spectre_charge_fr(int om_number )
{
    TFile *file = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_energie_435.root", "READ");
    gROOT->cd();
    TH2F* charge = (TH2F*)file->Get("histo_pm_charge");

    TH1D* spectre_charge = charge->ProjectionY(Form("charge%03d",om_number), om_number+1, om_number+1);
    file->Close();
    return spectre_charge;
  }

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

  TFile *file = new TFile(Form("Histo_simu/%s.root", name.c_str()), "READ");

  std::vector<int> *om_num = new std::vector<int>;
  std::vector<double> *energy = new std::vector<double>;
  std::vector<int> *om_num_new = new std::vector<int>;
  std::vector<double> *energy_ubc = new std::vector<double>;

  TRandom3 rando;

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_id",1);
  tree->SetBranchAddress("om_id", &om_num);
  tree->SetBranchStatus("energy",1);
  tree->SetBranchAddress("energy", &energy);
  tree->SetBranchStatus("om_id_new",1);
  tree->SetBranchAddress("om_id_new", &om_num_new);
  tree->SetBranchStatus("energy_ubc",1);
  tree->SetBranchAddress("energy_ubc", &energy_ubc);

  TH3D* MC_Simu = new TH3D(Form("MC_Simu_%s", name.c_str()), Form("MC_Simu_%s", name.c_str()),
                                eres_n_bin, eres_bin_min - (eres_bin_width/2), eres_bin_max + (eres_bin_width/2),
                                gain_n_bin, gain_bin_min - (gain_bin_width/2), gain_bin_max + (gain_bin_width/2),
                                charge_n_bin, charge_bin_min, charge_bin_max);
  TH3D* MC_Simu_ubc = new TH3D(Form("MC_%s_ubc", name.c_str()), Form("MC_Simu_%s_ubc", name.c_str()),
                                eres_n_bin, eres_bin_min - (eres_bin_width/2), eres_bin_max + (eres_bin_width/2),
                                gain_n_bin, gain_bin_min - (gain_bin_width/2), gain_bin_max + (gain_bin_width/2),
                                charge_n_bin, charge_bin_min, charge_bin_max);

  // for (int i = 0; i < 10000; i++) {
  // for (int i = 0; i < 1000000; i++) {
  // // for (int i = 0; i < tree->GetEntries(); i++) {
  //   double E_kolmo =0;
  //   tree->GetEntry(i);
  //
  //   for (int j = 0; j < energy->size(); j++) {
  //     if (om_num->at(j) < 520) {
  //       for (float eres = eres_bin_min; eres<= eres_bin_max; eres+= eres_bin_width) {
  //         float Evis = rando.Gaus(energy->at(j), (eres/235.482)*sqrt(energy->at(j)));
  //         for(float gain = gain_bin_min; gain<= gain_bin_max; gain+= gain_bin_width) {
  //             E_kolmo = gain * Evis;
  //             MC_Simu->Fill(eres, gain, E_kolmo);
  //           }
  //         }
  //       }
  //     }
  //   }
  for (int i = 0; i < 1000000; i++) {
    if (i%100000 == 0) {
      std::cout << "/* entry = " << i << '\n';
    }
    double E_kolmo =0;
    tree->GetEntry(i);

    for (int j = 0; j < energy_ubc->size(); j++) {
      // if (om_num_new->at(j) < 520 && (om_num_new->at(j)%13)!=0 && (om_num_new->at(j)%13)!=12) {
      if (om_num_new->at(j) < 520 && ((om_num_new->at(j)%13)==0 || (om_num_new->at(j)%13)==12)) {
      // if (om_num_new->at(j) < 648 && om_num_new->at(j) >519) {
      // if (om_num_new->at(j) > 647){
        for (float eres = eres_bin_min; eres<= eres_bin_max; eres+= eres_bin_width) {
          float Evis = rando.Gaus(energy_ubc->at(j), (eres/235.482)*sqrt(energy_ubc->at(j)));
          for(float gain = gain_bin_min; gain<= gain_bin_max; gain+= gain_bin_width) {
              E_kolmo = gain * Evis;
              MC_Simu_ubc->Fill(eres, gain, E_kolmo);
            }
          }
        }
      }
  }
  TFile *newfile = new TFile(Form("Histo_simu/MC_%s_eres_%d_gain_%d_OC_5p.root", name.c_str(), eres_n_bin, gain_n_bin), "RECREATE");
  newfile->cd();
  MC_Simu_ubc->Write();
  MC_Simu->Write();
  newfile->Close();

  return MC_Simu;
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

  // TTree tree_fit("tree_fit","");
  // tree_fit.Branch("n_evt", &n_evt););
  // tree_fit.Branch("mean", &mean);
  // tree_fit.Branch("sigma", &sigma);
  // tree_fit.Branch("nbg", &nbg);
  // tree_fit.Branch("chi", &chi);
  // tree_fit.Branch("ndf", &ndf);
  // tree_fit.Branch("chin", &chin);

  double* tab = new double[7];
  TCanvas* canvas = new TCanvas;
  canvas->SetLogy();

  // TH1D* spectre_om = spectre_charge(om);
  // spectre_om->Draw("same");
  TF1* f_ComptonEdgePoly = new TF1 ("f_ComptonEdgePoly","[0]/2.0*(1+TMath::Erf(([1]-x)/(TMath::Sqrt(2)*[2])))+[3]*x", 40000, 90000);
  f_ComptonEdgePoly->SetParNames("N_evt","Mean","Sigma","Nbg" );

if (om == 1111 )        //om multiple de (13)-1
  {
    TFile *histo_file = new TFile("histo_mystere/new_histo_mystere_1.root", "READ");
    TH1D* spectre_om = (TH1D*)histo_file->Get("mc_tot");
    spectre_om->Draw("same");
    spectre_om->SetLineColor(kBlack);
    f_ComptonEdgePoly->SetParameters(120, 72367, 6893, 3.91e-5);
    f_ComptonEdgePoly->SetRange(60000,100000);
    f_ComptonEdgePoly->Draw("same");
    spectre_om->Fit(f_ComptonEdgePoly, "RQ");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ");
  }
  else if (om == 1112 )        //om multiple de (13)-1
  {
    TFile *histo_file = new TFile("histo_mystere/new_histo_mystere_2.root", "READ");
    TH1D* spectre_om = (TH1D*)histo_file->Get("mc_tot");
    spectre_om->Draw("same");
    spectre_om->SetLineColor(kBlack);
    f_ComptonEdgePoly->SetParameters(120, 130000, 6893, 3.91e-5);
    f_ComptonEdgePoly->SetRange(120000,160000);
    f_ComptonEdgePoly->Draw("same");
    spectre_om->Fit(f_ComptonEdgePoly, "RQ");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ");
  }
  else if (om == 1113 )        //om multiple de (13)-1
  {
    TFile *histo_file = new TFile("histo_mystere/new_histo_mystere_3.root", "READ");
    TH1D* spectre_om = (TH1D*)histo_file->Get("mc_tot");
    spectre_om->Draw("same");
    spectre_om->SetLineColor(kBlack);
    f_ComptonEdgePoly->SetParameters(120, 90000, 6893, 3.91e-5);
    f_ComptonEdgePoly->SetRange(80000,120000);
    f_ComptonEdgePoly->Draw("same");
    spectre_om->Fit(f_ComptonEdgePoly, "RQ");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ");
  }
  else if (om == 1114 )        //om multiple de (13)-1
  {
    TFile *histo_file = new TFile("histo_mystere/new_histo_mystere_4.root", "READ");
    TH1D* spectre_om = (TH1D*)histo_file->Get("mc_tot");
    spectre_om->Draw("same");
    spectre_om->SetLineColor(kBlack);
    f_ComptonEdgePoly->SetParameters(120, 72367, 6893, 3.91e-5);
    f_ComptonEdgePoly->SetRange(60000,100000);
    f_ComptonEdgePoly->Draw("same");
    spectre_om->Fit(f_ComptonEdgePoly, "RQ");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ");
    f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
    spectre_om->Fit(f_ComptonEdgePoly, "RQ");
  }
  // if ((om % 13) == 12 )        //om multiple de (13)-1
  // {
  //   f_ComptonEdgePoly->SetParameters(120, 72367, 6893, 3.91e-5);
  //   f_ComptonEdgePoly->SetRange(6000,100000);
  //   f_ComptonEdgePoly->Draw("same");
  //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
  //   f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
  //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
  //   f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
  //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
  // }
  // else if ((om % 13) == 0)       //om multiple de 13
  // {
  //   f_ComptonEdgePoly->SetParameters(112, 68168, 5604, 1.2e-05);
  //   f_ComptonEdgePoly->SetRange(50000,100000);
  //   f_ComptonEdgePoly->Draw("same");
  //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
  //   f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
  //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
  //   f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-1.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+5*f_ComptonEdgePoly->GetParameter(2));
  //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
  // }
  // else         //om normaux (8pouces)
  // {
  //   f_ComptonEdgePoly->SetParameters(111, 60978, 3787, 4.19e-05);
  //   f_ComptonEdgePoly->SetRange(55000,100000);
  //   f_ComptonEdgePoly->Draw("same");
  //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
  //   f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-2.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+7.5*f_ComptonEdgePoly->GetParameter(2));
  //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
  //   f_ComptonEdgePoly->SetRange(f_ComptonEdgePoly->GetParameter(1)-3.5*f_ComptonEdgePoly->GetParameter(2),f_ComptonEdgePoly->GetParameter(1)+7.5*f_ComptonEdgePoly->GetParameter(2));
  //   spectre_om->Fit(f_ComptonEdgePoly, "RQ");
  // }

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

  tab[0] = mean;
  tab[1] = sigma;
  tab[2] = chin;
  tab[3] = chi;
  tab[4] = ndf;
  tab[5] = n_evt;
  tab[6] = nbg;

std::cout << tab[0] << '\n';
std::cout << tab[2] << '\n';
return tab;
}

void tuning_param()
{
  double* tab = new double[7];
  double mean_n = 0;
  double mean_13 = 0;
  double mean_0 = 0;
  double sigma_n = 0;
  double sigma_13 = 0;
  double sigma_0 = 0;
  double N_evt_n = 0;
  double N_evt_13 = 0;
  double N_evt_0 = 0;
  double NBG_n = 0;
  double NBG_13 = 0;
  double NBG_0 = 0;
  int compteur = 0;

  for (size_t om = 260; om < 520; om++) {
    tab = om_gain_fit(om);
    if ((om % 13) == 12 ) {               //om multiple de 13 - 1
      mean_0 += tab[0];
      sigma_0 += tab[1];
      N_evt_0 += tab[5];
      NBG_0 += tab[6];
      if (tab[0] > 10){
        compteur++;
      }
    }
    else if ((om % 13) == 0)       //om multiple de 13
    {
      mean_13 += tab[0];
      sigma_13 += tab[1];
      N_evt_13 += tab[5];
      NBG_13 += tab[6];
    }
    else{
      mean_n += tab[0];
      sigma_n += tab[1];
      N_evt_n += tab[5];
      NBG_n += tab[6];
    }
  }
  std::cout << "N_evt_n = " << (N_evt_0/compteur)*2.6 <<'\n';
  std::cout << "mean_n = " << (mean_0/compteur)*2.6 <<'\n';
  std::cout << "sigma_n = " << (sigma_0/compteur)*2.6 <<'\n';
  std::cout << "Nbg_n = " << (NBG_0/compteur)*2.6 <<'\n';
  std::cout << compteur << '\n';
}

void kolmo()
{
  std::ifstream  charge("/home/aguerre/Bureau/Thèse/Fit_Gain_Simu/gain_data/Resultat_mean_energie_charge_total.txt");
  float charge_valeur_fit [712];
  memset(charge_valeur_fit, -1, 712*sizeof(float));
  int charge_om_num;
  while (charge >> charge_om_num)
  {
    charge >> charge_valeur_fit[charge_om_num];
  }

  TFile *histo_file_Tl = new TFile("Histo_simu/MC_simu_Tl_208_recOC_eres_53_gain_221_OC.root", "READ");
  histo_file_Tl->cd();
  TH3D* MC_Tl_208 = (TH3D*)histo_file_Tl->Get("MC_Simu_simu_Tl_208_recOC_ubc");

  TFile *histo_file_Bi = new TFile("Histo_simu/MC_simu_Bi_214_recOC_eres_53_gain_221_OC_8p.root", "READ");
  histo_file_Bi->cd();
  TH3D* MC_Bi_214 = (TH3D*)histo_file_Bi->Get("MC_Simu_simu_Bi_214_recOC_ubc");

  TFile *histo_file_K = new TFile("Histo_simu/MC_simu_K_40_recOC_eres_53_gain_221_OC.root", "READ");
  histo_file_K->cd();
  TH3D* MC_K_40 = (TH3D*)histo_file_K->Get("MC_Simu_simu_K_40_recOC_ubc");

  // TFile *histo_file_Tl = new TFile("Histo_simu/MC_Simu_Tl_208_eres_53_gain_221.root", "READ");
  // histo_file_Tl->cd();
  // TH3D* MC_Tl_208 = (TH3D*)histo_file_Tl->Get("MC_Simu_Tl_208");
  //
  // TFile *histo_file_Bi = new TFile("Histo_simu/MC_Simu_Bi_214_eres_53_gain_221.root", "READ");
  // histo_file_Bi->cd();
  // TH3D* MC_Bi_214 = (TH3D*)histo_file_Bi->Get("MC_Simu_Bi_214");
  //
  // TFile *histo_file_K = new TFile("Histo_simu/MC_Simu_K_40_eres_53_gain_221.root", "READ");
  // histo_file_K->cd();
  // TH3D* MC_K_40 = (TH3D*)histo_file_K->Get("MC_Simu_K_40");

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();

  // TH3D* MC_Tl_208 = MC_Simu("Tl_208");
  // TH3D* MC_Bi_214 = MC_Simu("Bi_214");
  // TH3D* MC_K_40 = MC_Simu("K_40");

  int lim = 0;
  int om = 0;
  double param1 = 0;
  double param2 = 0;
  double param3 = 0;
  double_t error1 =0;
  double_t error2 =0;
  double_t error3 = 0;
  double Chi2NDF = 0;
  double activity_Tl = 0;
  double activity_Bi = 0;
  double activity_K = 0;
  double result_0_scale;
  double result_1_scale;
  double result_2_scale;
  double integrale_gauche_Tl;
  double integrale_droite_Tl;
  double integrale_gauche_Bi;
  double integrale_droite_Bi;
  double integrale_gauche_K;
  double integrale_droite_K;
  double mean_erf = 0;
  double sigma_erf = 0;

  double* tab = new double[7];
  float gain = 0;
  float eres = 0;
  double eff_tot;
  double eff_cut;

  TFile *eff_file = new TFile("histo_kolmo/histo_donee/histo_charge_amplitude_energie_435.root", "READ");

  // TTree* eff_tree = (TTree*)eff_file->Get("new_tree");
  // eff_tree->SetBranchStatus("*",0);
  // eff_tree->SetBranchStatus("eff_tot",1);
  // eff_tree->SetBranchAddress("eff_tot", &eff_tot);
  // eff_tree->SetBranchStatus("eff_cut",1);
  // eff_tree->SetBranchAddress("eff_cut", &eff_cut);
  TFile *newfile = new TFile("htest.root", "RECREATE");

  // TFile *newfile = new TFile("histo_kolmo/Simu_kolmo_MW8_new.root", "RECREATE");
  TH2D* Chi2 = new TH2D("Chi2", "Chi2", eres_n_bin-1, eres_bin_min, eres_bin_max, gain_n_bin-1, gain_bin_min, gain_bin_max);
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("Chi2NDF", &Chi2NDF);
  Result_tree.Branch("param1", &param1);
  Result_tree.Branch("param2", &param2);
  Result_tree.Branch("param3", &param3);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("eres", &eres);
  Result_tree.Branch("om", &om);
  Result_tree.Branch("lim", &lim);
  Result_tree.Branch("activity_Tl", &activity_Tl);
  Result_tree.Branch("activity_Bi", &activity_Bi);
  Result_tree.Branch("activity_K", &activity_K);
  Result_tree.Branch("integrale_droite_Tl", &integrale_droite_Tl);
  Result_tree.Branch("integrale_gauche_Tl", &integrale_gauche_Tl);
  Result_tree.Branch("integrale_droite_Bi", &integrale_droite_Bi);
  Result_tree.Branch("integrale_gauche_Bi", &integrale_gauche_Bi);
  Result_tree.Branch("integrale_droite_K", &integrale_droite_K);
  Result_tree.Branch("integrale_gauche_K", &integrale_gauche_K);
  Result_tree.Branch("mean_erf", &mean_erf);
  Result_tree.Branch("sigma_erf", &sigma_erf);

  for (om = 100; om < 101; om++) {

    if (charge_valeur_fit[om] != -1) {
      tab = om_gain_fit(om);
      float limi = 0.8;
      float lims = 1.3;
      if (om > 647) {
        mean_erf = 15000;
      }
       if (tab[0] == 0) {
        mean_erf = 30000;
        sigma_erf = 0;
        limi = 0.6;
        lims = 1.5;
        std::cout << "le fit erf n'a pas bien fonctionné !" << '\n';
      }
      else{
        sigma_erf = 1/tab[1];
        mean_erf = 1/tab[0];
      }
    TH1D* spectre_om = spectre_charge_it(om);
    if (om <520 && om >259) {
      TH1D* spectre_om = spectre_charge_fr(om);
    }
    if (om >519) {
      TH1D* spectre_om = spectre_charge_XW(om);
    }

      // lim = (1/(200000*mean_erf/1024));

      lim = 105;

      for (int bin =1; bin < lim; bin++) {
        spectre_om->SetBinContent(bin, 0);
      }

      for (int eres_count = 0; eres_count < 53; eres_count++) {
        for (int gain_count = 0; gain_count < 221; gain_count++) {
          // std::cout << (gain_bin_min + gain_bin_width*(gain_count-1))<< "   et    lim_inf = " << mean_erf*0.8 << "   sup  ="  << mean_erf*1.2 << '\n';
          // if (om> 647)
          if ((1/(gain_bin_min + gain_bin_width*(gain_count-1)) > charge_valeur_fit[om]*0.8) && (1/(gain_bin_min + gain_bin_width*(gain_count-1))<charge_valeur_fit[om]*1.2))
          // if (((gain_bin_min + gain_bin_width*(gain_count-1)) > (mean_erf*0.8)) && (1/(gain_bin_min + gain_bin_width*(gain_count-1)) < (mean_erf*1.5)))
          {
            std::cout << eres_count << '\n';
            std::cout << gain_count << '\n';
            TH1D *mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", eres_count, eres_count, gain_count, gain_count);    // first MC histogram
            TH1D *mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
            TH1D *mc2 = MC_K_40->ProjectionZ("Charge_K_40", eres_count, eres_count, gain_count, gain_count);    // second MC histogram

            integrale_gauche_Tl = mc0->Integral(0, lim);
            integrale_droite_Tl = mc0->Integral(lim+1, 1024);
            integrale_gauche_Bi = mc1->Integral(0, lim);
            integrale_droite_Bi = mc1->Integral(lim+1, 1024);
            integrale_gauche_K = mc2->Integral(0, lim);
            integrale_droite_K = mc2->Integral(lim+1, 1024);


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

              mc1->Draw();
              mc2->Draw("same");
              mc0->Draw("same");

              Chi2NDF = (fit->GetChisquare())/(fit->GetNDF());

              TH1D* result = (TH1D*) fit->GetPlot();
              TH1D* result_0 = (TH1D*) fit->GetMCPrediction(0);
              TH1D* result_1 = (TH1D*) fit->GetMCPrediction(1);
              TH1D* result_2 = (TH1D*) fit->GetMCPrediction(2);

              gain = 1/(gain_bin_min + gain_bin_width*(gain_count-1));
              eres = eres_bin_min + eres_bin_width*(eres_count-1);

              std::cout << "om : " << om << "   gain : " << gain << "    eres : " << eres << '\n';

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

              //
              // if (om == 680 && gain < 30000)
              // {
              //   canvas->SaveAs(Form("Fit_kolmo/GV/fit_kolmo_om_%d_eres_%d_gain_%d.png", om, eres_count, gain_count));
              // }
              // if (om < 390)
              // {
              //   Chi2->SetBinContent(eres_count, gain_count, Chi2NDF);
              // }
              std::cout << "kolmo = " << spectre_om->KolmogorovTest(result) << '\n';

              Result_tree.Fill();
              delete canvas;
              delete fit;
            }

            delete mc0;
            delete mc1;
            delete mc2;
          }
        }
      }
      // if ((om%13)==11  ) {
      //   om +=2;
      // }
    }
  }

 std::cout << "****************************************************" << '\n';
 std::cout << "The best X2 is : " <<  Result_tree.GetMinimum ("Chi2NDF") << '\n';
 std::cout << "****************************************************" << '\n';
 //
 // TH2D* parabole = new TH2D("parabole", "parabole", 500, 6, 21, 1000, 0, 8);
 // double good_gain = 0;
 //
 // for (int i = 0; i < Result_tree.GetEntries(); i++) {
 //   Result_tree.GetEntry(i);
 //
 //   if (Chi2NDF == Result_tree.GetMinimum ("Chi2NDF")) {
 //     good_gain = gain;
 //   }
 // }
 //
 // for (int i = 0; i < Result_tree.GetEntries(); i++) {
 //   Result_tree.GetEntry(i);
 //
 //   if (gain == good_gain) {
 //     parabole->Fill(eres, Chi2NDF);
 //   }
 // }
 newfile->cd();

 // parabole->Write();
 Result_tree.Write();

 newfile->Close();

}

void distrib() {
  TFile *eff_file = new TFile("histo_kolmo/Simu_kolmo_MW8_old.root", "READ");
  int om =0;
  double Chi2NDF = 0;
  TTree* eff_tree = (TTree*)eff_file->Get("Result_tree");
  eff_tree->SetBranchStatus("*",0);
  eff_tree->SetBranchStatus("om",1);
  eff_tree->SetBranchAddress("om", &om);
  eff_tree->SetBranchStatus("Chi2NDF",1);
  eff_tree->SetBranchAddress("Chi2NDF", &Chi2NDF);
  TH1D* dist = new TH1D ("distribution Chi2","distribution Chi2 OM 8p",100, 0, 2);

  double test =2;

  for (int j = 390; j < 400; j++) {
    double min =0;
    for (double i = 0; i < eff_tree->GetEntries(); i++) {
      eff_tree->GetEntry(i);

      if (om == j) {
        if (Chi2NDF < test) {
          test = Chi2NDF;
        }
      }

    }

  dist->Fill(test);
  test =2;
  }
  dist->Draw();

}

void kolmo_mystere()
{
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2();

  TFile *histo_file_Tl = new TFile("Histo_simu/MC_Simu_Tl_208_eres_53_gain_221.root", "READ");
  histo_file_Tl->cd();
  TH3D* MC_Tl_208 = (TH3D*)histo_file_Tl->Get("MC_Simu_Tl_208");

  TFile *histo_file_Bi = new TFile("Histo_simu/MC_Simu_Bi_214_eres_53_gain_221.root", "READ");
  histo_file_Bi->cd();
  TH3D* MC_Bi_214 = (TH3D*)histo_file_Bi->Get("MC_Simu_Bi_214");

  TFile *histo_file_K = new TFile("Histo_simu/MC_Simu_K_40_eres_53_gain_221.root", "READ");
  histo_file_K->cd();
  TH3D* MC_K_40 = (TH3D*)histo_file_K->Get("MC_Simu_K_40");

  // TH3D* MC_Tl_208 = MC_Simu("Tl_208");
  // std::cout << "ok Tl 208" << '\n';
  // TH3D* MC_Bi_214 = MC_Simu("Bi_214");
  // std::cout << "ok Bi 214" << '\n';
  // TH3D* MC_K_40 = MC_Simu("K_40");
  // std::cout << "ok K 40" << '\n';

  int lim = 0;
  double param1 = 0;
  double param2 = 0;
  double param3 = 0;
  double error1 =0;
  double error2 =0;
  double error3 = 0;
  double Chi2NDF = 0;
  double activity_Tl = 0;
  double activity_Bi = 0;
  double activity_K = 0;
  double integrale_gauche_Tl = 0;
  double integrale_droite_Tl = 0;
  double integrale_gauche_Bi = 0;
  double integrale_droite_Bi = 0;
  double integrale_gauche_K = 0;
  double integrale_droite_K = 0;
  double result_0_scale = 0;
  double result_1_scale = 0;
  double result_2_scale = 0;
  double result_0_scale_er = 0;
  double result_1_scale_er = 0;
  double result_2_scale_er = 0;
  double total_hit_Tl = 0;
  double total_hit_Bi = 0;
  double total_hit_K = 0;
  double error_gain = 0;
  double error_res = 0;
  double mean_erf = 0;
  double sigma_erf = 0;

  int lim_tree = 0;
  float gain = 0;
  float eres = 0;

  TFile *newfile = new TFile("histo_kolmo/Simu_mystere_4.root", "RECREATE");
  TH2D* Chi2 = new TH2D("Chi2", "Chi2", eres_n_bin-1, eres_bin_min, eres_bin_max, gain_n_bin-1, gain_bin_min, gain_bin_max);
  TTree Result_tree("Result_tree","");
  Result_tree.Branch("Chi2NDF", &Chi2NDF);
  Result_tree.Branch("param1", &param1);
  Result_tree.Branch("param2", &param2);
  Result_tree.Branch("param3", &param3);
  Result_tree.Branch("gain", &gain);
  Result_tree.Branch("eres", &eres);
  Result_tree.Branch("lim_tree", &lim_tree);
  Result_tree.Branch("activity_Tl", &activity_Tl);
  Result_tree.Branch("activity_Bi", &activity_Bi);
  Result_tree.Branch("activity_K", &activity_K);
  Result_tree.Branch("integrale_droite_Tl", &integrale_droite_Tl);
  // Result_tree.Branch("integrale_gauche_Tl", &integrale_gauche_Tl);
  Result_tree.Branch("integrale_droite_Bi", &integrale_droite_Bi);
  // Result_tree.Branch("integrale_gauche_Bi", &integrale_gauche_Bi);
  Result_tree.Branch("integrale_droite_K", &integrale_droite_K);
  // Result_tree.Branch("integrale_gauche_K", &integrale_gauche_K);
  Result_tree.Branch("total_hit_Tl", &total_hit_Tl);
  Result_tree.Branch("total_hit_Bi", &total_hit_Bi);
  Result_tree.Branch("total_hit_K", &total_hit_K);
  Result_tree.Branch("error_gain", &error_gain);
  Result_tree.Branch("error_res", &error_res);
  Result_tree.Branch("error1", &error1);
  Result_tree.Branch("error2", &error2);
  Result_tree.Branch("error3", &error3);
  Result_tree.Branch("result_0_scale_er", &result_0_scale_er);
  Result_tree.Branch("result_1_scale_er", &result_1_scale_er);
  Result_tree.Branch("result_2_scale_er", &result_2_scale_er);
  Result_tree.Branch("result_0_scale", &result_0_scale);
  Result_tree.Branch("result_1_scale", &result_1_scale);
  Result_tree.Branch("result_2_scale", &result_2_scale);
  Result_tree.Branch("mean_erf", &mean_erf);
  Result_tree.Branch("sigma_erf", &sigma_erf);

  double* tab = new double[7];

  for (int i = 4; i < 5; i++) {

  TFile *file = new TFile(Form("histo_mystere/new_histo_mystere_%i.root", i), "READ");
  TH1D* spectre_om = (TH1D*)file->Get("mc_tot");

  lim = (195);

  // int gain_count = 74;

    int plus = 1;
    tab = om_gain_fit(1110 + i);
      float limi = 0.8;
      float lims = 1.3;
      if (tab[0] == 0) {
        mean_erf = 40000;
        limi = 0.6;
        lims = 1.5;
        std::cout << "le fit erf n'a pas bien fonctionné !" << '\n';
      }
      else{
        mean_erf = tab[0];
      }
      std::cout << "mean erf = " << mean_erf <<'\n';

          // for ( lim = 0; lim <225; lim+=plus) {
          // std::cout << "lim = " << lim << '\n';

            for (int bin =1; bin < lim; bin++) {
              spectre_om->SetBinContent(bin, 0);
            }

        for (int gain_count = 1; gain_count <221; gain_count++) {
          for (int eres_count = 1; eres_count < 52; eres_count++) {

            // float P =(1/40000.0*0.5);
            // float r = (1/40000.0*1.5);
            float P =(38250.0*0.6);
            float r = (38250.0*1.4);
            lim_tree = lim;
            std::cout << (gain_bin_min + gain_bin_width*(gain_count-1))<< "   et    lim_inf = " << (mean_erf*0.8) << "   sup  ="  <<(mean_erf*1.5)<< '\n';

          if (((gain_bin_min + gain_bin_width*(gain_count-1)) < (mean_erf*1.2)) && ((gain_bin_min + gain_bin_width*(gain_count-1)) > (mean_erf*0.))){
          // if ((1/(gain_bin_min + gain_bin_width*(gain_count-1))> P) && (1/(gain_bin_min + gain_bin_width*(gain_count-1))< r))
          // {
            TH1D *mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", eres_count, eres_count, gain_count, gain_count);    // first MC histogram
            TH1D *mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
            TH1D *mc2 = MC_K_40->ProjectionZ("Charge_K_40", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
            TH1D *mc0_full = MC_Tl_208->ProjectionZ("Charge_Tl_208_full", eres_count, eres_count, gain_count, gain_count);    // first MC histogram
            TH1D *mc1_full = MC_Bi_214->ProjectionZ("Charge_Bi_214_full", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
            TH1D *mc2_full = MC_K_40->ProjectionZ("Charge_K_40_full", eres_count, eres_count, gain_count, gain_count);    // second MC histogram

            integrale_gauche_Tl = mc0->Integral(0, lim);
            integrale_droite_Tl = mc0->Integral(lim, 1024);
            integrale_gauche_Bi = mc1->Integral(0, lim);
            integrale_droite_Bi = mc1->Integral(lim, 1024);
            integrale_gauche_K = mc2->Integral(0, lim);
            integrale_droite_K = mc2->Integral(lim, 1024);

            param1 = 0;
            param2 = 0;
            param3 = 0;
            error1 =0;
            error2 =0;
            error3 = 0;

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

            TFractionFitter* fit = new TFractionFitter(spectre_om, mc); // initialise
            fit->Constrain(0, 0, 1);               // constrain fraction 1 to be between 0 and 1
            fit->Constrain(1, 0, 1);
            fit->Constrain(2, 0, 1);

            // fit->SetRangeX(0,307);                // use only the first 15 bins in the fit
            Int_t status = fit->Fit();               // perform the fit

            std::cout << "fit status: " << status << std::endl;

            if (status == 0) {                       // check on fit status

              fit->GetResult(0, param1, error1);
              fit->GetResult(1, param2, error2);
              fit->GetResult(2, param3, error3);

              std::cout << param1 << '\n';

              mc1->Draw();
              mc2->Draw("same");
              mc0->Draw("same");

              Chi2NDF = (fit->GetChisquare())/(fit->GetNDF());


              TH1D* result = (TH1D*) fit->GetPlot();
              TH1D* result_0 = (TH1D*) fit->GetMCPrediction(0);
              TH1D* result_1 = (TH1D*) fit->GetMCPrediction(1);
              TH1D* result_2 = (TH1D*) fit->GetMCPrediction(2);

              gain = (gain_bin_min + gain_bin_width*(gain_count-1));
              eres = eres_bin_min + eres_bin_width*(eres_count-1);
              std::cout << "eres = " << eres << " and gain = " << gain << '\n';

              TCanvas* canvas = new TCanvas;
              canvas->SetLogy();
              spectre_om->Draw("same");
              spectre_om->SetTitle(Form("Fit_simu Eres = %.2f, Gain = %.0f", eres, gain));
              // spectre_om->GetXaxis()->SetRangeUser(0, 40000);
              spectre_om->GetXaxis()->SetTitle("Charge (adc)");
              result_0->Draw("same");
              result_0_scale = param1/result_0->Integral()*spectre_om->Integral();
              result_0->Scale(result_0_scale);
              result_0_scale_er = error1/result_0->Integral()*spectre_om->Integral();
              result_0->SetLineColor(kGreen);
              result_1->Draw("same");
              result_1_scale = param2/result_1->Integral()*spectre_om->Integral();
              result_1_scale_er = error2/result_1->Integral()*spectre_om->Integral();
              result_1->Scale(result_1_scale);
              result_1->SetLineColor(kOrange);
              result_2->Draw("same");
              result_2_scale = param3/result_2->Integral()*spectre_om->Integral();
              result_2_scale_er = error3/result_2->Integral()*spectre_om->Integral();
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
              // return;
              activity_Tl = (mc0->Integral()/integrale_droite_Tl)/1800;
              activity_Bi = (mc1->Integral()/integrale_droite_Bi)/1800;
              activity_K = (mc2->Integral()/integrale_droite_K)/1800;

              error_gain = gain_bin_width;
              error_res = eres_bin_width;
              // return;
              if (eres < 17 && gain < 42000 )
              {
                // canvas->SaveAs(Form("Fit_kolmo/histo_mystere_4_cut/fit_kolmo_eres_%d_gain_%d_lim_%i_non_lin.png", eres_count, gain_count, lim));
                // canvas->SaveAs(Form("Fit_kolmo/histo_mystere/fit_kolmo_eres_%d_gain_%d.png", eres_count, gain_count));
              }
              Chi2->SetBinContent(eres_count, gain_count, Chi2NDF);

              mc0_full->Scale(result_0_scale);
              mc1_full->Scale(result_1_scale);
              mc2_full->Scale(result_2_scale);

              total_hit_Tl = mc0_full->Integral();
              total_hit_Bi = mc1_full->Integral();
              total_hit_K = mc2_full->Integral();




              Result_tree.Fill();
              delete canvas;
              delete fit;
            }
            delete mc0;
            delete mc1;
            delete mc2;
            delete mc0_full;
            delete mc1_full;
            delete mc2_full;

          }
        }
        if (lim >= 25) {
          plus = 20;
        }
      }
    }

 histo_file_Tl->Close();
 histo_file_Bi->Close();
 histo_file_K->Close();


 std::cout << "****************************************************" << '\n';
 std::cout << "The best X2 is : " <<  Result_tree.GetMinimum ("Chi2NDF") << '\n';
 std::cout << "****************************************************" << '\n';

 TH2D* parabole = new TH2D("parabole", "parabole", 500, 6, 21, 1000, 0, 8);
 double good_gain = 0;

 for (int i = 0; i < Result_tree.GetEntries(); i++) {
   Result_tree.GetEntry(i);

   if (Chi2NDF == Result_tree.GetMinimum ("Chi2NDF")) {
     good_gain = gain;
   }
 }

 for (int i = 0; i < Result_tree.GetEntries(); i++) {
   Result_tree.GetEntry(i);

   if (gain == good_gain) {
     parabole->Fill(eres, Chi2NDF);
   }
 }
 newfile->cd();

 parabole->Write();
 Result_tree.Write();
 Chi2->Write();

 newfile->Close();

}

void norm(string name) {

  int eres_count = 0;
  int gain_count = 0;
  double scale0 = 0;
  double scale1 = 0;
  double scale2 = 0;

  TFile* histo_file = new TFile(Form("histo_mystere/test%s_good_Tl.root", name.c_str()), "READ");
  TH3D* MC_Tl_208 = (TH3D*)histo_file->Get("MC_Simu_Tl_208");
  TH3D* MC_Bi_214 = (TH3D*)histo_file->Get("MC_Simu_Bi_214");
  TH3D* MC_K_40 = (TH3D*)histo_file->Get("MC_Simu_K_40");

  TFile *mystere_file = new TFile(Form("histo_mystere/new_histo_mystere_%s.root", name.c_str()), "RECREATE");

  TCanvas* canvas = new TCanvas;
  canvas->SetLogy();

  if (std::stoi(name) == 1){
    eres_count = 5;
    gain_count = 74;
    scale0 = 0.16228;
    scale1 = 0.95660;
    scale2 = 0.75965;
    TH1D* histo_1 = (TH1D*)histo_file->Get("histo_1");
    histo_1->Draw();
  }
  else if (std::stoi(name) == 3){
    eres_count = 99;
    gain_count = 99;
    scale0 = 0.99;
    scale1 = 0.99;
    scale2 = 0.99;
  }

  else {
    eres_count = 21;
    gain_count = 121;
    scale0 = 0.21384;
    scale1 = 0.68612;
    scale2 = 0.27334;
    TH1D* histo_2 = (TH1D*)histo_file->Get("histo_2");
    histo_2->Draw();
    }

  TH1D *mc0 = MC_Tl_208->ProjectionZ("Charge_Tl_208", eres_count, eres_count, gain_count, gain_count);    // first MC histogram
  TH1D *mc1 = MC_Bi_214->ProjectionZ("Charge_Bi_214", eres_count, eres_count, gain_count, gain_count);    // second MC histogram
  TH1D *mc2 = MC_K_40->ProjectionZ("Charge_K_40", eres_count, eres_count, gain_count, gain_count);

  mc0->Scale(scale0);
  mc1->Scale(scale1);
  mc2->Scale(scale2);


  TH1D *mc = MC_Tl_208->ProjectionZ("mc_tot", eres_count, eres_count, gain_count, gain_count);
  mc->Scale(scale0);
  mc->Add(mc1);
  mc->Add(mc2);

  mc0->Draw("same");
  mc0->SetLineColor(kGreen);
  mc1->Draw("same");
  mc1->SetLineColor(kOrange);
  mc2->Draw("same");
  mc2->SetLineColor(kBlue);
  mc->Draw("same");
  mc->SetLineColor(kPink);

  mystere_file->cd();

  mc->Write();

  mystere_file->Close();

}

void fusion() {

  TH1D* Simu_mystere_1 = new TH1D("Simu_mystere_1", "Simu_mystere_1", 1024, 0, 200000);
  TH1D* Simu_mystere_2 = new TH1D("Simu_mystere_2", "Simu_mystere_2", 1024, 0, 200000);
  TH1D* Simu_mystere_3 = new TH1D("Simu_mystere_3", "Simu_mystere_3", 1024, 0, 200000);
  TH1D* Simu_mystere_4 = new TH1D("Simu_mystere_4", "Simu_mystere_4", 1024, 0, 200000);

  TFile *mystere_file_1 = new TFile("histo_mystere/new_histo_mystere_1.root", "READ");
  TH1D* Simu_1 = (TH1D*)mystere_file_1->Get("mc_tot");

  TFile *mystere_file_2 = new TFile("histo_mystere/new_histo_mystere_2.root", "READ");
  TH1D* Simu_2 = (TH1D*)mystere_file_2->Get("mc_tot");

  TFile *mystere_file_3 = new TFile("histo_mystere/new_histo_mystere_3.root", "READ");
  TH1D* Simu_3 = (TH1D*)mystere_file_3->Get("mc_tot");

  TFile *mystere_file_4 = new TFile("histo_mystere/new_histo_mystere_4.root", "READ");
  TH1D* Simu_4 = (TH1D*)mystere_file_4->Get("mc_tot");

  TFile *new_file = new TFile("test.root", "RECREATE");
  new_file->cd();

  Simu_mystere_1->Write();
  Simu_mystere_2->Write();
  Simu_mystere_3->Write();
  Simu_mystere_4->Write();

  new_file->Close();
}

void fit_poly(string name) {

  double min = 0;
  double Chi2NDF = 0;
  float eres = 0;
  float gain = 0;
  int lim_tree = 0;
  double lim = 0;
  double x2 = 0;
  double x1 = 0;
  double x0 = 0;
  double error = 0;


  TFile *file = new TFile(Form("histo_kolmo/Simu_%s.root", name.c_str()), "READ");
  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("Chi2NDF",1);
  tree->SetBranchAddress("Chi2NDF", &Chi2NDF);
  tree->SetBranchStatus("eres",1);
  tree->SetBranchAddress("eres", &eres);
  tree->SetBranchStatus("gain",1);
  tree->SetBranchAddress("gain", &gain);
  tree->SetBranchStatus("lim_tree",1);
  tree->SetBranchAddress("lim_tree", &lim_tree);


  TFile *newfile = new TFile("test_2.root", "RECREATE");
  TTree test_tree("test_tree","");
  test_tree.Branch("lim", &lim);
  test_tree.Branch("min", &min);
  test_tree.Branch("error", &error);


  int plus = 1;
  for (lim = 195; lim < 385; lim+=plus) {

  TF1* f_poly = new TF1 ("f_poly","[p0]*pow(x,2)+[p1]*x+[p2]");
  f_poly->SetParNames("x2","x1","x0");
  f_poly->SetParameters(4, -0.4, 0.02);
  f_poly->SetRange(0,25);
  // gStyle->SetMarkerType(3);
  gStyle->SetMarkerSize(10.);
  tree->Draw("Chi2NDF:(eres) >> map(500,6,21,1000,0,8)",Form("gain > 1.998e-5 && gain < 2e-5 && lim_tree == %f", lim));


  TH2F *map = (TH2F*)gDirectory->Get("map");
  int i = lim;
  f_poly->Draw("same");
  gStyle->SetOptFit(1111);

  map->Fit(f_poly, "RQ");
  min = map->GetFunction("f_poly")->GetMinimumX();

  error = map->GetBinError(map->GetMinimumBin());
  std::cout << error << '\n';

  std::cout << "lim = " << lim << '\n';
  std::cout << "min = " << min  << '\n';

  delete f_poly;
  delete map;
  if (lim>=25) {
    plus = 20;
  }
  test_tree.Fill();
}

newfile->cd();
test_tree.Write();

newfile->Close();

}

void fit_poly2(string name) {

    double min = 0;
    double Chi2NDF = 0;
    float eres = 0;
    float gain = 0;
    int lim_tree = 0;
    double lim = 0;
    double x2 = 0;
    double x1 = 0;
    double x0 = 0;
    double error_plus = 0;
    double error_moins = 0;

    TFile *file = new TFile(Form("histo_kolmo/Simu_%s.root", name.c_str()), "READ");
    TH2D* map = (TH2D*)file->Get("parabole");

    TFile *newfile = new TFile("test_2.root", "RECREATE");
    TTree test_tree("test_tree","");
    test_tree.Branch("lim", &lim);
    test_tree.Branch("min", &min);

    map->Draw();
    // return;
      TF1* f_poly = new TF1 ("f_poly","[p0]*pow(x,2)+[p1]*x+[p2]");
      f_poly->SetParNames("x2","x1","x0");
      f_poly->SetParameters(3, -0.4, 0.02);
      f_poly->SetRange(0,25);

      int i = lim;
      f_poly->Draw("same");
      gStyle->SetOptFit(1111);

      map->Fit(f_poly, "RQ");
      f_poly->SetRange(map->GetFunction("f_poly")->GetMinimumX()-0.5*map->GetFunction("f_poly")->GetMinimumX(),map->GetFunction("f_poly")->GetMinimumX()+0.5*map->GetFunction("f_poly")->GetMinimumX());
      map->Fit(f_poly, "RQ");
      f_poly->SetRange(map->GetFunction("f_poly")->GetMinimumX()-0.5*map->GetFunction("f_poly")->GetMinimumX(),map->GetFunction("f_poly")->GetMinimumX()+0.5*map->GetFunction("f_poly")->GetMinimumX());
      map->Fit(f_poly, "RQ");

      min = map->GetFunction("f_poly")->GetMinimumX();
      std::cout << min << '\n';
      error_plus = f_poly->GetX((map->GetFunction("f_poly")->GetMinimum())+0.00225, map->GetFunction("f_poly")->GetMinimumX(), 200);
      error_moins = f_poly->GetX((map->GetFunction("f_poly")->GetMinimum())+0.00225,-200, map->GetFunction("f_poly")->GetMinimumX());

      std::cout << map->GetFunction("f_poly")->GetMinimum()-0.00225 << '\n';
      if ((error_plus - min) > (min - error_moins-0.001) && (error_plus - min) < (min - error_moins+0.001)) {
        std::cout << "min = "<< min << " ± " << (min - error_moins) << '\n';
      }
      else{
      std::cout << "min = "<< min << " between [" << error_moins << ", " << error_plus << "]" << '\n';
      }



    }

void eff_om_prep(string name) {

  std::vector<int> *om_id = new std::vector<int>;
  std::vector<double> *energy = new std::vector<double>;
  float vertex_position[3];


  TFile *file = new TFile(Form("Histo_simu/Simu_%s.root", name.c_str()), "READ");

  TTree* tree = (TTree*)file->Get("Result_tree");
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("om_id",1);
  tree->SetBranchAddress("om_id", &om_id);
  tree->SetBranchStatus("energy",1);
  tree->SetBranchAddress("energy", &energy);
  tree->SetBranchStatus("vertex_position",1);
  tree->SetBranchAddress("vertex_position", &vertex_position);

  TH1D energy_id("e_id", "e id ", 712, 0, 712);
  TH1D energy_id_cut("e_id_cut", "e id cut ", 712, 0, 712);

  for (int i = 0; i < tree->GetEntries(); i++) {
    memset(vertex_position, 0, 3*sizeof(float));
    tree->GetEntry(i);
    if (i % 100000 == 0) {
      std::cout << "entry = " << i << '\n';
    }

    if (vertex_position[1] > 4999) {
      for (int k = 0; k < om_id->size(); k++) {
        for (int j = 0; j < energy->size(); j++) {
          energy_id.Fill(om_id->at(k)-1);
          if (energy->at(j) > 1) {
            energy_id_cut.Fill(om_id->at(k)-1);
          }
        }
      }
    }
  }

  TFile *newfile = new TFile(Form("eff_om/eff_prep_%s.root", name.c_str()), "RECREATE");
  newfile->cd();
  energy_id.Write();
  energy_id_cut.Write();
  newfile->Close();
}

void eff_om(string name) {

  double lim =0;
  double eff_tot =0;
  double eff_cut = 0;
  double eff_tot_error =0;
  double eff_cut_error = 0;
  std::vector<int> *om_id = new std::vector<int>;
  std::vector<double> *energy = new std::vector<double>;

  TFile *file = new TFile(Form("eff_om/eff_prep_%s.root", name.c_str()), "READ");
  TH1D* eff_om_prep = (TH1D*)file->Get("e_id");
  TH1D* eff_om_prep_cut = (TH1D*)file->Get("e_id_cut");

  TH2F eff_tot_histo("eff", "efficacite mur ", 20, 0, 20, 13, 0, 13);
  TH2F eff_cut_histo("eff_cut", "efficacite mur cut", 20, 0, 20, 13, 0, 13);

  TFile *newfile = new TFile(Form("eff_om/eff_om_%s.root", name.c_str()), "RECREATE");
  TTree new_tree("new_tree","");
  new_tree.Branch("eff_tot", &eff_tot);
  new_tree.Branch("eff_cut", &eff_cut);

  // for (int i = 0; i < eff_om_prep->GetMaximumBin(); i++) {
  for (int i = 0; i < 260; i++) {

    int om_col = (i % 13);
    int om_row = (i / 13);

    eff_tot = (eff_om_prep->GetBinContent(i))/1.0e5;
    eff_tot_error = (sqrt(eff_om_prep->GetBinContent(i))/1.0e5);
    eff_cut = (eff_om_prep_cut->GetBinContent(i))/1.0e5;
    eff_cut_error = (sqrt(eff_om_prep_cut->GetBinContent(i))/1.0e5);

    eff_tot_histo.SetBinContent( om_row+1, om_col+1, eff_tot);
    eff_cut_histo.SetBinContent( om_row+1, om_col+1, eff_cut);

    new_tree.Fill();
  }

file->Close();

newfile->cd();

new_tree.Write();
eff_tot_histo.Write();
eff_cut_histo.Write();

newfile->Close();

}

int main(int argc, char const *argv[])
{
  kolmo_mystere();
  return 0;
}
