#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooProdPdf.h"
#include "RooAbsPdf.h"
#include "RooExponential.h"
#include "RooUniform.h"
#include "RooFitResult.h"
#include "RooConstVar.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "RooGlobalFunc.h"
#include "RooFit.h"
#include "RooCmdArg.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include "TChain.h"

// USAGE NOTE
//Expects signal files to be in ./signalfiles/



//#include "StandardHypoTestInvDemo.C"
#include "RooDCB.h"

// TODO: need to ask CO where he got this higher-precision lumi
#define LUMI 20.2769
#define FRAC_LUMI_UNC 0.028

// set to true if you want to use real (unblinded) data
// from data_ntuple.root
#define USE_DATA false
#define DATA_FILE "signalfiles/signal_189182.root"

// number of points to scan in StandardHypoTestInvDemo
#define N_POINTS 50

// struct to handle all the arguments to fitting function
struct fitArgs {
  float v_nbg;
  float bg_slope;
  float v_xsec_bsm;
  float v_xsec_sm;
  float in_unc_theory;
  float in_unc_eff;
  float in_unc_sig_theory;
  float v_reco;
};

fitArgs DEFAULT_FIT_ARGS = {
  11.67, -1.0/60.176, 0, 1.069/LUMI, 0.043, 0.029, 0.10, 0.63
};


double check_peak_shfit(fitArgs args, const char* datafile = DATA_FILE){

 // hard coded lumi stuff
  double v_lumi = LUMI;
  double v_unc_lumi = FRAC_LUMI_UNC*v_lumi; //CHECK, make these be absolute uncert, not rel

  double v_eff = 1.; // for symmetry
  double v_unc_eff = args.in_unc_eff*v_eff; // for symmetry

  double v_theory = 1.;
  double v_unc_theory = args.in_unc_theory*v_theory;
  double v_unc_sig_theory = args.in_unc_sig_theory*v_theory;

  // peak location /width uncertainties
  double v_mh = 125.0*1000;
  double v_unc_mh = 0.003*v_mh; //CHECK, make these be absolute uncert, not rel

  double v_sh = 1.354*1000;
  double v_unc_sh = 0.11*v_sh;  //CHECK, make these be absolute uncert, not rel

 
  // SM higgs 
  /*
  double v_alphaHi = 1.736;
  double v_nHi = 156.3;
  double v_alphaLo = 1.344;
  double v_nLo = 20.06;
  */

  double v_alphaHi = 1.736;
  double v_nHi = 156.3;
  double v_alphaLo = 1.344;
  double v_nLo = 20.06;


  RooWorkspace* wspace = new RooWorkspace("wspace");

  // the measured mass
  RooRealVar mgg("diphoton_m","diphoton_m",v_mh,105000,160000);
  RooRealVar diphoton_pt("diphoton_pt", "diphoton_pt", 0, 1000000);
  RooRealVar met_et("met_et", "met_et", 0, 1000000);
  RooRealVar fail_bch("fail_bch", "fail_bch", -1, 2);

  // signal shape parameters
  RooRealVar alphaHi("alphaHi", "alphaHi", v_alphaHi);
  RooRealVar nHi("nHi", "nHi", v_nHi);
  RooRealVar alphaLo("alphaLo", "alphaLo", v_alphaLo);
  RooRealVar nLo("nLo", "nLo", v_nLo);

  RooRealVar mh("mh", "mh", v_mh,120000,130000);
  RooRealVar sigma_h("sigma_h", "sigma_h", v_sh,v_sh-3*v_unc_sh,v_sh+3*v_unc_sh); // could change name to be consistent
  
  //  mh.setConstant(true);  // no longer true, since this is now a NP
  //  sigma_h.setConstant(true);// no longer true, since this is now a NP
  alphaHi.setConstant(true);
  nHi.setConstant(true);
  alphaLo.setConstant(true);
  nLo.setConstant(true);
  
  // signal PDF
  RooDCB G("G", "G", mgg, mh, sigma_h, alphaHi, nHi, alphaLo, nLo);
  wspace->import(G);

  TChain *data_ntuple = new TChain("hmet");
  data_ntuple->Add(Form("signalfiles/%s",datafile));

  RooDataSet*  data = new RooDataSet("data", "data", RooArgSet(mgg, diphoton_pt, met_et, fail_bch), RooFit::Import(*data_ntuple), 
         RooFit::Cut("diphoton_m>105000 && diphoton_m<160000 && met_et>90000 && diphoton_pt>90000 && fail_bch<0.5"));

  data->Print();

  RooRealVar* diphoton_m = wspace->var("diphoton_m");
  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  RooPlot* frame = diphoton_m->frame();
  data->plotOn(frame);
  wspace->pdf("G")->fitTo(*data);
  wspace->pdf("G")->plotOn(frame);
  wspace->pdf("G")->paramOn(frame);
  frame->Draw();
  c1->SaveAs(Form("signalfiles/%s.pdf",datafile));
  return wspace->var("mh")->getVal();
}

void doAll(){
  vector<string> files;

files.push_back("signal_189160.root");
files.push_back("signal_189161.root");
files.push_back("signal_189162.root");
files.push_back("signal_189163.root");
files.push_back("signal_189164.root");
files.push_back("signal_189165.root");
files.push_back("signal_189166.root");
files.push_back("signal_189167.root");
files.push_back("signal_189168.root");
files.push_back("signal_189169.root");
files.push_back("signal_189170.root");
files.push_back("signal_189171.root");
files.push_back("signal_189172.root");
files.push_back("signal_189173.root");
files.push_back("signal_189174.root");
files.push_back("signal_189175.root");
files.push_back("signal_189176.root");
files.push_back("signal_189177.root");
files.push_back("signal_189178.root");
files.push_back("signal_189179.root");
files.push_back("signal_189180.root");
files.push_back("signal_189181.root");
files.push_back("signal_189182.root");
files.push_back("signal_189183.root");
files.push_back("signal_189184.root");
files.push_back("signal_189185.root");
files.push_back("signal_189186.root");
files.push_back("signal_189187.root");
files.push_back("signal_189188.root");
files.push_back("signal_189189.root");
files.push_back("signal_189190.root");
files.push_back("signal_189191.root");
files.push_back("signal_189192.root");
files.push_back("signal_189193.root");
files.push_back("signal_189194.root");
files.push_back("signal_189195.root");
files.push_back("signal_189196.root");
files.push_back("signal_189197.root");
files.push_back("signal_189198.root");
files.push_back("signal_189199.root");

  vector<string>::const_iterator itr = files.begin();
  string results;
  for (; itr!=files.end(); ++itr){
    double mh = check_peak_shfit(DEFAULT_FIT_ARGS, itr->c_str());
    results += *itr+" "+std::to_string(mh) +"\n";
    cout << "XXX: " << itr->c_str() << " " << mh << endl;
  }
  cout << results << endl;
}

