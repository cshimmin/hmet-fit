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
#include "StandardHypoTestInvDemo.C"
#include "RooDCB.h"

TString fit_withsm( float v_nbg,float bg_slope, float v_xsec_bsm, float v_xsec_sm, float v_unc_xsec_sm, float v_unc_eff)
{
  // hard coded lumi stuff
  double v_lumi = 20.3;
  double v_unc_lumi = 0.028;//*v_lumi;

  // peak location /width uncertainties
  double v_unc_mh = 0.003;
  double v_unc_sh = 0.11;

  // signal shape params
  
  // SM higgs 
  //  expected limit (median) 0.262023  
  double v_mh = 125.0;
  double v_sh = 1.354;
  double v_alphaHi = 1.736;
  double v_nHi = 156.3;
  double v_alphaLo = 1.344;
  double v_nLo = 20.06;
  
  // BSM average
  //  expected limit (median) 0.260314
  /*
  double v_mh = 125.0;
  double v_sh = 1.14;
  double v_alphaHi = 1.53;
  double v_nHi = 152.95;
  double v_alphaLo = 1.32;
  double v_nLo = 45.48;


  // BSM average for simp models and low DM mass
  //  expected limit (median) 0.264728
  double v_mh = 125.0;
  double v_sh = 1.40;
  double v_alphaHi = 1.50;
  double v_nHi = 152.63;
  double v_alphaLo = 1.35;
  double v_nLo = 32.76;
  */

  /////  SIGNAL REGION
  // the measured mass
  RooRealVar mgg("mgg","mgg",v_mh,105,160);

  RooRealVar xsec_bsm("xsec_bsm","xsec_bsm",0,0,1);
  RooRealVar xsec_sm("xsec_sm","xsec_sm",0,0,100);

  RooRealVar eff("eff","eff",1,0.5,1.5);
  RooRealVar lumi("lumi","lumi",v_lumi,v_lumi-3*v_unc_lumi, v_lumi+3*v_unc_lumi); //KC
  RooRealVar lumiGlobs("lumiGlobs","lumiGlobs",v_lumi,v_lumi-3*v_unc_lumi, v_lumi+3*v_unc_lumi); //KC
  RooRealVar unc_lumi("unc_lumi","unc_lumi",v_unc_lumi,0.1,2*v_unc_lumi); //KC
  RooRealVar unc_mh("unc_mh","unc_mh",v_unc_mh,0.001,2*v_unc_mh); //DW copying KC
  RooRealVar unc_sh("unc_sh","unc_sh",v_unc_sh,0.001,2*v_unc_sh); //DW copying KC
  

  // bg model params
  RooRealVar nbg("nbg","nbg",0,0,1000);
  RooRealVar bgc("bgc","bgc",bg_slope,2*bg_slope,-2*bg_slope);
  //  bgc.setConstant(true);
  

  // signal shape parameters
  RooRealVar mh("mh", "mh", v_mh,120,130);
  RooRealVar sigma_h("sigma_h", "sigma_h", v_sh,v_sh*(1.0-3*v_unc_sh),v_sh*(1.0+3*v_unc_sh));
  RooRealVar alphaHi("alphaHi", "alphaHi", v_alphaHi);
  RooRealVar nHi("nHi", "nHi", v_nHi);
  RooRealVar alphaLo("alphaLo", "alphaLo", v_alphaLo);
  RooRealVar nLo("nLo", "nLo", v_nLo);

  //  mh.setConstant(true);  // no longer true, since this is now a NP
  //  sigma_h.setConstant(true);// no longer true, since this is now a NP
  alphaHi.setConstant(true);
  nHi.setConstant(true);
  alphaLo.setConstant(true);
  nLo.setConstant(true);
  
  // signal PDF
  RooDCB G("G", "G", mgg, mh, sigma_h, alphaHi, nHi, alphaLo, nLo);

  // SM PDF
  RooDCB S("S", "S", mgg, mh, sigma_h, alphaHi, nHi, alphaLo, nLo);

  // bg PDF
  RooExponential E("E","E",mgg,bgc);  // nominal model
  //RooUniform E("E","E",mgg);    // alternate model 
  
  RooWorkspace wspace("wspace");

  wspace.import(G);
  wspace.import(S);
  wspace.import(E);
  wspace.import(xsec_bsm);
  wspace.import(xsec_sm);
  wspace.import(eff);
  wspace.import(nbg);
  wspace.import(lumi);
  wspace.import(lumiGlobs);
  wspace.import(unc_lumi);
  wspace.import(unc_mh);
  wspace.import(unc_sh);

  // number of signal events
  wspace.factory("prod::nsig(lumi,eff,xsec_bsm)");  

  // number of SM higgs events
  wspace.factory("prod::nsm(lumi,eff,xsec_sm)");


  // the joint model: non-peaking BG + peaking BG (called "SM") and peaking signal
  wspace.factory("SUM:jointModel(nsig*G,nsm*S,nbg*E)");
  wspace.var("xsec_bsm")->setVal(v_xsec_bsm);
  wspace.var("nbg")->setVal(v_nbg);
  wspace.var("xsec_sm")->setVal(v_xsec_sm);
  wspace.var("eff")->setVal(1.0);
  wspace.var("lumi")->setVal(v_lumi);

  // unconstrained PDF
  RooAbsPdf *pdf = wspace.pdf("jointModel");

  // constraint PDF
  wspace.factory("Gaussian::theoryconstraint(theoryGlobs[0,100], xsec_sm, unc_xsec_sm[0,100])");
  wspace.factory(   "Gaussian::effconstraint(effGlobs[0,5]     , eff    , unc_eff[0,1.])");
  wspace.factory(  "Gaussian::peakconstraint(peakGlobs[120,130], mh     , unc_mh)");
  wspace.factory(  "Gaussian::widthconstraint(widthGlobs[0,3]  , sigma_h, unc_sh)");
  wspace.factory(  "Gaussian::lumiconstraint(lumiGlobs         , lumi   , unc_lumi)"); //KC: limits defined elsewhere
  wspace.var("theoryGlobs")->setVal(v_xsec_sm);
  wspace.var("effGlobs"   )->setVal(1.0);
  wspace.var("peakGlobs"  )->setVal(v_mh);
  wspace.var("widthGlobs" )->setVal(v_sh);
  wspace.var("lumiGlobs"  )->setVal(1.0);
  wspace.var("unc_xsec_sm")->setVal(v_unc_xsec_sm);
  wspace.var("unc_eff"    )->setVal(v_unc_eff);
  wspace.var("unc_mh"     )->setVal(v_unc_mh*v_mh);
  wspace.var("unc_sh"     )->setVal(v_unc_sh*v_sh);
  wspace.var("unc_lumi"   )->setVal(  v_unc_lumi );

  wspace.var("theoryGlobs")->setConstant(); // like it's data
  wspace.var("effGlobs"   )->setConstant(); // like it's data
  wspace.var("lumiGlobs"  )->setConstant(); // like it's data
  wspace.var("peakGlobs"  )->setConstant(); // like it's data
  wspace.var("widthGlobs" )->setConstant(); // like it's data

  // If unc_* not fixed, treats it as a parameter of fit. as unc->0 gaussian diverges
  wspace.var("unc_xsec_sm")->setConstant();
  wspace.var("unc_eff"    )->setConstant(); // added by kyle. 
  wspace.var("unc_mh"     )->setConstant(); //  copied by danielw
  wspace.var("unc_lumi"   )->setConstant(); // copied by danielw
  wspace.var("unc_sh"     )->setConstant(); // copied by danielw

  wspace.factory("PROD::jointModeld(jointModel, theoryconstraint, widthconstraint, peakconstraint, effconstraint,lumiconstraint)");
  RooAbsPdf *pdfc = wspace.pdf("jointModeld");

  pdfc->graphVizTree("debug.dot");

  /// generate  asimov dataset  
  RooRealVar *x_mgg = wspace.var("mgg");
  x_mgg->setBins(20);

  RooDataSet *data = pdfc->generate( *x_mgg, RooFit::ExpectedData() ); // asimov dataset

  data->SetName("data");
  data->Print();
  wspace.import(*data);
  wspace.Print();

  // try to plot it
  RooPlot *plot1 = x_mgg->frame();
  data->plotOn(plot1);

  // fit to the model
  RooArgSet cas(xsec_sm,eff,mh,sigma_h,lumi);
  RooFitResult *r = pdfc->fitTo(*data,RooFit::Constrain(cas),RooFit::Save(true));
  r->Print();

  //  pdf->plotOn(plot1,RooFit::ProjWData(*data));
  pdfc->plotOn(plot1);
  pdfc->plotOn(plot1, RooFit::Components("G"), RooFit::LineColor(kRed));
  pdfc->plotOn(plot1, RooFit::Components("S"), RooFit::LineColor(kGreen));
  pdfc->plotOn(plot1);
  pdfc->paramOn(plot1);

  TCanvas *tc = new TCanvas("tc","",400,400);
  plot1->Draw();


  TCanvas* c1 = new TCanvas("c1","",400,400);
  RooPlot *frame = wspace.var("xsec_bsm")->frame();
  RooAbsReal *nllJoint = pdfc->createNLL(*data, RooFit::Constrained()); // slice with fixed xsec_bsm
  RooAbsReal *profileJoint = nllJoint->createProfile(*wspace.var("xsec_bsm"));
  nllJoint->plotOn(frame, RooFit::LineColor(kRed), RooFit::ShiftToZero());
  profileJoint->plotOn(frame);
  frame->Draw();

  wspace.var("xsec_bsm")->setConstant(true);
  wspace.var("eff"     )->setConstant(true);
  wspace.var("mh"      )->setConstant(true);
  wspace.var("sigma_h" )->setConstant(true);
  wspace.var("lumi"    )->setConstant(true);
  wspace.var("xsec_sm" )->setVal(v_xsec_sm);
  wspace.var("eff"     )->setVal(1.0);
  wspace.var("lumi"    )->setVal(1.0);
  TH1* nllHist = profileJoint->createHistogram("xsec_bsm",100);
  wspace.import(*nllHist,"profLLeff");
  wspace.var("xsec_sm")->setConstant(false);
  wspace.var("eff"    )->setConstant(false);
  wspace.var("lumi"   )->setConstant(false);
  wspace.var("mh"     )->setConstant(false);
  wspace.var("sigma_h")->setConstant(false);

  RooStats::ModelConfig mc("ModelConfig",&wspace);
  //  mc.SetPdf(*pdf);
  mc.SetPdf(*pdfc);
  //  mc.SetParametersOfInterest(*wspace.var("nsig"));
  mc.SetParametersOfInterest(*wspace.var("xsec_bsm"));
  mc.SetObservables(*wspace.var("mgg"));
  wspace.defineSet("nuisParams","nbg,eff,mh,sigma_h,lumi,xsec_sm");
  
  mc.SetNuisanceParameters(*wspace.set("nuisParams"));
  wspace.import(mc);

  TString fname = Form("monoh_withsm_SRCR_bg%1.1f_bgslop%1.1f_nsig%1.1f.root",v_nbg,bg_slope,v_xsec_bsm);

  wspace.writeToFile(fname,true);
  std::cout << " Written WS to " << fname << std::endl;
  return fname;
 
}

void limit(TString fname)
{
  StandardHypoTestInvDemo( fname, "wspace","ModelConfig","","data",2,3,true,50,0,100);
}


void limit_bands(float v_nbg,float bg_slope, float v_xsec_bsm, float v_xsec_sm, float v_unc_xsec_sm, float v_unc_eff)
{

  TString fname = fit_withsm(v_nbg,bg_slope,v_xsec_bsm,v_xsec_sm,v_unc_xsec_sm,v_unc_eff);
  std::cout << "Reading WS from " << fname << std::endl;
  limit(fname);
}

void test(){
  limit_bands(11.67,-1.0/60.176,0,1.13/20.3,0.07,0.05);
}

/*
 Exponential:
Expected upper limits, using the B (alternate) model : 
 expected limit (median) 0.262295
 expected limit (-1 sig) 0.176799
 expected limit (+1 sig) 0.401923
 expected limit (-2 sig) 0.125894
 expected limit (+2 sig) 0.604977

 Uniform:
Expected upper limits, using the B (alternate) model : 
 expected limit (median) 0.258491
 expected limit (-1 sig) 0.174101
 expected limit (+1 sig) 0.395969
 expected limit (-2 sig) 0.123474
 expected limit (+2 sig) 0.596594


*/
