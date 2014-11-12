#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
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

class DoubleCB : public RooAbsPdf {
public:
  ClassDef(DoubleCB, 1);
	DoubleCB() {} ;
	DoubleCB(const char* name, const char *title, RooAbsReal &_m,
			RooAbsReal &_m0, RooAbsReal &_sigma,
			RooAbsReal &_alphaHi, RooAbsReal &_nHi,
			RooAbsReal &_alphaLo, RooAbsReal &_nLo);

	DoubleCB(const DoubleCB& other, const char* name = 0);
	virtual TObject* clone(const char* newname) const { return new DoubleCB(*this,newname); }

	virtual Double_t evaluate() const;

private:
	RooRealProxy m;
	RooRealProxy m0;
	RooRealProxy sigma;
	RooRealProxy alphaHi;
	RooRealProxy nHi;
	RooRealProxy alphaLo;
	RooRealProxy nLo;
};

ClassImp(DoubleCB);

DoubleCB::DoubleCB(const char* name, const char *title, RooAbsReal &_m,
		RooAbsReal &_m0, RooAbsReal &_sigma,
		RooAbsReal &_alphaHi, RooAbsReal &_nHi,
		RooAbsReal &_alphaLo, RooAbsReal &_nLo) :
  RooAbsPdf(name, title),
  m("m", "Dependent", this, _m),
  m0("m0", "M0", this, _m0),
  sigma("sigma", "Sigma", this, _sigma),
  alphaHi("alphaHi", "AlphaHi", this, _alphaHi),
  nHi("nHi", "NHi", this, _nHi),
  alphaLo("alphaLo", "AlphaLo", this, _alphaLo),
  nLo("nLo", "NLo", this, _nLo)
{
}

DoubleCB::DoubleCB(const DoubleCB& other, const char* name) :
  RooAbsPdf(other, name),
  m("m", this, other.m),
  m0("m0", this, other.m0),
  sigma("sigma", this, other.sigma),
  alphaHi("alphaHi", this, other.alphaHi),
  nHi("nHi", this, other.nHi),
  alphaLo("alphaLo", this, other.alphaLo),
  nLo("nLo", this, other.nLo)
{
}


Double_t DoubleCB::evaluate() const {
	Double_t t = (m-m0) / sigma;

	if (t < -alphaLo) {
		Double_t a = alphaLo/nLo;
		Double_t b = ( 1 - a * (alphaLo + t) );

		return exp(-0.5*alphaLo*alphaLo) / TMath::Power(b, nLo);
	}
	else if (t > alphaHi) {
		Double_t a = alphaHi/nHi;
		Double_t b = ( 1 - a * (alphaHi + t) );
		
		return exp(-0.5*alphaHi*alphaHi) / TMath::Power(b, nHi);
	}
	else {
		return exp(-0.5*t*t);
	}
}


const float v_unc_lumi = 0.028; // 2.8%


TString fit_withsm( float v_nbg,float bg_slope, float v_xsec_bsm, float v_xsec_sm, float v_unc_xsec_sm, float v_unc_eff)
{

  /////  SIGNAL REGION
  // the measured mass
  RooRealVar mgg("mgg","mgg",125,105,160);

  // signal model params
  RooRealVar mh("mh","mh",125,125,125);
  RooRealVar mbh("mbh","mbh",125,125,125);
  //RooRealVar sh("sh","sh",5,5,5);
  RooRealVar sh("sh","sh",1.39,1.39,1.39);
  RooRealVar alphaCB("alphaCB","alphaCB",1.58,1.58,1.58);
  RooRealVar nCB("nCB", "nCB", 18, 18, 18);
  mh.setConstant(true);
  mbh.setConstant(true);
  sh.setConstant(true);
  alphaCB.setConstant(true);
  nCB.setConstant(true);

  RooRealVar xsec_bsm("xsec_bsm","xsec_bsm",0,0,100);
  RooRealVar xsec_sm("xsec_sm","xsec_sm",0,0,100);

  RooRealVar eff("eff","eff",1,0.5,1.5);
  RooRealVar lumi("lumi","lumi",1,0.5,1.5);

  // bg model params
  RooRealVar nbg("nbg","nbg",0,0,1000);
  RooRealVar bgc("bgc","bgc",bg_slope,2*bg_slope,-2*bg_slope);
  //  bgc.setConstant(true);
  
  RooRealVar alphaHi("alphaHi", "alphaHi", 1.81, 1.81, 1.81);
  RooRealVar nHi("nHi", "nHi", 13.5, 13.5, 13.5);
  RooRealVar alphaLo("alphaLo", "alphaLo", 1.34, 1.34, 1.34);
  RooRealVar nLo("nLo", "nLo", 23.4, 23.4, 23.4);

  alphaHi.setConstant(true);
  nHi.setConstant(true);
  alphaLo.setConstant(true);
  nLo.setConstant(true);
  
  // signal PDF
  //RooGaussian G("G","G",mgg,mbh,sh);
  RooCBShape G("G", "G", mgg, mbh, sh, alphaCB, nCB);

  /*
    std::cout << "initializing doubleCB..." << std::endl;
  DoubleCB G("G", "G", mgg, mbh, sh, alphaHi, nHi, alphaLo, nLo);
  std::cout << "done!." << std::endl;
  std::cout << "val = " << G.evaluate() << std::endl;
  */

  // SM PDF
  //RooGaussian S("S","S",mgg,mh,sh);
  RooCBShape S("S", "S", mgg, mh, sh, alphaCB, nCB);
  /*
  std::cout << "initializing doubleCB(2)..." << std::endl;
  DoubleCB S("S", "S", mgg, mh, sh, alphaHi, nHi, alphaLo, nLo);
  std::cout << "done!." << std::endl;
  */

  // bg PDF
  RooExponential E("E","E",mgg,bgc);
  
  RooWorkspace wspace("wspace");

  wspace.import(G);
  wspace.import(S);
  wspace.import(E);
  wspace.import(xsec_bsm);
  wspace.import(xsec_sm);
  wspace.import(eff);
  wspace.import(lumi);
  wspace.import(nbg);

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
  wspace.var("lumi")->setVal(1.0);

  // unconstrained PDF
  RooAbsPdf *pdf = wspace.pdf("jointModel");

  // constraint PDF
  wspace.factory("Gaussian::theoryconstraint(theoryGlobs[0,100], xsec_sm, unc_xsec_sm[0,100])");
  wspace.factory(   "Gaussian::effconstraint(effGlobs[0,5]     , eff    , unc_eff[0,1.])");
  wspace.factory(  "Gaussian::lumiconstraint(lumiGlobs[0,5]    , lumi   , unc_lumi[0,1.])");
  wspace.var("theoryGlobs")->setVal(v_xsec_sm);
  wspace.var("effGlobs"   )->setVal(1.0);
  wspace.var("lumiGlobs"  )->setVal(1.0);
  wspace.var("unc_xsec_sm")->setVal(v_unc_xsec_sm);
  wspace.var("unc_eff"    )->setVal(v_unc_eff);
  wspace.var("unc_lumi"   )->setVal(v_unc_lumi);
  wspace.var("theoryGlobs")->setConstant(); // like it's data
  wspace.var("effGlobs"   )->setConstant(); // like it's data
  wspace.var("lumiGlobs"  )->setConstant(); // like it's data

  // If unc_* not fixed, treats it as a parameter of fit. as unc->0 gaussian diverges
  wspace.var("unc_xsec_sm")->setConstant();
  wspace.var("unc_eff"    )->setConstant(); // added by kyle. 
  wspace.var("unc_lumi"   )->setConstant(); // copied by danielw

  wspace.factory("PROD::jointModeld(jointModel, theoryconstraint, effconstraint,lumiconstraint)");
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
  RooArgSet cas(xsec_sm,eff,lumi);
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
  wspace.var("lumi"    )->setConstant(true);
  wspace.var("xsec_sm")->setVal(v_xsec_sm);
  wspace.var("eff"    )->setVal(1.0);
  wspace.var("lumi"   )->setVal(1.0);
  TH1* nllHist = profileJoint->createHistogram("xsec_bsm",100);
  wspace.import(*nllHist,"profLLeff");
  wspace.var("xsec_sm")->setConstant(false);
  wspace.var("eff"    )->setConstant(false);
  wspace.var("lumi"   )->setConstant(false);

  RooStats::ModelConfig mc("ModelConfig",&wspace);
  //  mc.SetPdf(*pdf);
  mc.SetPdf(*pdfc);
  //  mc.SetParametersOfInterest(*wspace.var("nsig"));
  mc.SetParametersOfInterest(*wspace.var("xsec_bsm"));
  mc.SetObservables(*wspace.var("mgg"));
  wspace.defineSet("nuisParams","nbg,eff,lumi,xsec_sm");
  
  mc.SetNuisanceParameters(*wspace.set("nuisParams"));
  wspace.import(mc);

  TString fname = Form("monoh_withsm_SRCR_bg%1.1f_bgslop%1.1f_nsig%1.1f.root",v_nbg,bg_slope,v_xsec_bsm);

  wspace.writeToFile(fname,true);
  std::cout << " Written WS to " << fname << std::endl;
  return fname;
 
}

void limit(TString fname)
{
  StandardHypoTestInvDemo( fname, "wspace","ModelConfig","","data",2,3,true,50,0,20);
}


void limit_bands(float v_nbg,float bg_slope, float v_xsec_bsm, float v_xsec_sm, float v_unc_xsec_sm, float v_unc_eff)
{

  TString fname = fit_withsm(v_nbg,bg_slope,v_xsec_bsm,v_xsec_sm,v_unc_xsec_sm,v_unc_eff);
  std::cout << "Reading WS from " << fname << std::endl;
  limit(fname);
}
