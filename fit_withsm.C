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

class DoubleCB : public RooAbsPdf {
public:
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

//ClassImp(DoubleCB)
//;

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


TString fit_withsm( float v_nbg,float bg_slope, float v_nsig, float v_nsm, float sig_nsm, float sig_eff)
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
  RooRealVar nsig("nsig","nsig",0,0,100);
  RooRealVar nsm("nsm","nsm",0,0,100);

  RooRealVar eff("eff","eff",1,0.5,1.5);

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
  //RooCBShape G("G", "G", mgg, mbh, sh, alphaCB, nCB);
  std::cout << "initializing doubleCB..." << std::endl;
  DoubleCB G("G", "G", mgg, mbh, sh, alphaHi, nHi, alphaLo, nLo);
  std::cout << "done!." << std::endl;
  std::cout << "val = " << G.evaluate() << std::endl;

  // SM PDF
  //RooGaussian S("S","S",mgg,mh,sh);
  //RooCBShape S("S", "S", mgg, mh, sh, alphaCB, nCB);
  std::cout << "initializing doubleCB(2)..." << std::endl;
  DoubleCB S("S", "S", mgg, mh, sh, alphaHi, nHi, alphaLo, nLo);
  std::cout << "done!." << std::endl;

  // bg PDF
  RooExponential E("E","E",mgg,bgc);
  
  RooWorkspace wspace("wspace");

  wspace.import(G);
  wspace.import(S);
  wspace.import(E);
  wspace.import(nsig);
  wspace.import(nsm);
  wspace.import(eff);
  wspace.import(nbg);
  wspace.factory("prod::tempsig(eff,nsig)");
  wspace.factory("prod::tempsm(eff,nsm)");
  wspace.factory("SUM:jointModel(tempsig*G,tempsm*S,nbg*E)");
  wspace.var("nsig")->setVal(v_nsig);
  wspace.var("nbg")->setVal(v_nbg);
  wspace.var("nsm")->setVal(v_nsm);
  wspace.var("eff")->setVal(1.0);

  // unconstrained PDF
  RooAbsPdf *pdf = wspace.pdf("jointModel");

  // constraint PDF
  wspace.factory("Gaussian::smconstraint(nsmGlobs[0,100], nsm, sigma[0,100])");
  wspace.factory("Gaussian::effconstraint(effGlobs[0,5], eff, sigmaeff[0,1.])");
  wspace.var("nsmGlobs")->setVal(v_nsm);
  wspace.var("effGlobs")->setVal(1.0);
  wspace.var("sigma")->setVal(sig_nsm);
  wspace.var("sigmaeff")->setVal(sig_eff);
  wspace.var("nsmGlobs")->setConstant(); // like it's data
  wspace.var("effGlobs")->setConstant(); // like it's data
  wspace.var("sigma")->setConstant();
  wspace.var("sigmaeff")->setConstant(); // added by kyle. 
  // If sigmaeff not fixed, treats it as a parameter of fit. as sigma->0 gaussian diverges

  wspace.factory("PROD::jointModelb(jointModel, smconstraint)");
  wspace.factory("PROD::jointModelc(jointModelb, effconstraint)");
  wspace.factory("PROD::jointModeld(jointModel, smconstraint, effconstraint)");
  RooAbsPdf *pdfc = wspace.pdf("jointModeld");

  pdfc->graphVizTree("debug.dot");

  /// generate  asimov dataset  
  RooRealVar *x_mgg = wspace.var("mgg");
  x_mgg->setBins(20);

  RooDataSet *data = pdfc->generate( *x_mgg, RooFit::ExpectedData() ); // asimov dataset
//  RooDataSet *data = pdfc->generate( *x_mgg)); // random data set
  data->SetName("data");
  data->Print();
  wspace.import(*data);
  wspace.Print();

  // try to plot it
  RooPlot *plot1 = x_mgg->frame();
  data->plotOn(plot1);

  // fit to the model
  RooArgSet cas(nsm,eff);
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
  RooPlot *frame = wspace.var("nsig")->frame();
//  nllJoint = pdfc->createNLL(*data, RooFit::Constrain(nsm)); // slice with fixed nsm
//  nllJoint = pdfc->createNLL(*data, RooFit::Constrain(*wspace.var("nsm"))); // slice with fixed nsm
  RooAbsReal *nllJoint = pdfc->createNLL(*data, RooFit::Constrained()); // slice with fixed nsm
  RooAbsReal *profileJoint = nllJoint->createProfile(*wspace.var("nsig"));
  nllJoint->plotOn(frame, RooFit::LineColor(kRed), RooFit::ShiftToZero());
  profileJoint->plotOn(frame);
  frame->Draw();

 wspace.var("nsm")->setConstant(true);
 wspace.var("eff")->setConstant(true);
  wspace.var("nsm")->setVal(v_nsm);
  wspace.var("eff")->setVal(1.0);
  TH1* nllHist = profileJoint->createHistogram("nsig",100);
  wspace.import(*nllHist,"profLLeff");
  wspace.var("nsm")->setConstant(false);
  wspace.var("eff")->setConstant(false);

  /*
  wspace.defineSet("two","nsig,nsm");
  RooAbsReal* profileContour = nllJoint->createProfile(*wspace.set("two"));
  RooRealVar* pnsig = wspace.var("nsig");

    pnsig->setMax(pnsig->getVal()+2*pnsig->getError());
  //  pnsig->setMin(pnsig->getVal()-2*pnsig->getErrorLo());
  RooRealVar* pnsm = wspace.var("nsm");
  pnsm->setMax(pnsm->getVal()+2*pnsm->getError());
  if(pnsm->getVal()-2*pnsm->getError() >0)
    pnsm->setMin(pnsm->getVal()-2*pnsm->getError());
  wspace.set("two")->Print("v");


  TH1* hist = profileContour->createHistogram("nsig,nsm", 30,30);
//  hist->SetMaximum(10.);
  hist->SetMinimum(0.);
  hist->Draw("colz");
  wspace.set("two")->Print("v");
  */

  RooStats::ModelConfig mc("ModelConfig",&wspace);
  //  mc.SetPdf(*pdf);
  mc.SetPdf(*pdfc);
  mc.SetParametersOfInterest(*wspace.var("nsig"));
  mc.SetObservables(*wspace.var("mgg"));
  wspace.defineSet("nuisParams","nbg,eff,nsm");
  
  mc.SetNuisanceParameters(*wspace.set("nuisParams"));
  wspace.import(mc);

  TString fname = Form("monoh_withsm_SRCR_bg%1.1f_bgslop%1.1f_nsig%1.1f.root",v_nbg,bg_slope,v_nsig);

  wspace.writeToFile(fname,true);
  std::cout << " Written WS to " << fname << std::endl;
  return fname;
 
}
