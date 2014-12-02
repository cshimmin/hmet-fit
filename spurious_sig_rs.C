#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TMath.h"
#include "TText.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TLegend.h"

#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooUniform.h"
#include "RooBernstein.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooDCB.h"
#include "TROOT.h"

/*

Script to calculate the "spurious signal" for mono-Higgs, h->gg.
1) Select sideband region in data, either "PT","MET" or "NOTIGHT"
2) Fit bg param function to the data, for a set of mH values
   bg param functions considered: exponential, flat, bernstein1, bernstein2
3) extract N_sig and uncertainty on N_bkg
4) compare N_sig to 0.2 * N_bkg

 */

// some variables for tree-> Roo conversion
  RooRealVar diphoton_m("diphoton_m","diphoton_m",105e3,160e3);
  RooRealVar met_et("met_et","met_et",0,200e3);
  RooRealVar both_tight("both_tight","both_tight",0,1);
  RooRealVar diphoton_pt("diphoton_pt","diphoton_pt",0,200e3);

// select side band region
#define SB_MET

#ifdef SB_NOTIGHT
char cuts[100]="met_et>50e3&&diphoton_pt>50e3&&both_tight<1&&diphoton_m>105e3&&diphoton_m<160e3";
char tag[100]="notight";
char fname[500] = "../param/data_withtight.root";
#endif

#ifdef SB_MET
char cuts[100]="met_et>30e3&&met_et<90e3&&diphoton_pt>80e3&&diphoton_m>105e3&&diphoton_m<160e3";
char tag[100]="met";
char fname[500] = "../param/data_unblinded.root";
#endif

#ifdef SB_PT
char cuts[100]="met_et>80e3&&diphoton_pt>30e3&&diphoton_pt<90e3&&diphoton_m>105e3&&diphoton_m<160e3";
char tag[100]="pt";
char fname[500] = "../param/data_unblinded.root";
#endif

// singleton class to hold the tree->Roo data conversion
class RootData
{
private:
  RootData()
  {
    std::cout << " Importing data Tree -> Roo" << std::endl;
    // import the data into a Roo DS
  TFile *tf = new TFile(fname);
  TNtuple *nt = (TNtuple*)tf->Get("hmet");


    ds = new RooDataSet("ds","ds",RooArgSet(met_et,diphoton_pt,diphoton_m,both_tight),RooFit::Import(*nt),RooFit::Cut(cuts));
    ds->Print();
    std::cout << " Done importing data Tree -> Roo: " << ds<<std::endl;
  }

  ~RootData() { }
public:
 
  static RooDataSet *ds;
  static RootData *_instance;

  static RootData *getInstance()
  {
    if (!_instance)
      {
	_instance = new RootData();
      }
    return _instance;
  }
 

};

RooDataSet *RootData::ds = 0;
RootData *RootData::_instance=0;

// fit the data (defined in global "fname" var) to Signal + background, where background is  RooAbsPdf &B
// mh_in :  assumed higgs mass
// nsig_out : measured n_sig
// bg_unc_out : measured uncertainty on bg normalization
void fit_rs(float mh_in,float &nsig_out, float &bg_unc_out, RooAbsPdf &B)
{


 // import the data into a Roo DS
  
  RootData *rd = RootData::getInstance();
  RooDataSet *ds = rd->ds;
  std::cout << "Got instance... " << ds << std::endl;
  ds->Print();

  ///// Signal PDF


  // shape parameters
  RooRealVar mh("mh", "mh", mh_in*1e3,mh_in*1e3,mh_in*1e3);
  RooRealVar mbh("mbh", "mbh", mh_in*1e3,mh_in*1e3,mh_in*1e3);
  RooRealVar sigma_h("sigma_h", "sigma_h", 1.354e3);
  RooRealVar alphaHi("alphaHi", "alphaHi", 1.736e3);
  RooRealVar nHi("nHi", "nHi", 156.3);
  RooRealVar alphaLo("alphaLo", "alphaLo", 1.344e3);
  RooRealVar nLo("nLo", "nLo", 20.06);

  mh.setConstant(true);
  mbh.setConstant(true);
  sigma_h.setConstant(true);
  alphaHi.setConstant(true);
  nHi.setConstant(true);
  alphaLo.setConstant(true);
  nLo.setConstant(true);
  
  // signal PDF
  RooDCB G("G", "G", diphoton_m, mbh, sigma_h, alphaHi, nHi, alphaLo, nLo);


  // Signal PDF: Gaussian for now, replace with DCB 
  /*
  RooRealVar mh("mh","mh",mh_in*1e3,mh_in*1e3,mh_in*1e3);
  RooRealVar sh("sh","sh",2e3,2e3,2e3);
  sh.setConstant(true);
  mh.setConstant(true);
  RooGaussian G("G","G",diphoton_m,mh,sh);
  */

  // scales
  RooRealVar nsig("nsig","number of signal events",1,0.,10000) ;
  RooRealVar nbkg("nbkg","number of background events",10,0,10000) ;

  // final model
  RooAddPdf  model("model","(g1+g2)+a",RooArgList(B,G),RooArgList(nbkg,nsig)) ;

  // fit the PDF
  RooFitResult *fitres = (RooFitResult*)model.fitTo(*ds,RooFit::Save(true));
 
  // Plot
  RooPlot* xframe = diphoton_m.frame(RooFit::Title("fit title"));
  ds->plotOn(xframe);
  model.plotOn(xframe,RooFit::Normalization(1.0));
  
  model.plotOn(xframe,RooFit::Components(B),RooFit::LineStyle(kDashed),RooFit::Normalization(1.0));
  model.plotOn(xframe,RooFit::Components(G),RooFit::LineStyle(kDotted),RooFit::Normalization(1.0));

  TCanvas *tc = new TCanvas(Form("%s_%1.1f",B.GetName(),mh_in),Form("%s_%1.1f",B.GetName(),mh_in),400,400);
  xframe->Draw();
  
  tc->Print(Form("spurious_%s_%s_%1.1f.pdf",tag,B.GetName(),mh_in));

  // get the result out
  RooArgList ral_fitresult = fitres->floatParsFinal();
  RooRealVar *nbkg_fitresult = (RooRealVar*)ral_fitresult.find("nbkg");
  if (nbkg_fitresult)
    {
      B.Print();
    std::cout << " nsig = " << nsig.getVal() << " nbg = " << nbkg.getVal() << " +- " << nbkg_fitresult->getError() << std::endl;
    }
  // put it in the output vars
  nsig_out = nsig.getVal();
  bg_unc_out = nbkg_fitresult->getError();


}


// index = where in the TGraph arrays to write the results
// arrays for nsignal fitted, and +  and - bg uncertainty bounds
// RooAbsPdf &B is the bg parameterization

void spurious_signal(int index, TGraph *tgs[],TGraph *tgb[],TGraph*tgbd[], RooAbsPdf &B)
{

  int np=0;
  float mhs[1000];
  float nsig[1000];
  float bgu[1000];
  float bgud[1000];

  // loop over higgs masses
  for (float mh=115;mh<135;mh += 1)
    {
      mhs[np] = mh;

      // fit the data to the bg parameter + gaussian at mh, getting out
      //    nsig[np] = number of signal events fit
      //    bgu[np] = uncertainty on background
      fit_rs(mh,nsig[np],bgu[np],B);

      std::cout << " mh = " << mh << " nsig = " << nsig[np] << " unc = " << bgu[np] << std::endl;
      // 20% criterion
      bgu[np] *= 0.2;

      // make nevative unvelope
      bgud[np] = -bgu[np];
    
      np++;
    }

  std::cout << " making graphs " << std::endl;
  //  TCanvas *tc = new TCanvas(funcname,funcname,400,400);
  tgs[index] = new TGraph(np,mhs,nsig);
  tgb[index] = new TGraph(np,mhs,bgu);
  tgbd[index] = new TGraph(np,mhs,bgud);
  
  tgs[index]->SetLineWidth(2);
  tgb[index]->SetLineStyle(2);
  tgbd[index]->SetLineStyle(2);


}

// main routine, runs all the tests
void spurios_sig()
{

  gROOT->ProcessLine(".L RooDCB.C+");
  // number of functions to consider for bg
  const int nf=6;

  // names of functions, for labelling
  char names[nf][100] = { "flat","poly1","poly2","exp","bern1","bern2"};

  // graphs of number of signal events fit, upper edge of envelope (20% of bg uncertainty), and lower edge
  TGraph *tgs[nf];
  TGraph *tgb[nf];
  TGraph *tgbd[nf];

  // each call below tests one function, putting the result into the arrays tgX

  // 
  RooUniform U("U","U",diphoton_m);
  spurious_signal(0,tgs,tgb,tgbd,U);

  RooRealVar pcoef1a("pcoef1a","pcoef1a",0.5,-10,10);
  RooPolynomial P1("P1","P1",diphoton_m,RooArgList(pcoef1a));
  spurious_signal(1,tgs,tgb,tgbd,P1);

  RooRealVar pcoef1("pcoef1","pcoef1",0.5,-10,10);
  RooRealVar pcoef2("pcoef2","pcoef2",0.5,-10,10);
  RooPolynomial P2("P2","P2",diphoton_m,RooArgList(pcoef1,pcoef2));
  spurious_signal(2,tgs,tgb,tgbd,P2);

  //
  float bg_slope = -1/(60e3);
  RooRealVar bgc("bgc","bgc",bg_slope,2*bg_slope,-2*bg_slope);
  RooExponential E("E","E",diphoton_m,bgc);
  spurious_signal(3,tgs,tgb,tgbd,E);

  RooRealVar coef1("coef1","coef1",0.5,-10,10);
  RooRealVar coef2("coef2","coef2",0.5,-10,10);

  RooBernstein B2("B2","B2",diphoton_m,RooArgList(coef1,coef2));
  spurious_signal(4,tgs,tgb,tgbd,B2);

  RooRealVar bcoef1("bcoef1","bcoef1",0.5,-10,10);
  RooRealVar bcoef2("bcoef2","bcoef2",0.5,-10,10);
  RooRealVar bcoef3("bcoef3","bcoef3",0.5,-10,10);
  RooBernstein B3("B3","B3",diphoton_m,RooArgList(bcoef1,bcoef2,bcoef3));
  spurious_signal(5,tgs,tgb,tgbd,B3);





  
  
  // summarize the results
  TCanvas *tc = new TCanvas(tag,tag,400,400);

  tgs[0]->SetMaximum(40);
  tgs[0]->SetMinimum(-10);

  tgs[0]->SetTitle(Form("%s control region",tag));
  tgs[0]->Draw("ALP");
  tgs[0]->GetXaxis()->SetTitle("m_{H} [GeV]");
  tgs[0]->GetYaxis()->SetTitle("Spurious Signal");

  TLegend *tl = new TLegend(0.6,0.6,0.87,0.87);
  tl->SetBorderSize(0);
  tl->SetFillColor(10);
  tl->SetTextFont(42);
  int cols[nf] = { 1,2,4,6,8,9};

  for (int i=0;i<nf;i++)
    {
      tl->AddEntry(tgs[i],names[i],"L");
      tgs[i]->SetLineColor(cols[i]);
      tgb[i]->SetLineColor(cols[i]);
      tgbd[i]->SetLineColor(cols[i]);
      tgs[i]->Draw("LP");
      tgb[i]->Draw("LP");
      tgbd[i]->Draw("LP");
    }
  tl->Draw();
  tc->Print(Form("spurious_%s.pdf",tag)); 

}


void test_fit()
{
  float nsig, bgu;

  float bg_slope = -1/(60e3);
  RooRealVar bgc("bgc","bgc",bg_slope,2*bg_slope,-2*bg_slope);
  RooExponential E("E","E",diphoton_m,bgc);
  fit_rs(125.0,nsig,bgu,E);
}

