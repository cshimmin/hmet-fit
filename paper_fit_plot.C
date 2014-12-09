#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
#include "TCanvas.h"
#include "TAxis.h"


#include "TLegend.h"
#include "TFrame.h"
#include "TPaveText.h"

//////////////////
//////////////////
//////////////////
/// Code borrowed from https://svnweb.cern.ch/trac/atlasusr/browser/yanght/publications_plotting_macro_2014

TPaveText* CreatePaveText(double x1,double y1,double x2,double y2,vector<TString> text,double textsize){
	  TPaveText *tex=new TPaveText();
	  tex->SetFillColor(0);tex->SetTextSize(0.05);
	  tex->SetFillStyle(0);tex->SetBorderSize(0);
	  int n=text.size();
	  // for(int i=n-1;i>=0;i--) tex->AddText(text[i].Data());
	  for(int i=0;i<n;i++) tex->AddText(text[i].Data());
	  tex->SetX1NDC(x1);
	  tex->SetY1NDC(y1);
	  tex->SetX2NDC(x2);
	  tex->SetY2NDC(y2);
	  tex->SetTextSize(textsize);
	  return tex;
	}
	

void GetX1Y1X2Y2(TVirtualPad *c,double &x1,double &y1,double &x2,double &y2)
{
  x1=c->GetFrame()->GetX1()+c->GetLeftMargin();
  y1=c->GetFrame()->GetY1()+c->GetBottomMargin();
  x2=c->GetFrame()->GetX2()-c->GetRightMargin();
  y2=c->GetFrame()->GetY2()-c->GetTopMargin();
}


TLegend* FastLegend(double x_low, double y_low, double x_up, double y_up, double textsize)
{
  TLegend *legend = new TLegend(x_low,y_low,x_up,y_up);
  legend->SetBorderSize(0);
  legend->SetTextFont(42);
  legend->SetTextSize(textsize);
  legend->SetFillStyle(0);
  legend->SetFillColor(0);
  legend->SetLineColor(0);
  return legend;
 }
 
TStyle* AtlasStyle() 
  {
  TStyle *atlasStyle = new TStyle("ATLAS","Atlas style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  atlasStyle->SetFrameBorderMode(icol);
  atlasStyle->SetFrameFillColor(icol);
  atlasStyle->SetCanvasBorderMode(icol);
  atlasStyle->SetCanvasColor(icol);
  atlasStyle->SetPadBorderMode(icol);
  atlasStyle->SetPadColor(icol);
  atlasStyle->SetStatColor(icol);
  //atlasStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  atlasStyle->SetPaperSize(20,26);

  // set margin sizes
  atlasStyle->SetPadTopMargin(0.05);
  atlasStyle->SetPadRightMargin(0.05);
  atlasStyle->SetPadBottomMargin(0.16);
  atlasStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  atlasStyle->SetTitleXOffset(1.1);
  atlasStyle->SetTitleYOffset(1.3);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.05; // originally 0.05
  atlasStyle->SetTextFont(font);

  atlasStyle->SetTextSize(tsize);
  atlasStyle->SetLabelFont(font,"x");
  atlasStyle->SetTitleFont(font,"x");
  atlasStyle->SetLabelFont(font,"y");
  atlasStyle->SetTitleFont(font,"y");
  atlasStyle->SetLabelFont(font,"z");
  atlasStyle->SetTitleFont(font,"z");
  
  atlasStyle->SetLabelSize(tsize,"x");
  atlasStyle->SetTitleSize(tsize,"x");
  atlasStyle->SetLabelSize(tsize,"y");
  atlasStyle->SetTitleSize(tsize,"y");
  atlasStyle->SetLabelSize(tsize,"z");
  atlasStyle->SetTitleSize(tsize,"z");

  // use bold lines and markers
  atlasStyle->SetMarkerStyle(20);
  atlasStyle->SetMarkerSize(1.2);
  atlasStyle->SetHistLineWidth((Width_t)3.0);
  atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  //atlasStyle->SetErrorX(0.001);
  // get rid of error bar caps
  atlasStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  atlasStyle->SetOptTitle(0);
  //atlasStyle->SetOptStat(1111);
  atlasStyle->SetOptStat(0);
  //atlasStyle->SetOptFit(1111);
  atlasStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  atlasStyle->SetPadTickX(1);
  atlasStyle->SetPadTickY(1);

  return atlasStyle;

}

void SetAtlasStyle()
{
  std::cout << "\nApplying ATLAS style settings...\n" << std::endl ;
  TStyle* atlasStyle = AtlasStyle();
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
}

#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooPlot.h"

//////////////////////////////////////////////////////
//////////////////////////////////////////////////////


void paper_fit_plot()
{
  SetAtlasStyle();

  // get the stuff from the file
  char fname_postfit[200] = "monoh_withsm_SRCR_bg11.7_bgslop-0.0_nsig0.0.root";
  TFile *tf = new TFile(fname_postfit);
  RooWorkspace *ws_post = (RooWorkspace*)tf->Get("wspace");

  // data and pdf
  RooAbsData *data = ws_post->data("data");
  RooAbsPdf *pdfc = ws_post->pdf("jointModeld");
  
  // variables to be constrained
  RooRealVar *mh = ws_post->var("mh");
  RooRealVar *sigma_h = ws_post->var("sigma_h");
  RooRealVar *eff = ws_post->var("eff");
  RooRealVar *theory = ws_post->var("theory");
  RooRealVar *lumi = ws_post->var("lumi");
  RooRealVar *x_mgg = ws_post->var("mgg");
  RooArgSet cas(*mh,*sigma_h,*eff,*theory,*lumi);  

  // redo the fit
  RooFitResult *r = pdfc->fitTo(*data,RooFit::Constrain(cas),RooFit::Save(true));

  // make the frame
  RooPlot *frame = x_mgg->frame();
  TCanvas *tc = new TCanvas("tc","",700,500);
  frame->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  frame->GetYaxis()->SetTitle("Events / 1 GeV");

  // add the data
  data->plotOn(frame,RooFit::Binning(55),RooFit::Name("xdata"),RooFit::DataError(RooAbsData::Poisson));

  frame->SetMinimum(1e-7);
  frame->SetMaximum(8);


  ///  PDFs
  pdfc->plotOn(frame, RooFit::Components("E"), RooFit::LineColor(kRed),   RooFit::Name("xbackground"));
  pdfc->plotOn(frame, RooFit::LineColor(kRed), RooFit::LineStyle(kDashed),RooFit::Name("xtotal"));

  frame->SetMinimum(1e-7);
  frame->SetMaximum(8);

  // Legend
  double x1,y1,x2,y2;
  GetX1Y1X2Y2(tc,x1,y1,x2,y2);
  TLegend *legend_sr=FastLegend(x1+0.02,y2-0.3,x1+0.35,y2-0.02,0.045);
  legend_sr->AddEntry(frame->findObject("xdata"),"Data","LEP");
  legend_sr->AddEntry(frame->findObject("xbackground"),"Background fit","L");
  legend_sr->AddEntry(frame->findObject("xtotal"),"Higgs","L");
  frame->Draw();
  legend_sr->Draw("SAME");

  // descriptive text
  vector<TString> pavetext11;
  pavetext11.push_back("#bf{#it{ATLAS Internal}}");
  pavetext11.push_back("#sqrt{#it{s}} = 8 TeV #scale[0.6]{#int}Ldt = 20.3 fb^{-1}");
  pavetext11.push_back("#it{H}+MET, #it{H #rightarrow #gamma#gamma}, #it{m}_{#it{H}} = 125.4 GeV");

  TPaveText* text11=CreatePaveText(x2-0.47,y2-0.25,x2-0.05,y2-0.05,pavetext11,0.045);
  text11->Draw();
  tc->Print("paper_fit_plot.pdf");
 }
