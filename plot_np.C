#include "TFile.h"
#include "RooRealVar.h"
#include "RooWorkSpace.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLatex.h"

const int np=5;
char names[np][20] =  {"lumi",    "eff"    ,"theory"    ,"mh"    ,"sigma_h"};
char unames[np][20] = {"unc_lumi","unc_eff","unc_theory","unc_mh","unc_sh"};

void plot_np()
{
  char fname_prefit[200] = "monoh_ws_prefit.root";
  char fname_postfit[200] = "monoh_withsm_SRCR_bg11.7_bgslop-0.0_nsig0.0.root";

  TFile *tf_pre = new TFile(fname_prefit);
  RooWorkspace *ws_pre = (RooWorkspace*)tf_pre->Get("wspace");

  TFile *tf_post = new TFile(fname_postfit);
  RooWorkspace *ws_post = (RooWorkspace*)tf_post->Get("wspace");

  double pull_mean[np];
  double pull_width[np];
  double pull_pos[np];
  double pull_pos_err[np];

  for (int i=0;i<np;i++)
    {
      double pre_val = ws_pre->var(names[i])->getVal();
      double pre_err = ws_pre->var(unames[i])->getVal();
      double post_val = ws_post->var(names[i])->getVal();
      double post_err = ws_post->var(names[i])->getError();

      pull_mean[i] = (post_val-pre_val)/(pre_err);

      pull_width[i] = post_err/pre_err;

      std::cout << "var " << names[i] 
		<< " lumi | pre = " <<  pre_val
		<< " +- " << pre_err
		<< " | post = " << post_val
		<< " +- " << post_err 
		<< " | pull = " << pull_mean[i] << " +- " << pull_width[i]
		<< std::endl;
      pull_pos[i] = i;
      pull_pos_err[i] = 0;
    }

  TCanvas *tc = new TCanvas("tc","",400,400);
  TGraphErrors *tg = new TGraphErrors(np,pull_pos,pull_mean,pull_pos_err,pull_width);
  tg->SetMarkerStyle(20);
  tg->SetMaximum(1.5);
  tg->SetMinimum(-1.5);
  tg->SetTitle("");
  tg->Draw("AP");
  tg->GetYaxis()->SetTitle("Pull");
  tg->GetXaxis()->SetLabelOffset(99);
  tg->GetXaxis()->SetNdivisions(np+1);
  tg->GetXaxis()->SetRangeUser(-1.5,np-0.5);

  for (int i=0;i<np;i++)
    {
      TLatex *t1 = new TLatex(pull_pos[i]+0.15,1.1,names[i]);
      t1->SetTextSize(0.03);
      t1->SetTextFont(42);
      t1->SetTextAngle(90);
      t1->Draw();
      TLatex *t2 = new TLatex(pull_pos[i]+0.15,-1.2,Form("%1.2f #pm %1.2f",pull_mean[i],pull_width[i]));
      t2->SetTextSize(0.03);
      t2->SetTextFont(42);
      t2->SetTextAngle(90);
      t2->Draw();
    }
  tc->Print("pulls.pdf");
  

}
