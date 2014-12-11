#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include <iostream>
#include "TAxis.h"

const int nmass=5;
float masses[nmass] = {1,65,100,500,1000};

const int nmod=8;

char mods[nmod][20] = { "xxhhs","xxhhg5","xgxFhDh","xdxhDh",
			"zp100","zp1000","s100","s1000"};

float fid_eff[nmod][nmass] = 
  {
    { 0.611,0.631,0.674,0.750,0.761},
    { 0.610,0.619,0.663,0.753,0.763},
    { 0.755,0.761,0.765,0.784,0.787},
    { 0.511,0.623,0.656,0.744,0.749},
    { 0.161,0.255,0.356,0.622,0.675},
    { 0.513,0.510,0.507,0.502,0.665},
    { 0.239,0.315,0.382,0.598,0.619},
    { 0.597,0.612,0.581,0.484,0.639}
  };

float reco_eff[nmod][nmass] = 
  {
    { 0.601,0.602,0.593,0.566,0.556},
    { 0.610,0.603,0.599,0.563,0.556},
    { 0.625,0.620,0.623,0.625,0.622},
    { 0.602,0.595,0.592,0.566,0.564},
    { 0.635,0.636,0.627,0.628,0.622},
    { 0.641,0.632,0.638,0.644,0.654},
    { 0.610,0.611,0.620,0.600,0.598},
    { 0.612,0.605,0.611,0.601,0.606}
  };


void plot_reco_effs()
{
  TGraph *tg[nmod];

  TLegend *tl = new TLegend(0.15,0.15,0.4,0.4);
  tl->SetFillColor(10);
  tl->SetTextFont(42);
  tl->SetBorderSize(0);
  

  for (int imod=0;imod<nmod;imod++)
    {


      tg[imod] = new TGraph(nmass,masses,reco_eff[imod]);

      //      tg[imod]->SetLineColor(2);
      //tg[imod]->SetMarkerStyle(200);

      std::cout << " model " << 
	imod << " graph " << tg[imod] << std::endl<< std::endl<< std::endl<< std::endl<< std::endl;
      tl->AddEntry(tg[imod],mods[imod],"LP");

    }
  
  TCanvas *tc = new TCanvas("reco","",400,400);
  gPad->SetLogx(1);
  
  tg[0]->SetTitle("");
  tg[0]->Draw("ALP");
  tg[0]->SetMaximum(0.7);
  tg[0]->SetMinimum(0.5);
  tg[0]->GetXaxis()->SetTitle("m_{#chi} [GeV]");
  tg[0]->GetYaxis()->SetTitle("Reco Efficiency");

  for (int imod=0;imod<nmod;imod++)
    {
      tg[imod]->Draw("LP");
      tg[imod]->SetMarkerStyle(20);
      tg[imod]->SetMarkerColor(imod+1);
      tg[imod]->SetLineColor(imod+1);
    }
  tl->Draw();
  tc->Print("reco_eff.pdf");
}


void plot_fid_effs()
{
  TGraph *tg[nmod];

  TLegend *tl = new TLegend(0.15,0.15,0.4,0.4);
  tl->SetFillColor(10);
  tl->SetTextFont(42);
  tl->SetBorderSize(0);
  

  for (int imod=0;imod<nmod;imod++)
    {


      tg[imod] = new TGraph(nmass,masses,fid_eff[imod]);

      //      tg[imod]->SetLineColor(2);
      //tg[imod]->SetMarkerStyle(200);

      std::cout << " model " << 
	imod << " graph " << tg[imod] << std::endl<< std::endl<< std::endl<< std::endl<< std::endl;
      tl->AddEntry(tg[imod],mods[imod],"LP");

    }
  
  TCanvas *tc = new TCanvas("reco","",400,400);
  gPad->SetLogx(1);
  
  tg[0]->SetTitle("");
  tg[0]->Draw("ALP");
  tg[0]->SetMaximum(0.8);
  tg[0]->SetMinimum(0.0);
  tg[0]->GetXaxis()->SetTitle("m_{#chi} [GeV]");
  tg[0]->GetYaxis()->SetTitle("Fiducial Efficiency");

  for (int imod=0;imod<nmod;imod++)
    {
      tg[imod]->Draw("LP");
      tg[imod]->SetMarkerStyle(20);
      tg[imod]->SetMarkerColor(imod+1);
      tg[imod]->SetLineColor(imod+1);
    }
  tl->Draw();
  tc->Print("fid_eff.pdf");
}

void plot_lims()
{
  const float sigma_fid = 0.70; // fb

TGraph *tg[nmod];

  TLegend *tl = new TLegend(0.65,0.65,0.87,0.87);
  tl->SetFillColor(10);
  tl->SetTextFont(42);
  tl->SetBorderSize(0);
  

  for (int imod=0;imod<nmod;imod++)
    {
      float xs[nmass];
      for (int im=0;im<nmass;im++)
	xs[im] = sigma_fid/(fid_eff[imod][im]);

      tg[imod] = new TGraph(nmass,masses,xs);

      tl->AddEntry(tg[imod],mods[imod],"LP");

    }
  
  TCanvas *tc = new TCanvas("reco","",400,400);
  gPad->SetLogx(1);
  
  tg[0]->SetTitle("");
  tg[0]->Draw("ALP");
  tg[0]->SetMaximum(5);
  tg[0]->SetMinimum(0.0);
  tg[0]->GetXaxis()->SetTitle("m_{#chi} [GeV]");
  tg[0]->GetYaxis()->SetTitle("95% CL Limit on #sigma_{BSM}");

  for (int imod=0;imod<nmod;imod++)
    {
      tg[imod]->Draw("LP");
      tg[imod]->SetMarkerStyle(20);
      tg[imod]->SetMarkerColor(imod+1);
      tg[imod]->SetLineColor(imod+1);
    }
  tl->Draw();
  tc->Print("xs_lims.pdf");
  
}
