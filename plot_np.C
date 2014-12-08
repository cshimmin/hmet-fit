

const int np=1;
char names[np] = {"lumi"};

void plot_np()
{
  char fname_prefit[200] = "monoh_ws_prefit.root";
  char fname_postfit[200] = "monoh_withsm_SRCR_bg11.7_bgslop-0.0_nsig0.0.root";

  TFile *tf_pre = new TFile(fname_prefit);
  RooWorkSpace *ws_pre = (RooWorkSpace*)tf_pre->Get("wspace");

  TFile *tf_post = new TFile(fname_postfit);
  RooWorkSpace *ws_post = (RooWorkSpace*)tf_post->Get("wspace");

  for (int i=0;i<np;i++)
    {
      std::cout << "var " << names[i] << " lumi | pre = " << ws_pre->var(names[i])->getVal() << " | post = " << ws_pre->var(names[i])->getVal() << std::endl;
    }

}
