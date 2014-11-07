
#include "TFile.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/HypoTestResult.h"

#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooWorkspace.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooExponential.h"
#include "RooUniform.h"
#include "RooFitResult.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "RooGlobalFunc.h"
#include "RooFit.h"
#include "RooCmdArg.h"
#include "TH1.h"

#include "StandardHypoTestInvDemo.C"

void limit(TString fname)
{
  StandardHypoTestInvDemo( fname, "wspace","ModelConfig","","data",2,3,true,50,0,20);
}


#include "fit_withsm.C"
// Calculate the expected limit on N_events given the following parameters:
//
// 
// v_nbg : bg normalization, the 'A' in  A exp(-x/t)
// bg_slope : bg slope, the 't' in A exp(-x/t)
// v_nsig : number of signal events
// v_nsm  : number of SM peaking events
// s_nsm  : constraint on the number of signal events
//
// Example usage 1:
//
//  limit_CR(11.67,-1.0/60.176,0,1.13,0.05)
//
//    v_nbg and bg_slope come from fitting the bg parameterization to MC, for MET>90, pt>90
//    v_nsig = 0 for the bg-only hypothesis
//    v_nsm = 1.13 for the estimate of SM higgs in our window
//    s_nsm = 0.05 means 5% relative uncertainty on 1.13

//  gives this result:
// Expected upper limits, using the B (alternate) model : 
// expected limit (median) 6.81151
// expected limit (-1 sig) 4.64374
// expected limit (+1 sig) 10.2337
// expected limit (-2 sig) 3.33272
// expected limit (+2 sig) 15.0239

// Example usage 2:
//
//  limit_CR(11.67,-1.0/60.176,0,0,0.01) 
//
//    v_nbg and bg_slope come from fitting the bg parameterization to MC, for MET>90, pt>90
//    v_nsig = 0 for the bg-only hypothesis
//    v_nsm = 0 for no peaking SM 
//    s_nsm = 0.01 beacuse it is unstable with 0.00

//  gives this result:
// expected limit (median) 6.27288
// expected limit (-1 sig) 4.25024
// expected limit (+1 sig) 9.50661
// expected limit (-2 sig) 3.04664
// expected limit (+2 sig) 14.0772

void limit_bands( float v_nbg, float bg_slope,float v_nsig,float v_nsm, float s_nsm, float s_eff)
{
  TString fname = fit_withsm(v_nbg,bg_slope,v_nsig,v_nsm,s_nsm,s_eff);
  std::cout << "Reading WS from " << fname << std::endl;
  limit(fname);
}

