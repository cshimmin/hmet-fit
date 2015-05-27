

#include "StandardFrequentistDiscovery.C"

void pvalue(TString fname)
{
  gROOT->ProcessLine(".L RooDCB.C+");
  StandardFrequentistDiscovery(fname,"wspace","ModelConfig","data",10000);
}
