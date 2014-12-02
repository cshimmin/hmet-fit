/*
void StandardTestStatDistributionDemo(const char* infile = "",
				      const char* workspaceName = "combined",
				      const char* modelConfigName = "ModelConfig",
				      const char* dataName = "obsData") */

void check_asympt(char file[500])
{
  gROOT->ProcessLine(".L RooDCB.C+");
  gROOT->ProcessLine(".L StandardTestStatDistributionDemo.C+");
  StandardTestStatDistributionDemo(file,"wspace","ModelConfig","data");

}
