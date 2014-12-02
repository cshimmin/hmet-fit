#include "StandardTestStatDistributionDemo.C"


/*
void StandardTestStatDistributionDemo(const char* infile = "",
				      const char* workspaceName = "combined",
				      const char* modelConfigName = "ModelConfig",
				      const char* dataName = "obsData") */

void check_asympt(char file[500])
{
  StandardTestStatDistributionDemo(file,"wspace","ModelConfig","data");

}
