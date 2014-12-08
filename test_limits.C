// to run from command line:
// $ root -q "test_limits.C();"

#define LUMI 20.2769

void test_limits() {
	// note that RooDCB has to be loaded first!!
	gROOT->ProcessLine(".L RooDCB.C+");
	gROOT->ProcessLine(".L fit_withsm.C+");
	//limit_bands(11.67,-1.0/60.176,0,1.069/LUMI,0.07,0.05);
	limit_bands(11.67,-1.0/60.176,0,1.069/LUMI,0.042,0.029);
}
