// to run from command line:
// $ root -q "test_limits.C();"

void test_limits() {
	// note that RooDCB has to be loaded first!!
	gROOT->ProcessLine(".L RooDCB.C+");
	gROOT->ProcessLine(".L fit_withsm.C+");
	test();
}
