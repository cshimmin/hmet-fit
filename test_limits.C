// to run from command line:
// $ root -q "test_limits.C();"

#define LUMI 20.2769

void test_limits() {
	// note that RooDCB has to be loaded first!!
	gROOT->ProcessLine(".L RooDCB.C+");
	gROOT->ProcessLine(".L fit_withsm.C+");

	fitArgs args = DEFAULT_FIT_ARGS;

	limit_level_t lim = MODEL;

	switch (lim) {
	case VISIBLE:
		args.limit_level = VISIBLE;
		break;
	case FIDUCIAL:
		args.limit_level = FIDUCIAL;
		args.v_reco = 0.61;
		break;
	case MODEL:
		args.limit_level = MODEL;
		args.in_unc_sig_theory = 0.40;
		args.v_reco = 0.355;
		break;
	}

	limit_bands(args);
}
