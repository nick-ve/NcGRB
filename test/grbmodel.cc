///////////////////////////////////////////////////////////////////
// ROOT macro to process IceCube data in plain TTree format.     //
// The actual C++ code to be processed is contained in a file    //
// called grbmodel.cxx with "void user()" as main function.      //
//                                                               //
// To execute the user code, just invoke ($ = prompt)            //
// $ root -b -q grbmodel.cc                                      //
//                                                               //
//--- Nick van Eijndhoven 03-sep-2013 IIHE-VUB, Brussels         //
///////////////////////////////////////////////////////////////////
{
 gSystem->Load("ncfspack");
 gSystem->Load("icepack");

 // Load and compile the user code
 gROOT->LoadMacro("NcGRB.cxx+");

 // The (transient) cosmic coincidence analysis
 NcGRB grb;
 grb.SetExperiment("IceCube");
 grb.SetUT("11-04-2015","12:00:00.0",0); // Fixed fictative analysis date
 grb.SetRandomiser(-1); // Use the UT timestamp to generate a seed
 grb.Data();

 grb.SetParameter("Declmin",5);
 grb.SetParameter("Declmax",85);
 grb.SetParameter("Grbnu",-0.03);
// grb.SetParameter("Tbint90",1.5);
// grb.SetParameter("Nmax",75);

// grb.SetParameter("Nrandom",1e7);
// grb.SetParameter("Ncut",10);

 cout << endl;
 grb.PrintSettings();

 grb.MakeZdist("../grbweb/GRB-z-Swift.root","T","z",200,0,20);
 grb.MakeT90dist("../grbweb/GRB-t90-Fermi.root","T","t90");

 grb.LoadBurstData("../grbweb/GRBweb.root","T");

 cout << endl;
 grb.ListSignals("equ","J",1,"T",10);
 cout << endl;

/***********
// grb.LoadBurstData("GRB-IC86-*.root","T");
 grb.LoadBurstData("GRB-IC86-2011.root","T");
 grb.LoadBurstData("GRB-IC86-2012.root","T");
 grb.LoadBurstData("GRB-IC86-2013.root","T");
 grb.LoadBurstData("GRB-IC86-2014.root","T");
 grb.LoadBurstData("GRB-IC86-2015.root","T");
 grb.LoadBurstData("GRB-IC86-2016+.root","T");
************/
// grb.GenBurstData(50);
// grb.GenBurstData(50);

 grb.ExecuteTask();

 cout << endl;
 grb.PrintSettings();

 Double_t rlow,rup;
 grb.GetBayesianSignalRate(90,rlow,rup);

 grb.GetLiMaSignificance();

 grb.GetBayesianPsiStatistics("time",2,1e4);
 grb.GetBayesianPsiStatistics("angle",2,1e4);
 grb.GetBayesianPsiStatistics("cosa",2,1e4);
 grb.GetBayesianPsiStatistics("dt",2,1e4);

 grb.GetChi2Statistics("time",2);
 grb.GetChi2Statistics("angle",2);
 grb.GetChi2Statistics("cosa",2);
 grb.GetChi2Statistics("dt",2);

/*************
 // Create a multi-task job environment
 NcJob job;
 job->Add(&grb);
 job->ListEnvironment();

 // Execute the analysis
 job->ExecuteJob(10);
**********/

 grb.WriteHistograms("output.root");

}
