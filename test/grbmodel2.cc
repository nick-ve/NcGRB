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

 // Load and compile the user code
 gROOT->LoadMacro("NcAstrolab2.cxx+");

 // The (transient) cosmic coincidence analysis
 NcAstrolab2 grb;
 grb.SetExperiment("IceCube");
 grb.SetUT("11-04-2015","12:00:00.0",0); // Fixed fictative analysis date
 grb.SetRandomiser(-1); // Use the UT timestamp to generate a seed
 grb.Data();

// grb.SetBurstParameter("Nmax",20);
 grb.SetBurstParameter("Declmin",5);
 grb.SetBurstParameter("Declmax",85);
// grb.SetParameter("Grbnu",-0.03);
// grb.SetParameter("Tbint90",1.5);

 NcDevice* params=grb.GetBurstParameters();
 if (params) params->Data();

 cout << endl;
 grb.ListBurstParameters();

 grb.MakeBurstZdist("../grbweb/GRB-z-Swift.root","T","z",200,0,20);
 grb.MakeBurstT90dist("../grbweb/GRB-t90-Fermi.root","T","t90");

 grb.LoadBurstGCNdata("../grbweb/GRBweb.root","T");

/***********
// grb.LoadBurstGCNdata("GRB-IC86-*.root","T");
 grb.LoadBurstGCNdata("GRB-IC86-2011.root","T");
 grb.LoadBurstGCNdata("GRB-IC86-2012.root","T");
 grb.LoadBurstGCNdata("GRB-IC86-2013.root","T");
 grb.LoadBurstGCNdata("GRB-IC86-2014.root","T");
 grb.LoadBurstGCNdata("GRB-IC86-2015.root","T");
 grb.LoadBurstGCNdata("GRB-IC86-2016+.root","T");
************/

/// grb.GenBurstGCNdata(5,"GRB");
/// grb.GenBurstGCNdata(5,"GW");

 cout << endl;
 grb.ListSignals("equ","J",1,"T",10);
 cout << endl;

 grb.GenBurstSignals();

 cout << endl;
 grb.ListBurstParameters();

 Double_t rlow,rup;
 grb.GetBurstBayesianSignalRate(90,rlow,rup);

/**************
 grb.GetLiMaSignificance();

 grb.GetBayesianPsiStatistics("time",2,1e4);
 grb.GetBayesianPsiStatistics("angle",2,1e4);
 grb.GetBayesianPsiStatistics("cosa",2,1e4);
 grb.GetBayesianPsiStatistics("dt",2,1e4);

 grb.GetChi2Statistics("time",2);
 grb.GetChi2Statistics("angle",2);
 grb.GetChi2Statistics("cosa",2);
 grb.GetChi2Statistics("dt",2);
***************/

 grb.WriteBurstHistograms("tessie.root");
}
