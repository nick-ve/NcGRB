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
 grb.SetBurstParameter("Grbnu",-0.03);
// grb.SetBurstParameter("Tbint90",1.5);
 grb.SetBurstParameter("Kinangle",3);

 NcDevice* params=grb.GetBurstParameters();
 if (params) params->Data();

 cout << endl;
 grb.ListBurstParameters();

 grb.MakeBurstZdist("../grbweb/GRB-z-Swift.root","T","z",200,0,20);
 grb.MakeBurstT90dist("../grbweb/GRB-t90-Fermi.root","T","t90");
 grb.MakeBurstBkgEdist("IC86*data.root","tree","logE","dec","rad",200,1e7,1000);

 Double_t emin=100;
 Double_t emax=1e7;
 Int_t nbins=10000;
 TF1 spec("spec","pow(x,-2.)");
// grb.MakeBurstEdist(spec,emin,emax,nbins);

 Double_t gamma=2;
 grb.MakeBurstEdist(gamma,emin,emax,nbins);

/***************
 TF1 fx;
 grb.GetNeutrinoAngle(1000,"deg",2,&fx);
 grb.GetNeutrinoAngle(1000,"deg",2,&fx);
 grb.GetNeutrinoAngle(1000,"deg",2,&fx);
fx->Print();
 grb.GetNeutrinoAngle(10000,"deg",2,&fx);
 grb.GetNeutrinoAngle(10000,"deg",2,&fx);
 grb.GetNeutrinoAngle(10000,"deg",2,&fx);
fx->Print();
 grb.GetNeutrinoAngle(100000,"deg",2,&fx);
 grb.GetNeutrinoAngle(100000,"deg",2,&fx);
 grb.GetNeutrinoAngle(100000,"deg",2,&fx);
fx->Print();
 grb.GetNeutrinoAngle(100,"deg",2,&fx);
 grb.GetNeutrinoAngle(100,"deg",2,&fx);
 grb.GetNeutrinoAngle(100,"deg",2,&fx);
fx->Print();
*************/

/*********
 Double_t E=0;
 Double_t ang=0;
 Nc3Vector v;
 v.SetVector(1,90,0,"sph","deg");
 Nc3Vector w;
 w.Load(v);
 cout << " --- Original vector" << endl;
 w.Data("sph","deg");
 cout << " --- Shifted vector" << endl;
 v.Data("sph","deg");
 cout << " Opening angle : " << v.GetOpeningAngle(w,"deg") << " degrees." << endl;
 for (Int_t i=0; i<50; i++)
 {
  E=grb.GetBurstBackgroundEnergy();
  ang=grb.GetNeutrinoAngle(E,"deg");
  cout << " Energy : " << E << " GeV --> Angle : " << ang << " degrees." << endl;
///  grb.ShiftPosition(v,ang);
  grb.SmearPosition(v,ang);
  cout << " --- Original vector" << endl;
  w.Data("sph","deg");
  cout << " --- Shifted vector" << endl;
  v.Data("sph","deg");
  ang=ang-v.GetOpeningAngle(w,"deg");
  cout << " *** Opening angle : " << v.GetOpeningAngle(w,"deg") << " degrees." << endl;
  cout << " *** Kinematic-opening angle : " << ang << " degrees." << endl;
//@@@@  w.Load(v);
  v.Load(w);
 }
*********/

 grb.LoadBurstGCNdata("../grbweb/GRBweb.root","T");

// grb.GenBurstGCNdata(100,"GRB");
// grb.GenBurstGCNdata(5,"GW");

 cout << endl;
 grb.ListSignals("equ","J",1,"T",10);
 cout << endl;

 grb.GenBurstSignals();

 cout << endl;
 grb.ListBurstParameters();

 Double_t rlow,rup;
 grb.GetBurstBayesianSignalRate(90,rlow,rup);

 grb.GetBurstLiMaSignificance();

 grb.GetBurstBayesianPsiStatistics("time",2,1e4);
 grb.GetBurstBayesianPsiStatistics("angle",2,1e4);
 grb.GetBurstBayesianPsiStatistics("cosa",2,1e4);
 grb.GetBurstBayesianPsiStatistics("dt",2,1e4);

 grb.GetBurstChi2Statistics("time",2);
 grb.GetBurstChi2Statistics("angle",2);
 grb.GetBurstChi2Statistics("cosa",2);
 grb.GetBurstChi2Statistics("dt",2);

 grb.WriteBurstHistograms("tessie.root");
}
