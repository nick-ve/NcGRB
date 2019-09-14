//////////////////////////////////////////////////////////
// Code (ROOT) to model the GRB neutrino analysis method.
// To run this macro, it is necessary to have the ROOT
// framework (http://root.cern.ch) installed.
// In addition the NCFS-Pack (astro)physics analysis
// framework (http://sites.google.com/site/nickveweb)
// has to be installed as well.
// Both the ROOT and NCFS packages are freely available.
//
// To run the macro and obtain the output in a file called
// "output.log" just issue the command ($ means cmd prompt):
//
// $root -b -q grbmodel.cc >output.log
//
// Author   : Nick van Eijndhoven 17-feb-2006 Utrecht University
// Modified : NvE 01-feb-2007 Utrecht University
// Modified : NvE 14-feb-2007 Utrecht University
// Modified : NvE 26-jul-2007 Utrecht University
// Modified : NvE 07-oct-2007 Utrecht University
//            Variable bin size option implemented.
// Modified : NvE 18-jul-2008 Utrecht University
//            Flat background Psi expectation calculated
//            via the new PsiExtreme facility.
// Modified : NvE 30-mar-2010 IIHE-VUB
//            Adopted the NCFS-Pack namings.
// Modified : NvE 08-nov-2010 IIHE-VUB
//            Explicit Gaussian smearing on theta and phi
//            separately instead of only on opening angle
//            to incorporate detector resolution.
// Modified : NvE 22-apr-2011 IIHE-VUB
//            Various functionalities in the macro more
//            clearly separated and comments added.
//            Attempt to judge Psi quality by using the
//            new features of NcMath::PsiExtreme.            
// Modified : NvE 06-sep-2012 IIHE-VUB
//            Fix to prevent statistical overfluctuation
//            of number of GRB neutrinos in case no 
//            stat. fluctuations are requested (fGrbnu<0).
//            Thanks to Lionel Brayeur for identifying
//            this problem.
// Modified : NvE 25-sep-2014 IIHE-VUB
//            Interfaced with real redshift and T90 measurements
//            and introduced option to use actual GCN data
//            of observed GRBs instead of random sample.
// Modified : NvE 28-oct-2014 IIHE-VUB
//            Corrections for histogram underflow/overflow
//            introduced and also Poisson based dt statistics
//            added to allow an unbinned time series analysis.
// Modified : NvE 15-jun-2015 IIHE-VUB
//            Combination of the GRB location error with the
//            IceCube angular resolution introduced together
//            with the possibility of using GRB specific angular
//            search windows based on this combined spatial resolution. 
//            Also location randomisation and angular smearing
//            performed via the NcAstrolab facilities.
// Modified : NvE 01-jul-2015 IIHE-VUB
//            Number of events (nevtdt) parameter introduced
//            to specify the number of events in a dt cell for which the inter-muon
//            dt statistics will be performed via the Erlang distribution.
//            In this way the dt statistics and Erlang distribution will
//            always be consistent.
//            Thanks Martin Casier and Lionel Brayeur for pointing this out. 
// Modified : NvE 15-jul-2015 IIHE-VUB
//            Usage of the new facility NcAstrolab::GetDxHistogram() introduced.
//            This replaces all the code to create the delta t histograms.
// Modified : NvE 15-feb-2016 IIHE-VUB
//            Posterior background and signal rate PDFs included using the
//            newly introduced NcAstrolab facilities.
// Modified : NvE 21-nov-2016 IIHE-VUB
//            Selection of [zmin,zmax] introduced to enable study of narrow
//            redshift ranges in order to circumvent redshift corrections.
//            Also correct J2000 coordinates stored for GRBs from GCN data.
// Modified : NvE 21-nov-2016 IIHE-VUB
//            Separate filling of "off source" background patches to mimick
//            background observations in corresponding "opposite RA" patches.
//////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>

#include "TStyle.h"
#include "TString.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TChain.h"

#include "NcAstrolab.h"
#include "NcSample.h"
#include "NcMath.h"

void user()
{

 /////////////////////////////////////////////////////
 // User settings to reflect the physical situation //
 /////////////////////////////////////////////////////

 Float_t fDeclmin=5;       // Minimal declination (in degrees) for GRB position acceptance
 Float_t fDeclmax=85;      // Maximal declination (in degrees) for GRB position acceptance
 Float_t fT90min=1e-6;     // Minimal duration (t90 in sec) for GRB acceptance
 Float_t fT90max=1e6;      // Maximal duration (t90 in sec) for GRB acceptance
 Float_t fZmin=-1e-6;         // Minimal redshift for GRB acceptance (<0 : [|fZmin|,fZmax] with random z when redshift is unknown)
 Float_t fZmax=99;         // Maximal redshift for GRB acceptance
 Float_t fSigmagrb=-2.5;    // Angular uncertainty (sigma in degrees) on GRB position (<0 : use GCN data)
 Float_t fMaxsigma=999;    // Maximal combined GRB position and track angular uncertainty (sigma in degrees) for acceptance
 Int_t fNgrbs=-9999;       // Number of (fictative) GRB triggers within GRB declination acceptance (<0 : use GCN data)
 Float_t fGrbnu=-0.15;     // Maximum number of detectable neutrinos per GRB in full (IC86) IceCube (<0 : no stat. fluct.)
 Float_t fAvgrbz=-2;       // Average GRB redshift (<0 : determine from observations)
 Float_t fAvgrbt90=-30;    // Average GRB duration (T90) in seconds (<0 : determine from observations)
 Int_t fInburst=0;         // Flag to indicate that neutrinos are produced coupled (1) or not (0) to the gamma flare duration
 Float_t fDtnu=-60;        // Mean time diff. (in sec) between gamma and nu production (decoupled) or in T90 units w.r.t. trigger (coupled)
 Float_t fDtnus=-0.5;      // Sigma of time difference (in sec) between gamma and nu production (<0 is in T90 units)
 Float_t fAngres=0.5;      // Detector angular resolution (degrees)
 Float_t fTimres=1e-5;     // Detector time resolution (sec)
 Float_t fBkgrate=0.003;   // Mean rate (in Hz) of upgoing bkg muons in full IceCube
 Int_t fNevtdt=2;          // Number of events within a dt cell for which the inter-muon dt statistics will be performed 
 Float_t fDtwin=7200.;     // Total search time window (in sec) centered at GRB trigger
 Float_t fDawin=5;         // Angular search circle (<0 is decl. band) in degrees or sigma around (above/below) GRB position
 Int_t fDatype=0;          // Type of angular window specification (0=in degrees 1=in units of combined GRB/track sigma) 
 Float_t fNbkg=0.5;        // Mean number of counts per bin for auto-binning @@@ auto-renaming clashes with nenbkg and nenbkgdt 
 Int_t fTbint90=1;         // Use time bin size in units of average T90 (1) or not (0) 
 Float_t fTbin=1;          // Time bin size in sec or in av. T90 units (0=variable bins  <0 will result in a mean fNbkg counts/bin)
 Float_t fVarTbin=10;      // Size (in sec) of the first time bin in case of variable time bins
 Float_t fAbin=1;          // Angular bin size in degrees (<0 will result in a mean fNbkg counts/bin)
 Int_t fFreq=0;            // Use frequentist's approximation (1) or exact Bayesian expression (0)
 Int_t fNpsi=0;            // Number of psi entries for bkg psi-value distributions (<0 : time shuffling)
 Int_t fUsetott=1;         // Use the observed tott number of entries in case of time shuffling 
 Int_t fGrbpos=1;          // Use the original grb locations (1) or random ones (0) for bkg studies
 Double_t fNrandom=-1e7;    // (Maximum) number of randomised configurations for direct psi P-value determination (0 means n=1e19)
 Int_t fNcut=10;           // Number of obtained psi>psi0 entries at which randomisations will be terminated (0= no early termination)

 if (fT90min<0) fT90min=0;

 // The data file for observed GRB redshifts
 ifstream zdata;
 zdata.clear();
 zdata.open("GRB-t90-z-Swift.txt");
 if (!zdata.good())
 {
  cout << " *** Data file for observed GRB redshifts not found ***" << endl;
  return;
 }
 zdata.seekg(0); // Position at begin of file

 // The data file for observed GRB t90 durations
 ifstream t90data;
 t90data.clear();
 t90data.open("GRB-t90-Fermi.txt");
 if (!t90data.good())
 {
  cout << " *** Data file for observed GRB T90 durations not found ***" << endl;
  return;
 }
 t90data.seekg(0); // Position at begin of file

 // The GCN TTree data for the IceCube GRBs
 TChain* gcn=new TChain("T");
 gcn->Add("GRB-IC86-*.root");

 Int_t date;
 Float_t ra,decl,err,t90,tstart,tstop,t100,fluence,z;
 Double_t trigmjd;

 // The variables from the Tree
 gcn->SetBranchAddress("date",&date);
 gcn->SetBranchAddress("ra",&ra);
 gcn->SetBranchAddress("decl",&decl);
 gcn->SetBranchAddress("err",&err);
 gcn->SetBranchAddress("t90",&t90);
 gcn->SetBranchAddress("mjd",&trigmjd);
 gcn->SetBranchAddress("t1",&tstart);
 gcn->SetBranchAddress("t2",&tstop);
 gcn->SetBranchAddress("t100",&t100);
 gcn->SetBranchAddress("fluence",&fluence);
 gcn->SetBranchAddress("z",&z);

 // The internal storage array and output file for the produced histograms
 TObjArray histos;
 TFile* fout=new TFile("output.root","RECREATE","GRB analysis results");

 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 // Creation of Swift observed GRB redshift and redshift corrected T90 histos
 TH1F* hz=new TH1F("hz","Swift observed GRB redshifts",200,0,20);
 histos.Add(hz);
 hz->GetXaxis()->SetTitle("GRB redshift");
 hz->GetYaxis()->SetTitle("Counts");

 TH1F* ht90z=new TH1F("ht90z","Swift observed GRB durations with known redshift",100,-5,5);
 histos.Add(ht90z);
 ht90z->GetXaxis()->SetTitle("GRB duration ^{10}log(T90) in sec.");
 ht90z->GetYaxis()->SetTitle("Counts");

 TH2F* ht90z2=new TH2F("ht90z2","Swift observed GRB redshift vs. duration",100,-5,5,200,0,20);
 histos.Add(ht90z2);
 ht90z2->GetXaxis()->SetTitle("GRB duration ^{10}log(T90) in sec.");
 ht90z2->GetYaxis()->SetTitle("GRB redshift");

 TH1F* ht90zc=new TH1F("ht90zc","Swift observed GRB durations redshift corrected",100,-5,5);
 histos.Add(ht90zc);
 ht90zc->GetXaxis()->SetTitle("GRB duration ^{10}log(T90/(z+1)) in sec.");
 ht90zc->GetYaxis()->SetTitle("Counts");

 TH2F* ht90zc2=new TH2F("ht90zc2","Swift observed GRB redshift vs. z corrected duration",100,-5,5,200,0,20);
 histos.Add(ht90zc2);
 ht90zc2->GetXaxis()->SetTitle("GRB z corrected duration ^{10}log(T90) in sec.");
 ht90zc2->GetYaxis()->SetTitle("GRB redshift");

 NcSample zsample;
 zsample.SetStoreMode();
 NcSample t90sample;
 t90sample.SetStoreMode();
 
 // Read the title line in the redshift data file
 string line;
 getline(zdata,line);

 TString name;
 while (zdata >> name >> t90 >> z)
 {
  if (z>=0)
  {
   if (fAvgrbz<0) zsample.Enter(z);
   hz->Fill(z);
   if (t90>0)
   {
    ht90z->Fill(log10(t90));
    ht90z2->Fill(log10(t90),z);
    ht90zc->Fill(log10(t90/(z+1)));
    ht90zc2->Fill(log10(t90/(z+1)),z);
   }
  }
 }

 // Creation of observed GRB t90 duration histo
 TH1F* ht90=new TH1F("ht90","Observed GRB durations",100,-5,5);
 histos.Add(ht90);
 ht90->GetXaxis()->SetTitle("GRB duration ^{10}log(T90) in sec.");
 ht90->GetYaxis()->SetTitle("Counts");

 // Read the title line in the T90 data file
 getline(t90data,line);

 while (t90data >> name >> t90)
 {
  if (t90>0)
  {
   if (fAvgrbt90<0) t90sample.Enter(t90);
   ht90->Fill(log10(t90));
  }
 }

 // Obtaining the GRB locations in IceCube coordinates from the GCN data
 NcAstrolab lab;
 lab.SetExperiment("IceCube");
 lab.SetUT("11-04-2015","12:00:00.0",0); // Fixed fictative analysis date
 Double_t mjd;
 NcTimestamp ts;
 NcSignal* sx=0;
 NcPosition rgrb;
 Float_t t90grb=0;
 Double_t zgrb=0;
 Int_t ngcn=0;
 Int_t jgrb=0;
 TString grbname;
 Float_t sigmagrb=fSigmagrb;
 Float_t totsigma=0;
 Float_t  maxtotsigma=-1.*(sigmagrb*sigmagrb+fAngres*fAngres); // Maximum of the encountered totsigma values

 // Creation of the observed GRB position uncertainty histo
 TH1F* hsigmagrb=new TH1F("hsigmagrb","Observed GRB position uncertainty",450,0,90);
 histos.Add(hsigmagrb);
 hsigmagrb->GetXaxis()->SetTitle("GRB position uncertainty (sigma in degrees)");
 hsigmagrb->GetYaxis()->SetTitle("Counts");

 // Creation of the combined GRB position and track resolution uncertainty histo
 TH1F* htotsigma=new TH1F("htotsigma","Combined GRB position and track resolution uncertainty",450,0,90);
 histos.Add(htotsigma);
 htotsigma->GetXaxis()->SetTitle("Combined GRB position and track resolution uncertainty (sigma in degrees)");
 htotsigma->GetYaxis()->SetTitle("Counts");

 if (fNgrbs<0) // GCN locations requested @@@ Put this in a memberfunction GcnGRBs()
 {
  zsample.Reset();
  t90sample.Reset();

  for (Int_t ient=0; ient<gcn->GetEntries(); ient++)
  {
   gcn->GetEntry(ient);

   // Convert the GCN  error box value into a 1 sigma uncertainty if needed
   if (fSigmagrb<0) sigmagrb=fabs(err);

   // Determine the combined GRB position and track resolution uncertainty.
   totsigma=sigmagrb*sigmagrb+fAngres*fAngres;
   totsigma=sqrt(totsigma);
   
   if (decl<fDeclmin || decl>fDeclmax || totsigma>fMaxsigma) continue;

   t90grb=t90;
   if (t90grb<=0) t90grb=t100;

   if (t90grb<fT90min || t90grb>fT90max) continue;

   zgrb=z;
   if (fZmin<0 && zgrb<0) zgrb=hz->GetRandom();

   if (zgrb<fabs(fZmin) || zgrb>fZmax) continue;

   mjd=trigmjd;

   if (mjd<0) continue;

   hsigmagrb->Fill(sigmagrb);
   htotsigma->Fill(totsigma);

   if (totsigma>maxtotsigma) maxtotsigma=totsigma;

   jgrb++;
   grbname="GRB";
   grbname+=date;
   ts.SetMJD(mjd);
   sx=lab.SetSignal(zgrb,ra,"deg",decl,"deg","equ",&ts,-1,"J",grbname);

   if (!sx) continue;

   sx->AddNamedSlot("t90");
   sx->SetSignal(t90grb,"t90");
   sx->AddNamedSlot("sigmagrb");
   sx->SetSignal(sigmagrb,"sigmagrb");
   sx->AddNamedSlot("totsigma");
   sx->SetSignal(totsigma,"totsigma");

   if (fAvgrbz<0) zsample.Enter(zgrb);
   if (fAvgrbt90<0) t90sample.Enter(t90grb);

   ngcn++;
   if (ngcn>=abs(fNgrbs)) break;
  }
  fNgrbs=lab.GetNRefSignals();
  fNgrbs*=-1;
 }

 // Determine median redshift if requested
 if (fAvgrbz<0)
 {
  fAvgrbz=zsample.GetMedian(1);
  fAvgrbz*=-1.;
 }

 // Determine median T90 duration if requested
 if (fAvgrbt90<0)
 {
  fAvgrbt90=t90sample.GetMedian(1);
  fAvgrbt90*=-1.;
 }

 //////////////////////////////////////////////
 // The implementation of the actual program //
 //////////////////////////////////////////////

 Float_t pi=acos(-1.);

 // Mean number of upgoing bkg muons in search time window
 Float_t mup=fBkgrate*fDtwin;

 // Mean number of upgoing bkg muons per day
 Float_t nmupday=fBkgrate*24.*3600.;

 cout << endl;
 cout << " ============================ User provided settings ==================================" << endl;
 cout << " Declination interval (in degrees) for GRB position acceptance : [" << fDeclmin << "," << fDeclmax << "]" << endl;
 cout << " Duration interval (t90 in sec) for GRB acceptance : [" << fT90min << "," << fT90max << "]" << endl;
 cout << " Redshift interval for GRB acceptance : [" << fabs(fZmin) << "," << fZmax << "]" << endl;
 if (fZmin<0) cout << " Random redshift values taken from z-distribution in case of unknown redshift" << endl;
 if (fSigmagrb>=0) cout << " Fixed GRB position uncertainty (sigma in degrees) : " << sigmagrb << endl;
 cout << " Maximal combined GRB position and track angular uncertainty (sigma in degrees) for acceptance : " << fMaxsigma << endl;
 if (fNgrbs>=0) cout << " Number of randomly generated GRBs within declination acceptance : " << fNgrbs << endl;
 if (fAvgrbz>=0) cout << " User defined average GRB redshift : " << fAvgrbz << endl;
 if (fAvgrbt90>=0) cout << " User defined average GRB T90 duration : " << fAvgrbt90 << endl;
 cout << " Neutrino production coupled (1) or not (0) to the gamma flare duration : " << fInburst << endl;
 if (!fInburst)
 {
  cout << " Mean decoupled time difference (in sec) between GRB gammas and nus : " << fDtnu << endl;
 }
 else
 {
 cout << " Mean coupled time difference (in units of T90 w.r.t. trigger) between GRB gammas and nus : " << fDtnu << endl;
 }
 if (fDtnus>=0)
 {
 cout << " Sigma of mean time difference (in sec) between GRB gammas and nus : " << fDtnus << endl;
 }
 else
 {
 cout << " Sigma of mean time difference (in units of T90) between GRB gammas and nus : " << fabs(fDtnus) << endl;
 }
 cout << " Maximum number (<0 is without stat. fluct.) of detectable neutrinos per GRB : " << fGrbnu << endl;
 cout << " Angular resolution (degrees) of the detector : " << fAngres << endl;
 cout << " Time resolution (sec) of the detector : " << fTimres << endl;
 cout << " Mean rate (Hz) of upgoing bkg muons in the detector : " << fBkgrate << endl;
 cout << " Number of events within a dt cell for which the inter-muon dt statistics will be performed : " << fNevtdt << endl; 
 cout << " Total search time window (in sec) centered at GRB trigger : " << fDtwin << endl;
 if (fDawin>=0)
 {
  if (!fDatype)
  {
   cout << " Angular search circle (in degrees) around the GRB position : " << fDawin << endl;
  }
  else
  {
   cout << " Angular search circle (in combined GRB/track sigma) around the GRB position : " << fDawin << endl;
  }
 }
 else
 {
  if (!fDatype)
  {
   cout << " Angular declination band (in degrees) above/below the GRB position : " << fabs(fDawin) << endl;
  }
  else
  {
   cout << " Angular declination band (in combined GRB/track sigma) above/below the GRB position : " << fabs(fDawin) << endl;
  }
 }
 if (fTbin<0) cout << " Automatic time binning with as mean number of bkg counts/bin : " << fNbkg << endl;
 if (!fTbint90 && fTbin>0) cout << " Time bin size in sec : " << fTbin << endl;
 if (fTbint90 && fTbin>0) cout << " Time bin size in average T90 units : " << fTbin << " (=" << fabs(fTbin*fAvgrbt90) << " sec)" << endl;
 if (!fTbin) cout << " Variable time binning with as size (in sec) for the first time : " << fVarTbin << endl;
 if (fAbin<0)
 {
  cout << " Automatic angular binning with as mean number of bkg counts per bin : " << fNbkg << endl;
 }
 else
 {
  cout << " Angular bin size in degrees : " << fAbin << endl;
 }
 cout << " Use Frequentist's approximation for psi determination : " << fFreq << endl;
 cout << " Number of psi entries for bkg psi-value distributions (<0 means time shuffling) : " << fNpsi << endl;
 cout << " Usage of the observed tott number of entries in case of time shuffling : " << fUsetott << endl; 
 cout << " Usage of actually observed GRB positions for bkg studies : " << fGrbpos << endl;
 if (!fNcut)
 {
  cout << " Number of randomised configurations for direct psi P-value determination : " << fNrandom << endl;
 }
 else
 {
  cout << " Maximum number of randomised configurations for direct psi P-value determination : " << fNrandom << endl;
  cout << " Randomisation will be terminated when " << fNcut << " entries with psi>psi0 have been obtained." << endl;
 }
 cout << " ============================== Derived parameters ====================================" << endl;
 cout << " Mean number of upgoing bkg muons in the time window in the detector : " << mup << endl;
 cout << " Mean number of upgoing bkg muons per day in the detector : " << nmupday << endl;
 if (fNgrbs<0) cout << " Number of accepted GCN observed GRBs : " << abs(fNgrbs) << endl;
 cout << " Number of detected neutrinos per GRB (<0 is without stat. fluct.) : " << fGrbnu << endl;
 if (fAvgrbz<0) cout << " Median GRB redshift from the data sample : " << fabs(fAvgrbz) << endl;
 if (fAvgrbt90<0) cout << " Median GRB T90 duration from the data sample : " << fabs(fAvgrbt90) << endl;
 cout << " ======================================================================================" << endl;
 cout << endl;

 fNgrbs=abs(fNgrbs);

 Float_t danglow=0;    // Lower value (in degrees) of angular difference histo
 Float_t dangup=fDawin; // Upper value (in degrees) of angular difference histo
 if (fDatype) dangup=fDawin*fabs(maxtotsigma);
 if (dangup<0 || dangup>180) dangup=180;

 //////////////////////////////////////////////////////////////////////////
 // Automatic definition of the various signal and background histograms //
 // based on the provided user settings                                  //
 //////////////////////////////////////////////////////////////////////////

 Int_t ntbins=0;
 Double_t* binarr=0;
 if (fabs(fTbin)>0) // Fixed time bins
 {
  if (fTbin>0) // User specified time binning
  {
   if (fTbint90) fTbin=fTbin*fabs(fAvgrbt90); 
   ntbins=int(fDtwin/fTbin);
  }
  else // Automatic time binning to get the specified maximal bkg counts per bin
  {
   ntbins=int(mup*float(fNgrbs)/fNbkg);
  }
 }
 else // Variable time bins
 {
  Int_t nbx=int(fDtwin/fVarTbin);
  Float_t gamma=fabs(fAvgrbz)+1.;
  Float_t* bins=new Float_t[nbx];
  ntbins=0;
  Float_t xlow=0,xup=0,size=fVarTbin;
  for (Int_t i=0; i<nbx-1; i++)
  {
   xup=xlow+size;
   if (xup>fDtwin/2.) // Store the last lowerbound
   {
    bins[i]=xlow;
    ntbins++;
    break;
   }
   bins[i]=xlow;
   ntbins++;
   xlow=xup;
   size=xlow*gamma;
  }
  binarr=new Double_t[2*ntbins-1];
  for (Int_t j=ntbins; j>0; j--)
  {
   binarr[ntbins-j]=-bins[j-1];
   binarr[ntbins+j-2]=bins[j-1];
  }
  ntbins=2*ntbins-2;
 }

 // Binning for the opening angle histo
 Int_t nabins=int((dangup-danglow)/fAbin);
 if (fAbin<0) nabins=int(((dangup-danglow)/180.)*mup*float(fNgrbs)/fNbkg);

 // Binning for the cos(opening angle) histo
 Float_t upcos=cos(danglow*pi/180.);
 Float_t lowcos=cos(dangup*pi/180.);
 Int_t nabins2=int((upcos-lowcos)/(1.-cos(fAbin*pi/180.)));
 if (fAbin<0) nabins2=int(((upcos-lowcos)/2.)*mup*float(fNgrbs)/fNbkg);
 if (nabins2<0) nabins2*=-1;

 if (ntbins<2) ntbins=2;
 if (nabins<2) nabins=2;
 if (nabins2<2) nabins2=2;

 gStyle->SetOptStat("e"); // Only display number of entries in stats box

 Float_t tbinfine=0.1; // Bin size (in sec) for the fine binned histos
 Int_t ntbinsfine=int(fDtwin/tbinfine);

 TString title,s;
 title="t of bkg mu-up in twindow";
 title+=";Upgoing #mu arrival time (in sec) w.r.t. GRB #gamma trigger;Counts per bin of size %-10.3g";
 s=title.Format(title.Data(),tbinfine);
 TH1F* bkgtfine=new TH1F("bkgtfine",s.Data(),ntbinsfine,-fDtwin/2.,fDtwin/2.);
 histos.Add(bkgtfine);

 title="t of all mu-up in twindow";
 title+=";Upgoing #mu arrival time (in sec) w.r.t. GRB #gamma trigger;Counts per bin of size %-10.3g";
 s=title.Format(title.Data(),tbinfine);
 TH1F* tottfine=new TH1F("tottfine",s.Data(),ntbinsfine,-fDtwin/2.,fDtwin/2.);
 histos.Add(tottfine);

 TH1F* bkgt=0;
 TH1F* tott=0;
 TH2F* bkg2=0;
 TH2F* tot2=0;
 if (fabs(fTbin)>0) // Fixed time bins
 {
  bkgt=new TH1F("bkgt","t of bkg mu-up in twindow",ntbins,-fDtwin/2.,fDtwin/2.);
  tott=new TH1F("tott","t of all mu-up in twindow",ntbins,-fDtwin/2.,fDtwin/2.);
  bkg2=new TH2F("bkg2","t vs. dang of bkg mu-up in twindow",
                nabins/10,danglow,dangup,ntbins,-fDtwin/2.,fDtwin/2.);
  tot2=new TH2F("tot2","t vs. dang of all mu-up in twindow",
                nabins/10,danglow,dangup,ntbins,-fDtwin/2.,fDtwin/2.);
 }
 else // Variable time bins
 {
  bkgt=new TH1F("bkgt","t of bkg mu-up in twindow",ntbins,binarr);
  tott=new TH1F("tott","t of all mu-up in twindow",ntbins,binarr);
  bkg2=new TH2F("bkg2","t vs. dang of bkg mu-up in twindow",
                nabins/10,danglow,dangup,ntbins,binarr);
  tot2=new TH2F("tot2","t vs. dang of all mu-up in twindow",
                nabins/10,danglow,dangup,ntbins,binarr);
 }
 histos.Add(bkgt);
 histos.Add(tott);
 histos.Add(bkg2);
 histos.Add(tot2);

 // The opening angle histo
 TH1F* bkga=new TH1F("bkga","dang of bkg mu-up in twindow",nabins,danglow,dangup);
 TH1F* tota=new TH1F("tota","dang of all mu-up in twindow",nabins,danglow,dangup);
 histos.Add(bkga);
 histos.Add(tota);

 // The cos(opening angle) histo
 TH1F* bkgcosa=new TH1F("bkgcosa","cos(dang) of bkg mu-up in twindow",nabins2,lowcos,upcos);
 TH1F* totcosa=new TH1F("totcosa","cos(dang) of all mu-up in twindow",nabins2,lowcos,upcos);
 histos.Add(bkgcosa);
 histos.Add(totcosa);

 Int_t itbin=int(fTbin);
 if (fTbin<0) itbin=int(fDtwin/float(ntbins));
 s="Counts per ";
 if (fabs(fTbin)>0)
 {
  s+=itbin;
  s+=" seconds";
 }
 else
 {
  s+="time bin";
 }
 bkgt->GetXaxis()->SetTitle("Upgoing #mu arrival time (in sec) w.r.t. GRB #gamma trigger");
 bkgt->GetYaxis()->SetTitle(s.Data());
 tott->GetXaxis()->SetTitle("Upgoing #mu arrival time (in sec) w.r.t. GRB #gamma trigger");
 tott->GetYaxis()->SetTitle(s.Data());

 //////////////////////////////////////////////////////////
 // Generation of the signal and background observations //
 // based on the provided user settings                  //
 //////////////////////////////////////////////////////////

 Int_t nmup;
 Double_t thetagrb,phigrb;
 Float_t thetamu,phimu,cost;
 NcRandom rndm;
 Float_t dt=0;
 NcPosition rgrb2; // Actual GRB position from which the neutrinos/muons arrive
 NcPosition rmu;
 Float_t dang;
 Float_t ranlow,ranup;
 Float_t thlow,thup;
 Int_t nmugrb=0;
 NcTimestamp* tx=0;
 Float_t solidangle=0; // Total stacked solid angle 

 // Obtain the (fictative) GRB space-time positions in the declination acceptance
 for (Int_t igrb=0; igrb<fNgrbs; igrb++)
 {
  if (!ngcn) // Generate a random GRB position @@@ Put this in a memberfunction RandomGRBs()
  {
   thlow=fDeclmin+90.;
   thup=fDeclmax+90.;
   if (thup>180) thup=180;
   zgrb=hz->GetRandom();
   rgrb.SetPosition(zgrb,0,0,"sph","deg");
   lab.RandomPosition(rgrb,thlow,thup,0,360);
   thetagrb=rgrb.GetX(2,"sph","deg");
   phigrb=rgrb.GetX(3,"sph","deg");

   t90grb=-1;
   while (t90grb<fT90min || t90grb>fT90max)
   {
    t90grb=ht90->GetRandom();
    t90grb=pow(float(10),t90grb);
   }

   // Determine the combined GRB position and track resolution uncertainty.
   totsigma=sigmagrb*sigmagrb+fAngres*fAngres;
   totsigma=sqrt(totsigma);
   
   if (totsigma>fMaxsigma) continue;

   grbname="Random-GRB";
   grbname+=(igrb+1);
   sx=lab.SetSignal(zgrb,thetagrb,"deg",phigrb,"deg","loc",0,-1,"M",grbname);

   if (!sx) continue;

   sx->AddNamedSlot("t90");
   sx->SetSignal(t90grb,"t90");
   sx->AddNamedSlot("sigmagrb");
   sx->SetSignal(sigmagrb,"sigmagrb");
   sx->AddNamedSlot("totsigma");
   sx->SetSignal(totsigma,"totsigma");
  }
  else // Get GRB position from the GCN data
  {
   sx=lab.GetSignal(igrb+1);

   if (!sx) continue;

   tx=sx->GetTimestamp();
   t90grb=sx->GetSignal("t90");
   sigmagrb=sx->GetSignal("sigmagrb");
   totsigma=sx->GetSignal("totsigma");
   lab.GetSignal(zgrb,thetagrb,"deg",phigrb,"deg","loc",tx,igrb+1);
   rgrb.SetPosition(zgrb,thetagrb,phigrb,"sph","deg");
  } 

   // Update the total stacked solid angle that is probed
   if (fDawin<0) // Declination band
   {
    if (!fDatype)
    {
     thlow=thetagrb-0.5*fabs(fDawin);
     thup=thetagrb+0.5*fabs(fDawin);
    }
    else
    {
     thlow=thetagrb-0.5*fabs(fDawin*totsigma);
     thup=thetagrb+0.5*fabs(fDawin*totsigma);
    }
   }
   else // Circle around GRB position
   {
    if (!fDatype)
    {
     thlow=0;
     thup=fabs(fDawin);
    }
    else
    {
     thlow=0;
     thup=fabs(fDawin*totsigma);
    }
   }

   solidangle+=lab.GetSolidAngle(thlow,thup,"deg",0,360,"deg");    

  // Generate the upgoing bkg muons in the search time window 
  // for both this GRB angular cone and the corresponding "opposite RA" bkg patch
  for (Int_t bkgpatch=0; bkgpatch<=1; bkgpatch++)
  {
   nmup=int(rndm.Poisson(mup));
   for (Int_t imup=0; imup<nmup; imup++)
   {
    ranlow=-fDtwin/2.;
    ranup=fDtwin/2.;
    dt=rndm.Uniform(ranlow,ranup);
    // Smear the time difference with the Gaussian time resolution
    if (fTimres>0) dt=rndm.Gauss(dt,fTimres); //@@@ Is this needed ?
    lab.RandomPosition(rmu,90,180,0,360);
    // Smear the direction of the upgoing bkg muon according to  the detector resolution
    lab.SmearPosition(rmu,fAngres); //@@@ Is this needed
    thetamu=rmu.GetX(2,"sph","deg");
    phimu=rmu.GetX(3,"sph","deg");

    if (fDawin<0) // Declination band
    {
     dang=fabs(thetagrb-thetamu);
    }
    else // Circle around GRB position
    {
     dang=rgrb.GetOpeningAngle(rmu,"deg");
    }

    if ((!fDatype && dang>fabs(fDawin)) || (fDatype && dang>fabs(fDawin*totsigma))) continue;

    if (!bkgpatch)
    {
     tottfine->Fill(dt);
     tott->Fill(dt);
     tota->Fill(dang);
     totcosa->Fill(cos(dang*pi/180.));
     tot2->Fill(dang,dt);
    }
    else
    {
     bkgtfine->Fill(dt);
     bkgt->Fill(dt);
     bkga->Fill(dang);
     bkgcosa->Fill(cos(dang*pi/180.));
     bkg2->Fill(dang,dt);
    }
   }
  }

  // Generate the GRB related upgoing muon(s) in the search window.
  // The GRB position gets Gaussian smeared to reflect the actual position.
  // The time difference between the gammas and the neutrinos gets corrected
  // for the GRB redshift and smeared by the detector time resolution.
  // The muon direction gets Gaussian smeared by the detector angular resolution.

  // Prevent statistical overfluctuation in number of GRB muons if requested by fGrbnu<0
  if (fGrbnu<0 && nmugrb>=int(fabs(fGrbnu)*float(fNgrbs))) continue;

  // Obtain actual GRB position
  rgrb2.Load(rgrb);
  lab.SmearPosition(rgrb2,sigmagrb);

  nmup=int(fabs(fGrbnu));
  if (!nmup && rndm.Uniform()<fabs(fGrbnu)) nmup=1;
  for (Int_t imup2=0; imup2<nmup; imup2++)
  {
   nmugrb++;
   if (!fInburst) // Neutrino and gamma production decoupled
   {
    if (fDtnus<0) // Sigma in units of T90
    {
     dt=rndm.Gauss(fDtnu,fabs(fDtnus)*t90grb);
    }
    else
    {
     dt=rndm.Gauss(fDtnu,fDtnus);
    }
    dt=dt*(zgrb+1.);
   }
   else // Coupled neutrino and gamma production
   {
    if (fDtnus<0) // Sigma in units of T90
    {
     dt=rndm.Gauss(fDtnu*t90grb,fabs(fDtnus)*t90grb);
    }
    else
    {
     dt=rndm.Gauss(fDtnu*t90grb,fDtnus);
    }
   }
   if (fTimres>0) dt=rndm.Gauss(dt,fTimres);
   rmu.Load(rgrb2);
   // Smear the direction of the upgoing GRB muon according to the detector resolution
   if (fAngres>0) lab.SmearPosition(rmu,fAngres);

   // Determine angular difference w.r.t. the presumed GRB position
   dang=rgrb.GetOpeningAngle(rmu,"deg");

   if ((!fDatype && dang>fabs(fDawin)) || (fDatype && dang>fabs(fDawin*totsigma))) continue;

   tottfine->Fill(dt);
   tott->Fill(dt);
   tota->Fill(dang);
   totcosa->Fill(cos(dang*pi/180.));
   tot2->Fill(dang,dt);
  }
 } // End of loop over the individual GRBs

 // Compensate statistical underfluctuation in number of GRB muons if requested by fGrbnu<0
 if (fGrbnu<0)
 {
  nmup=int(fabs(fGrbnu)*float(fNgrbs));
  while (nmugrb<nmup)
  {
   // Pick randomly one of the stored GRBs
   jgrb=int(rndm.Uniform(0.,float(fNgrbs)));
   if (jgrb==0) jgrb=1;
   sx=lab.GetSignal(jgrb);

   if (!sx) continue;

   tx=sx->GetTimestamp();
   t90grb=sx->GetSignal("t90");
   lab.GetSignal(zgrb,thetagrb,"deg",phigrb,"deg","loc",tx,jgrb);
   rgrb.SetPosition(zgrb,thetagrb,phigrb,"sph","deg");
   sigmagrb=sx->GetSignal("sigmagrb");

   // Obtain actual GRB position
   rgrb2.Load(rgrb);
   lab.SmearPosition(rgrb2,sigmagrb); //@@@

   nmugrb++;
   if (!fInburst) // Neutrino and gamma production decoupled
   {
    if (fDtnus<0) // Sigma in units of T90
    {
     dt=rndm.Gauss(fDtnu,fabs(fDtnus)*t90grb);
    }
    else
    {
     dt=rndm.Gauss(fDtnu,fDtnus);
    }
    dt=dt*(zgrb+1.);
   }
   else // Coupled neutrino and gamma production
   {
    if (fDtnus<0) // Sigma in units of T90
    {
     dt=rndm.Gauss(fDtnu*t90grb,fabs(fDtnus)*t90grb);
    }
    else
    {
     dt=rndm.Gauss(fDtnu*t90grb,fDtnus);
    }
   }
   if (fTimres>0) dt=rndm.Gauss(dt,fTimres);
   rmu.Load(rgrb2);
   // Smear the direction of the upgoing GRB muon according to the detector resolution
   if (fAngres>0) lab.SmearPosition(rmu,fAngres);

   // Determine angular difference w.r.t. the presumed GRB position
   dang=rgrb.GetOpeningAngle(rmu,"deg");

   if ((!fDatype && dang>fabs(fDawin)) || (fDatype && dang>fabs(fDawin*totsigma))) continue;

   tottfine->Fill(dt);
   tott->Fill(dt);
   tota->Fill(dang);
   totcosa->Fill(cos(dang*pi/180.));
   tot2->Fill(dang,dt);
  }
 } 

 lab.ListSignals("equ","J",1,"T",10);
 cout << endl;

 /////////////////////////////////////////////////////////////////
 // Creation of the inter-muon time histogram                   //
 // to reflect possible concentration(s) of muon arrival times. //
 // The fine time histos are used for the basic input.          //
 /////////////////////////////////////////////////////////////////

 NcMath math;
 
 ///////////////////////////////////
 // Create the delta t histograms //
 ///////////////////////////////////

 TH1F* htotdt=(TH1F*)(lab.GetDxHistogram(tottfine,fNevtdt,-1,0,-1).Clone("htotdt"));
 histos.Add(htotdt);
 Double_t deltatbin=htotdt->GetXaxis()->GetBinWidth(1);
 Double_t deltatmin=htotdt->GetXaxis()->GetXmin();
 Double_t deltatmax=htotdt->GetXaxis()->GetXmax();
cout << " htotdt. deltatbin : " << deltatbin << " deltatmin : " << deltatmin << " deltatmax : " << deltatmax << endl;
 TH1F* hbkgdt=(TH1F*)(lab.GetDxHistogram(bkgtfine,fNevtdt,deltatbin,deltatmin,deltatmax).Clone("hbkgdt"));
 histos.Add(hbkgdt);
 deltatbin=hbkgdt->GetXaxis()->GetBinWidth(1);
 deltatmin=hbkgdt->GetXaxis()->GetXmin();
 deltatmax=hbkgdt->GetXaxis()->GetXmax();
cout << " hbkgdt. deltatbin : " << deltatbin << " deltatmin : " << deltatmin << " deltatmax : " << deltatmax << endl;

 // Create titles and labels for the delta t histograms
 title="dt to contain ";
 title+=fNevtdt;
 title+=" mu-up in bkg twindow";
 title+=";Upgoing #mu dt (in sec);Counts per bin of size %-10.3g";
 s=title.Format(title.Data(),deltatbin);

 hbkgdt->SetTitle(s.Data());

 title="dt to contain ";
 title+=fNevtdt;
 title+=" mu-up in total twindow";
 title+=";Upgoing #mu dt (in sec);Counts per bin of size %-10.3g";
 s=title.Format(title.Data(),deltatbin);

 htotdt->SetTitle(s.Data());

 // Determination of total and background event rates
 Int_t nbt=tott->GetNbinsX();
 Int_t nba=tota->GetNbinsX();
 Int_t nba2=totcosa->GetNbinsX();
 Float_t underflow, overflow;
 Float_t nentot=tott->GetEntries();
 underflow=tott->GetBinContent(0);
 overflow=tott->GetBinContent(nbt+1);
 nentot=nentot-(underflow+overflow);
 Float_t nenbkg=bkgt->GetEntries();
 underflow=bkgt->GetBinContent(0);
 overflow=bkgt->GetBinContent(nbt+1);
 nenbkg=nenbkg-(underflow+overflow);

 Float_t ratetot=nentot/(fDtwin);
 Float_t ratebkg=nenbkg/(fDtwin);

 // The Bayesian posterior background and signal rate PDFs
 Double_t Non=nentot;
 Double_t Ton=fDtwin*float(fNgrbs);
 Double_t Noff=nenbkg;
 Double_t Toff=Ton;
 TF1 fbkgrpdf=lab.GetBackgroundRatePDF(Noff,Toff);
 TF1 fsigrpdf=lab.GetSignalRatePDF(Non,Ton,Noff,Toff);

 // Provide the 90% credible interval for the signal rate
 Double_t rlow,rup;
 Float_t frac;
 frac=lab.GetCredibleInterval(fsigrpdf,90,rlow,rup);
 cout << " === The " << frac << "% credible interval from the Bayesian signal pdf :"
      << " rlow=" << rlow << " rup=" << rup << endl;

 // Provide the rate PDFs as histograms in the output file
 fbkgrpdf.SetRange(0,3.*Noff/Toff);
 fbkgrpdf.SetNpx(10000);
 TH1* hpdfbkgr=fbkgrpdf.GetHistogram();
 histos.Add(hpdfbkgr);
 hpdfbkgr->SetName("hpdfbkgr");
 hpdfbkgr->Write();
 fsigrpdf.SetRange(0,3.*Non/Ton);
 fsigrpdf.SetNpx(10000);
 TH1* hpdfsigr=fsigrpdf.GetHistogram();
 histos.Add(hpdfsigr);
 hpdfsigr->SetName("hpdfsigr");
 hpdfsigr->Write();

 // Provide the 90% credible interval for the signal rate
 frac=lab.GetCredibleInterval(hpdfsigr,90,rlow,rup);
 cout << " === The " << frac << "% credible interval from the Bayesian signal pdf :"
      << " rlow=" << rlow << " rup=" << rup << endl;

/***********
@@@@@@@@@@@@@@@@@@@@@@@@
// Just for checking
 Int_t nbdt=hbkgdt->GetNbinsX();
 Float_t nenbkgdt=hbkgdt->GetEntries();
 underflow=hbkgdt->GetBinContent(0);
 overflow=hbkgdt->GetBinContent(nbdt+1);
 nenbkgdt=nenbkgdt-(underflow+overflow);
 nbdt=htotdt->GetNbinsX();
 Float_t nentotdt=htotdt->GetEntries();
 underflow=htotdt->GetBinContent(0);
 overflow=htotdt->GetBinContent(nbdt+1);
 nentotdt=nentotdt-(underflow+overflow);

cout << " hbkgdt entries : " << hbkgdt->GetEntries() << " nenbkgdt : " << nenbkgdt << endl;
cout << " htotdt entries : " << htotdt->GetEntries() << " nentotdt : " << nentotdt << endl;
@@@@@@@@@@@@@@@@@@@@@@@@
******************/

 // Creation of the Poisson based dt PDFs for a background only hypothesis 
 TF1 fdttot=math.PoissonDtDist(ratetot,fNevtdt); // Only for reference, not used in the analysis
 TF1 fdtbkg=math.PoissonDtDist(ratebkg,fNevtdt);

 // Provide the dt PDFs as histograms in the output file
 fdttot.SetRange(0,deltatmax);
 fdttot.SetNpx(10000);
 TH1* hpdftotdt=fdttot.GetHistogram();
 histos.Add(hpdftotdt);
 hpdftotdt->SetName("hpdftotdt");
 hpdftotdt->Write();
 fdtbkg.SetRange(0,deltatmax);
 fdtbkg.SetNpx(10000);
 TH1* hpdfbkgdt=fdtbkg.GetHistogram();
 histos.Add(hpdfbkgdt);
 hpdfbkgdt->SetName("hpdfbkgdt");
 hpdfbkgdt->Write();

 ////////////////////////////////////////////////////////////////////////////////
 // Statistical evaluation of the generated signal and background observations //
 //                                                                            //
 // Determination of the Bayesian psi value for the time and angular histos    //
 // under the assumption that there is no GRB signal.                          //
 // This corresponds to searching out the Bernoulli class B_m                  //
 // with m=nbins of the histogram.                                             //
 // An orthodox chi-squared analysis is also performed.                        //
 ////////////////////////////////////////////////////////////////////////////////

/*****
 NcMath math;

 Int_t nbt=tott->GetNbinsX();
 Int_t nba=tota->GetNbinsX();
 Int_t nba2=totcosa->GetNbinsX();
 Float_t underflow, overflow;
 Float_t nentot=tott->GetEntries();
 underflow=tott->GetBinContent(0);
 overflow=tott->GetBinContent(nbt+1);
 nentot=nentot-(underflow+overflow);
 Float_t nenbkg=bkgt->GetEntries();
 underflow=bkgt->GetBinContent(0);
 overflow=bkgt->GetBinContent(nbt+1);
 nenbkg=nenbkg-(underflow+overflow);
 Float_t ratetot=nentot/(fDtwin);
 Float_t ratebkg=nenbkg/(fDtwin);
*****/

 // Statistics of the stacked event samples
 cout << " *** Statistics of the stacked observed event samples ***" << endl;
 cout << " Total stacked solid angle (in sr) : " << solidangle << endl; 
 cout << " *On source* Number of entries : " << nentot << " Number of time bins : " << nbt << " Number of angular bins : " << nba 
      << " Number of cos(angle) bins : " << nba2 << endl;
 cout << " Stacked event rate (Hz) : " << ratetot << " --> Event rate (Hz) per GRB : " << ratetot/float(fNgrbs) << endl;
 cout << " --- (Unknown) Number of GRB muons : " << nmugrb << " Number of \"on source\" bkg entries : " << (nentot-nmugrb) << endl;
 cout << " *Off source* Number of bkg entries : " << nenbkg << endl;
 cout << " Stacked bkg event rate (Hz) : " << ratebkg << " --> Bkg event rate (Hz) per GRB : " << ratebkg/float(fNgrbs) << endl;
 cout << endl; 

 // Time histo Bayesian statistics
 Double_t psitott=math.PsiValue(tott,0,0,fFreq);
 Double_t psibkgt=math.PsiValue(bkgt,0,0,fFreq);

 // Angular histo Bayesian statistics
 TF1 pdfa("pdfa","sin(x*acos(-1.)/180.)");
 Double_t psitota=math.PsiValue(tota,0,&pdfa,fFreq);
 Double_t psibkga=math.PsiValue(bkga,0,&pdfa,fFreq);
 Double_t psitotcosa=math.PsiValue(totcosa,0,0,fFreq);
 Double_t psibkgcosa=math.PsiValue(bkgcosa,0,0,fFreq);

 Float_t psitdif=psitott-psibkgt;
 Float_t psiadif=psitota-psibkga;
 Float_t psicosadif=psitotcosa-psibkgcosa;

 // Extreme Psi values for a pure background hypothesis of the recorded time and angular entries
 Float_t psimintott=math.PsiExtreme(tott,0,0,-2);
 Float_t psimaxtott=math.PsiExtreme(tott,0,0,-1);
 Float_t psifractott=(psimaxtott-psitott)/(psimaxtott-psimintott);
 Float_t psiminbkgt=math.PsiExtreme(bkgt,0,0,-2);
 Float_t psimaxbkgt=math.PsiExtreme(bkgt,0,0,-1);
 Float_t psifracbkgt=(psimaxbkgt-psibkgt)/(psimaxbkgt-psiminbkgt);
 Float_t psimintota=math.PsiExtreme(tota,0,&pdfa,-2);
 Float_t psimaxtota=math.PsiExtreme(tota,0,&pdfa,-1);
 Float_t psifractota=(psimaxtota-psitota)/(psimaxtota-psimintota);
 Float_t psiminbkga=math.PsiExtreme(bkga,0,&pdfa,-2);
 Float_t psimaxbkga=math.PsiExtreme(bkga,0,&pdfa,-1);
 Float_t psifracbkga=(psimaxbkga-psibkga)/(psimaxbkga-psiminbkga);
 Float_t psimintotcosa=math.PsiExtreme(totcosa,0,0,-2);
 Float_t psimaxtotcosa=math.PsiExtreme(totcosa,0,0,-1);
 Float_t psifractotcosa=(psimaxtotcosa-psitotcosa)/(psimaxtotcosa-psimintotcosa);
 Float_t psiminbkgcosa=math.PsiExtreme(bkgcosa,0,0,-2);
 Float_t psimaxbkgcosa=math.PsiExtreme(bkgcosa,0,0,-1);
 Float_t psifracbkgcosa=(psimaxbkgcosa-psibkgcosa)/(psimaxbkgcosa-psiminbkgcosa);

 // Direct determination of the P-value of the various psi values
 Double_t nrxtott=-1, nrxtota=-1, nrxtotcosa=-1;
 Double_t nrxbkgt=-1, nrxbkga=-1, nrxbkgcosa=-1;
 Double_t pvaluettot=-1, pvalueatot=-1, pvaluecosatot=-1;
 Double_t pvaluetbkg=-1, pvalueabkg=-1, pvaluecosabkg=-1;

 TH1F* hrpsitott=new TH1F("hrpsitott","Random #psi distr. for bkg hypothesis of tott",100,psimintott-1.,psimaxtott+1.);
 TH1F* hrpsibkgt=new TH1F("hrpsibkgt","Random #psi distr. for bkg hypothesis of bkgt",100,psiminbkgt-1.,psimaxbkgt+1.);
 TH1F* hrpsitota=new TH1F("hrpsitota","Random #psi distr. for bkg hypothesis of tota",100,psimintota-1.,psimaxtota+1.);
 TH1F* hrpsibkga=new TH1F("hrpsibkga","Random #psi distr. for bkg hypothesis of bkga",100,psiminbkga-1.,psimaxbkga+1.);
 TH1F* hrpsitotcosa=new TH1F("hrpsitotcosa","Random #psi distr. for bkg hypothesis of totcosa",100,psimintotcosa-1.,psimaxtotcosa+1.);
 TH1F* hrpsibkgcosa=new TH1F("hrpsibkgcosa","Random #psi distr. for bkg hypothesis of bkgcosa",100,psiminbkgcosa-1.,psimaxbkgcosa+1.);
 histos.Add(hrpsitott);
 histos.Add(hrpsibkgt);
 histos.Add(hrpsitota);
 histos.Add(hrpsibkga);
 histos.Add(hrpsitotcosa);
 histos.Add(hrpsibkgcosa);

 pvaluettot=math.PsiPvalue(-1,fNrandom,tott,0,0,fFreq,0,hrpsitott,fNcut,&nrxtott);
 pvalueatot=math.PsiPvalue(-1,fNrandom,tota,0,&pdfa,fFreq,0,hrpsitota,fNcut,&nrxtota);
 pvaluecosatot=math.PsiPvalue(-1,fNrandom,totcosa,0,0,fFreq,0,hrpsitotcosa,fNcut,&nrxtotcosa);
 pvaluetbkg=math.PsiPvalue(-1,fNrandom,bkgt,0,0,fFreq,0,hrpsibkgt,fNcut,&nrxbkgt);
 pvalueabkg=math.PsiPvalue(-1,fNrandom,bkga,0,&pdfa,fFreq,0,hrpsibkga,fNcut,&nrxbkga);
 pvaluecosabkg=math.PsiPvalue(-1,fNrandom,bkgcosa,0,0,fFreq,0,hrpsibkgcosa,fNcut,&nrxbkgcosa);

 cout << " *** Observed Bayesian psi values (in dB) for the hyp. of no GRB signal ***" << endl;
 cout << " --- Observed psi values (in dB) for the \"on source\" stacked patches ---" << endl;
 cout << " psi for tott : " << psitott << " tota : " << psitota << " totcosa : " << psitotcosa << endl;
 cout << " ==> P-value of the observed tott psi : " << pvaluettot << " Used # of randomisations : " << nrxtott << endl;
 cout << " ==> P-value of the observed tota psi : " << pvalueatot << " Used # of randomisations : " << nrxtota << endl;
 cout << " ==> P-value of the observed totcosa psi : " << pvaluecosatot << " Used # of randomisations : " << nrxtotcosa << endl;
 cout << " --- Observed psi values (in dB) for the corresponding \"opposite RA\" stacked \"off source\" bkg patches ---" << endl;
 cout << "     psi for bkgt : " << psibkgt << " bkga : " << psibkga << " bkgcosa : " << psibkgcosa << endl;
 cout << " ==> P-value of the background bkgt psi : " << pvaluetbkg << " Used # of randomisations : " << nrxbkgt << endl;
 cout << " ==> P-value of the background bkga psi : " << pvalueabkg << " Used # of randomisations : " << nrxbkga << endl;
 cout << " ==> P-value of the background bkgcosa psi : " << pvaluecosabkg << " Used # of randomisations : " << nrxbkgcosa << endl;
 cout << " --- Difference between observed \"on source\" and \"off source\" psi values (in dB) ---" << endl;
 cout << "     Delta psi for tott-bkgt : " << psitdif << " tota-bkga : " << psiadif
               << " totcosa-bkgcosa : " << psicosadif << endl;

 cout << " === Extreme Psi values for the case of pure background ===" << endl;
 cout << " *** tott psimin : " << psimintott << " psimax : " << psimaxtott << " (psimax-psi)/range : " << psifractott << endl;
 cout << " --- bkgt psimin : " << psiminbkgt << " psimax : " << psimaxbkgt << " (psimax-psi)/range : " << psifracbkgt << endl;
 cout << " *** tota psimin : " << psimintota << " psimax : " << psimaxtota << " (psimax-psi)/range : " << psifractota << endl;
 cout << " --- bkga psimin : " << psiminbkga << " psimax : " << psimaxbkga << " (psimax-psi)/range : " << psifracbkga << endl;
 cout << " *** totcosa psimin : " << psimintotcosa << " psimax : " << psimaxtotcosa << " (psimax-psi)/range : " << psifractotcosa << endl;
 cout << " --- bkgcosa psimin : " << psiminbkgcosa << " psimax : " << psimaxbkgcosa << " (psimax-psi)/range : " << psifracbkgcosa << endl;
 cout << endl;

 // Poisson statistics for the observed dt spectrum
 // Note :
 // ------
 // Only the bkg dt PDF is used, since this may be obtained from off-source measurements.
 // Using the total dt PDF would artificially lower the sensitivity due to possible signal events.
 Double_t nrxdttot=-1, nrxdtbkg=-1;
 Double_t pvaluedttot=-1, pvaluedtbkg=-1;

/******
 TF1 fdttot=math.PoissonDtDist(ratetot,2);
 TF1 fdtbkg=math.PoissonDtDist(ratebkg,2);

 // Provide the dt PDF's as histograms in the output file
 fdttot.SetRange(0,deltatmax);
 fdttot.SetNpx(10000);
 TH1* hpdftotdt=fdttot.GetHistogram();
 hpdftotdt->SetName("hpdftotdt");
 hpdftotdt->Write();
 fdtbkg.SetRange(0,deltatmax);
 fdtbkg.SetNpx(10000);
 TH1* hpdfbkgdt=fdtbkg.GetHistogram();
 hpdfbkgdt->SetName("hpdfbkgdt");
 hpdfbkgdt->Write();
*****/

 Double_t psitotdt=math.PsiValue(htotdt,0,&fdtbkg,fFreq);
 Double_t psibkgdt=math.PsiValue(hbkgdt,0,&fdtbkg,fFreq);

 Float_t psimintotdt=math.PsiExtreme(htotdt,0,&fdtbkg,-2);
 Float_t psimaxtotdt=math.PsiExtreme(htotdt,0,&fdtbkg,-1);
 Float_t psifractotdt=(psimaxtotdt-psitotdt)/(psimaxtotdt-psimintotdt);
 Float_t psiminbkgdt=math.PsiExtreme(hbkgdt,0,&fdtbkg,-2);
 Float_t psimaxbkgdt=math.PsiExtreme(hbkgdt,0,&fdtbkg,-1);
 Float_t psifracbkgdt=(psimaxbkgdt-psibkgdt)/(psimaxbkgdt-psiminbkgdt);

 TH1F* hrpsitotdt=new TH1F("hrpsitotdt","Random #psi distr. for bkg hypothesis of tot dt",100,psimintotdt-1.,psimaxtotdt+1.);
 TH1F* hrpsibkgdt=new TH1F("hrpsibkgdt","Random #psi distr. for bkg hypothesis of bkg dt",100,psiminbkgdt-1.,psimaxbkgdt+1.);
 histos.Add(hrpsitotdt);
 histos.Add(hrpsibkgdt);

 pvaluedttot=math.PsiPvalue(-1,fNrandom,htotdt,0,&fdtbkg,fFreq,0,hrpsitotdt,fNcut,&nrxdttot);
 pvaluedtbkg=math.PsiPvalue(-1,fNrandom,hbkgdt,0,&fdtbkg,fFreq,0,hrpsibkgdt,fNcut,&nrxdtbkg);
 
 cout << " *** Psi analysis for the Poisson Dt distributions ***" << endl;
 cout << " Stacked bkg only event rate (Hz) : " << ratebkg << " --> psi for totdt : " << psitotdt << endl;
 cout << " ==> P-value of the observed totdt psi : " << pvaluedttot << " Used # of randomisations : " << nrxdttot << endl;
 cout << " --- psi for bkgdt : " << psibkgdt << endl;
 cout << "     ==> P-value of the bkgdt psi : " << pvaluedtbkg << " Used # of randomisations : " << nrxdtbkg << endl;

 cout << endl;

 cout << " === Poisson Dt extreme Psi values for the case of pure background ===" << endl;
 cout << " *** totdt psimin : " << psimintotdt << " psimax : " << psimaxtotdt << " (psimax-psi)/range : " << psifractotdt << endl;
 cout << " ---  bkgdt psimin : " << psiminbkgdt << " psimax : " << psimaxbkgdt << " (psimax-psi)/range : " << psifracbkgdt << endl;
 cout << endl;
 
 // The time conventional chi-squared evaluation
 Int_t ndftott;
 Float_t chitott=math.Chi2Value(tott,0,0,&ndftott);
 Int_t ndfbkgt;
 Float_t chibkgt=math.Chi2Value(bkgt,0,0,&ndfbkgt);

 // The angular conventional chi-squared evaluation
 Int_t ndftota;
 Float_t chitota=math.Chi2Value(tota,0,&pdfa,&ndftota);
 Int_t ndfbkga;
 Float_t chibkga=math.Chi2Value(bkga,0,&pdfa,&ndfbkga);
 Int_t ndftotcosa;
 Float_t chitotcosa=math.Chi2Value(totcosa,0,0,&ndftotcosa);
 Int_t ndfbkgcosa;
 Float_t chibkgcosa=math.Chi2Value(bkgcosa,0,0,&ndfbkgcosa);

 cout << endl;
 cout << " *** Observed chi-squared values for the hypothesis of no GRB signal ***" << endl;
 cout << " chi2 for tott : " << chitott << " ndf : " << ndftott
      << " P-value : " << math.Chi2Pvalue(chitott,ndftott)
      << " (" << math.Chi2Pvalue(chitott,ndftott,0,1) << " * sigma)" << endl;
 cout << " chi2 for tota : " << chitota << " ndf : " << ndftota
      << " P-value : " << math.Chi2Pvalue(chitota,ndftota)
      << " (" << math.Chi2Pvalue(chitota,ndftota,0,1) << " * sigma)" << endl;
 cout << " chi2 for totcosa : " << chitotcosa << " ndf : " << ndftotcosa
      << " P-value : " << math.Chi2Pvalue(chitotcosa,ndftotcosa)
      << " (" << math.Chi2Pvalue(chitotcosa,ndftotcosa,0,1) << " * sigma)" << endl;
 cout << endl;
 cout << " --- (Unknown) Corresponding chi-squared values for the bkg upgoing muons ---" << endl;
 cout << " chi2 for bkgt : " << chibkgt << " ndf : " << ndfbkgt
      << " P-value : " << math.Chi2Pvalue(chibkgt,ndfbkgt)
      << " (" << math.Chi2Pvalue(chibkgt,ndfbkgt,0,1) << " * sigma)" << endl;
 cout << " chi2 for bkga : " << chibkga << " ndf : " << ndfbkga
      << " P-value : " << math.Chi2Pvalue(chibkga,ndfbkga)
      << " (" << math.Chi2Pvalue(chibkga,ndfbkga,0,1) << " * sigma)" << endl;
 cout << " chi2 for bkgcosa : " << chibkgcosa << " ndf : " << ndfbkgcosa
      << " P-value : " << math.Chi2Pvalue(chibkgcosa,ndfbkgcosa)
      << " (" << math.Chi2Pvalue(chibkgcosa,ndfbkgcosa,0,1) << " * sigma)" << endl;

 ////////////////////////////////////////////// 
 // Statistics of bkg psi-value distribution //
 ////////////////////////////////////////////// 

 if (fNpsi)
 {
  cout << endl;
  cout << " +++ Simulating GRB background measurements +++" << endl;
  if (fNpsi>0)
  {
   cout << " The above analysis will be repeated (off-burst) " << fNpsi << " times." << endl;
  }
  else
  {
   if (!fUsetott)
   {
    cout << " The above GRB samples will be taken (off-burst) only once more." << endl;
    cout << " By random re-filling the obtained bkg time entries we construct " << abs(fNpsi) << " bkg samples." << endl;
   }
   else
   {
    cout << " By random re-filling the original tott time entries we construct " << abs(fNpsi) << " bkg samples." << endl;
   }
  }
  cout << endl;

  TH1F* hpsibkgt=0;
  if (fabs(fTbin)>0)
  {
   hpsibkgt=new TH1F("hpsibkgt","t of bkg mu-up in twindow",ntbins,-fDtwin/2.,fDtwin/2.);
  }
  else
  {
   hpsibkgt=new TH1F("hpsibkgt","t of bkg mu-up in twindow",ntbins,binarr);
  }
  TH1F* hpsibkga=new TH1F("hpsibkga","dang of bkg mu-up in twindow",nabins,danglow,dangup);
  TH1F* hpsit=new TH1F("hpsit","time bkg psi-value distribution",100,0,2.*psitott);
  TH1F* hpsia=new TH1F("hpsia","angular bkg psi-value distribution",100,0,2.*psitota);
  TObjArray* bkgthists=new TObjArray();
  TObjArray* bkgahists=new TObjArray();
  Float_t bkgpsit=0,bkgchit=0;
  Float_t bkgpsia=0,bkgchia=0;
  NcSample psit,psia,chit,chia;
  psit.SetStoreMode(1);
  psia.SetStoreMode(1);
  chit.SetStoreMode(1);
  chia.SetStoreMode(1);
  Int_t nloop=fNpsi,nshuffle=1;
  if (fNpsi<0)
  {
   nloop=1;
   nshuffle=abs(fNpsi);
  }
  for (Int_t ipsi=0; ipsi<nloop; ipsi++)
  {
   hpsibkgt->Reset();
   hpsibkga->Reset();

   // Generate the (fictative) GRB space-time positions in the declination acceptance.
   // In case the fGrbpos flag was activated, the above GRB positions will be used.
   for (Int_t igrb=0; igrb<fNgrbs; igrb++)
   {
    if (fGrbpos)
    {
     sx=lab.GetSignal(igrb+1);

     if (!sx) continue;

     tx=sx->GetTimestamp();
     totsigma=sx->GetSignal("totsigma");
     lab.GetSignal(zgrb,thetagrb,"deg",phigrb,"deg","loc",tx,igrb+1);
     rgrb.SetPosition(zgrb,thetagrb,phigrb,"sph","deg");
    }
    else
    {
     thlow=fDeclmin+90.;
     thup=fDeclmax+90.;
     if (thup>180) thup=180;
     zgrb=hz->GetRandom();
     rgrb.SetPosition(zgrb,0.,0.,"sph","deg");
     lab.RandomPosition(rgrb,thlow,thup,0,360);
     thetagrb=rgrb.GetX(2,"sph","deg");
     phigrb=rgrb.GetX(3,"sph","deg");

     // Determine the combined GRB position and track resolution uncertainty.
     totsigma=sigmagrb*sigmagrb+fAngres*fAngres;
     totsigma=sqrt(totsigma);
    }

    // Generate the upgoing bkg muons in the search window for this GRB
    nmup=int(rndm.Poisson(mup));
    for (Int_t jmup=0; jmup<nmup; jmup++)
    {
     ranlow=-fDtwin/2.;
     ranup=fDtwin/2.;
     dt=rndm.Uniform(ranlow,ranup);
     // Smear the time difference with the Gaussian time resolution
     if (fTimres>0) dt=rndm.Gauss(dt,fTimres);
     lab.RandomPosition(rmu,90,180,0,360);
     // Smear the direction of the upgoing bkg muon according to the detector resolution
     lab.SmearPosition(rmu,fAngres);
     thetamu=rmu.GetX(2,"sph","deg");
     phimu=rmu.GetX(3,"sph","deg");

     if (fDawin<0) // Declination band
     {
      dang=fabs(thetagrb-thetamu);
     }
     else // Circle around GRB position
     {
      dang=rgrb.GetOpeningAngle(rmu,"deg");
     }

     if ((!fDatype && dang>fabs(fDawin)) || (fDatype && dang>fabs(fDawin*totsigma))) continue;

     hpsibkgt->Fill(dt);
     hpsibkga->Fill(dang);
    }
   } // End of nburst loop

   // Get the corresponding psi values for this bkg simulation

   // The angular location histos
   bkgpsia=math.PsiValue(hpsibkga,0,&pdfa,fFreq);
   bkgchia=math.Chi2Value(hpsibkga,0,&pdfa);

   bkgahists->Add(new TH1F(*hpsibkga));

   // The arrival time histos

   // Refill the time histo by a random re-distribution
   // of the same number of entries as tott
   if (fUsetott && fNpsi<0)
   {
    hpsibkgt->Reset();
    nenbkg=tott->GetEntries();
    for (Int_t ient=0; ient<nenbkg; ient++)
    {
     ranlow=-fDtwin/2.;
     ranup=fDtwin/2.;
     dt=rndm.Uniform(ranlow,ranup);
     hpsibkgt->Fill(dt);
    }
   }

   for (Int_t ishuffle=0; ishuffle<nshuffle; ishuffle++)
   {
    bkgpsit=math.PsiValue(hpsibkgt,0,0,fFreq);
    bkgchit=math.Chi2Value(hpsibkgt);

    bkgthists->Add(new TH1F(*hpsibkgt));
         
    psit.Enter(bkgpsit);
    chit.Enter(bkgchit);
    hpsit->Fill(bkgpsit);

    // Refill the time histo by a random re-distribution of the same number of entries
    nenbkg=hpsibkgt->GetEntries();
    hpsibkgt->Reset();
    for (Int_t ient=0; ient<nenbkg; ient++)
    {
     ranlow=-fDtwin/2.;
     ranup=fDtwin/2.;
     dt=rndm.Uniform(ranlow,ranup);
     hpsibkgt->Fill(dt);
    }
   } // End of shuffle loop

   psia.Enter(bkgpsia);
   chia.Enter(bkgchia);
   hpsia->Fill(bkgpsia);
  } // End of psi loop

  Float_t psitmean=psit.GetMean(1);
  Float_t psitsigma=psit.GetSigma(1);
  Float_t psitmedian=psit.GetMedian(1);
  Float_t psitspread=psit.GetSpread(1);
  Float_t psitdiff=psitott-psitmean;
  Float_t psitdiff2=psitott-psitmedian;
  Float_t psiamean=psia.GetMean(1);
  Float_t psiasigma=psia.GetSigma(1);
  Float_t psiamedian=psia.GetMedian(1);
  Float_t psiaspread=psia.GetSpread(1);
  Float_t psiadiff=psitota-psiamean;
  Float_t psiadiff2=psitota-psiamedian;
  cout << " *** Observed Bayesian bkg psi-value (in dB) statistics for " << abs(fNpsi) << " entries ***" << endl;
  cout << " Time bkg psi distr. Mean : " << psitmean << " Sigma : " << psitsigma
       << " Median : " << psitmedian << " Spread : " << psitspread << endl;
  if (fNpsi>0)
  {
   cout << " Ang. bkg psi distr. Mean : " << psiamean << " Sigma : " << psiasigma
        << " Median : " << psiamedian << " Spread : " << psiaspread << endl;
  }
  cout << " *** Comparison with GRB observed psi-values (in dB) ***" << endl;
  cout << " Time psi-psimean : " << psitdiff;
  if (psitdiff && psitsigma>0) cout << " (" << fabs(psitdiff/psitsigma) << " * sigma)";
  cout << " psi-psimedian : " << psitdiff2;
  if (psitdiff2 && psitspread>0) cout << " (" << fabs(psitdiff2/psitspread) << " * spread)";
  cout << endl;
  if (psitdiff && psitsigma>0)
  {
   cout << " ===> Two sided Gaussian P-value of psi w.r.t. bkg psimean  : " << math.GaussPvalue(fabs(psitdiff/psitsigma)) << endl;
  }
  if (fNpsi>0)
  {
   cout << " Ang. psi-psimean : " << psiadiff;
   if (psiadiff && psiasigma>0) cout << " (" << fabs(psiadiff/psiasigma) << " * sigma)";
   cout << " psi-psimedian : " << psiadiff2;
   if (psiadiff2 && psiaspread>0) cout << " (" << fabs(psiadiff2/psiaspread) << " * spread)";
   cout << endl;
   if (psiadiff && psiasigma>0)
   {
    cout << " ===> Two sided Gaussian P-value of psi w.r.t. bkg psimean  : " << math.GaussPvalue(fabs(psiadiff/psiasigma)) << endl;
   }
  }

  Float_t chitmean=chit.GetMean(1);
  Float_t chitsigma=chit.GetSigma(1);
  Float_t chitmedian=chit.GetMedian(1);
  Float_t chitspread=chit.GetSpread(1);
  Float_t chitdiff=chitott-chitmean;
  Float_t chitdiff2=chitott-chitmedian;
  Float_t chiamean=chia.GetMean(1);
  Float_t chiasigma=chia.GetSigma(1);
  Float_t chiamedian=chia.GetMedian(1);
  Float_t chiaspread=chia.GetSpread(1);
  Float_t chiadiff=chitota-chiamean;
  Float_t chiadiff2=chitota-chiamedian;
  cout << " *** Observed bkg chi-squared statistics for " << abs(fNpsi) << " entries ***" << endl;
  cout << " Time bkg chi-squared values Mean : " << chitmean << " Sigma : " << chitsigma
       << " Median : " << chitmedian << " Spread : " << chitspread << endl;
  if (fNpsi>0)
  {
   cout << " Ang. bkg chi-squared values Mean : " << chiamean << " Sigma : " << chiasigma
        << " Median : " << chiamedian << " Spread : " << chiaspread << endl;
  }
  cout << " *** Comparison with GRB observed chi-squared values ***" << endl;
  cout << " Time chi-chimean : " << chitdiff;
  if (chitdiff && chitsigma>0) cout << " (" << fabs(chitdiff/chitsigma) << " * sigma)";
  cout << " chi-chimedian : " << chitdiff2;
  if (chitdiff2 && chitspread>0) cout << " (" << fabs(chitdiff2/chitspread) << " * spread)";
  cout << endl;
  if (fNpsi>0)
  {
   cout << " Ang. chi-chimean : " << chiadiff;
   if (chiadiff && chiasigma>0) cout << " (" << fabs(chiadiff/chiasigma) << " * sigma)";
   cout << " chi-chimedian : " << chiadiff2;
   if (chiadiff2 && chiaspread>0) cout << " (" << fabs(chiadiff2/chiaspread) << " * spread)";
   cout << endl;
  }

  Float_t psibkgtdiff=psibkgt-psitmean;
  Float_t psibkgtdiff2=psibkgt-psitmedian;
  Float_t psibkgadiff=psibkga-psiamean;
  Float_t psibkgadiff2=psibkga-psiamedian;
  Float_t chibkgtdiff=chibkgt-chitmean;
  Float_t chibkgtdiff2=chibkgt-chitmedian;
  Float_t chibkgadiff=chibkga-chiamean;
  Float_t chibkgadiff2=chibkga-chiamedian;
  cout << endl;
  cout << " --- Comparison with (unknown) GRB bkg psi-values (in dB) ---" << endl;
  cout << " Time psibkg-psimean : " << psibkgtdiff;
  if (psibkgtdiff && psitsigma>0) cout << " (" << fabs(psibkgtdiff/psitsigma) << " * sigma)";
  cout << " psibkg-psimedian : " << psibkgtdiff2;
  if (psibkgtdiff2 && psitspread>0) cout << " (" << fabs(psibkgtdiff2/psitspread) << " * spread)";
  cout << endl;
  if (fNpsi>0)
  {
   cout << " Ang. psibkg-psimean : " << psibkgadiff;
   if (psibkgadiff && psiasigma>0) cout << " (" << fabs(psibkgadiff/psiasigma) << " * sigma)";
   cout << " psibkg-psimedian : " << psibkgadiff2;
   if (psibkgadiff2 && psiaspread>0) cout << " (" << fabs(psibkgadiff2/psiaspread) << " * spread)";
   cout << endl;
  }
  cout << " --- Comparison with (unknown) GRB bkg chi-squared values ---" << endl;
  cout << " Time chi-chimean : " << chibkgtdiff;
  if (chibkgtdiff && chitsigma>0) cout << " (" << fabs(chibkgtdiff/chitsigma) << " * sigma)";
  cout << " chi-chimedian : " << chibkgtdiff2;
  if (chibkgtdiff2 && chitspread>0) cout << " (" << fabs(chibkgtdiff2/chitspread) << " * spread)";
  cout << endl;
  if (fNpsi>0)
  {
   cout << " Ang. chi-chimean : " << chibkgadiff;
   if (chibkgadiff && chiasigma>0) cout << " (" << fabs(chibkgadiff/chiasigma) << " * sigma)";
   cout << " chi-chimedian : " << chibkgadiff2;
   if (chibkgadiff2 && chiaspread>0) cout << " (" << fabs(chibkgadiff2/chiaspread) << " * spread)";
   cout << endl;
  }
 }

  cout << endl;
  Int_t nhistos=histos.GetEntries();
  cout << " ====== The following " << nhistos << " histograms have been generated in the output file ======" << endl;
  for (Int_t ih=0; ih<nhistos; ih++)
  {
   TH1* hx=(TH1*)histos.At(ih);
   if (!hx) continue;
   cout << hx->GetName() << " : " << hx->GetTitle() << endl;
  }
  if (fNpsi)
  {
   cout << " Background studies : hpsibkgt hpsibkga hpsit hpsia" << endl;
   cout << " and all hpsibkgt and hpsibkga histos in the TObjArrays bkgthists and bkgahists" << endl;
  }

 fout->Write();
}
