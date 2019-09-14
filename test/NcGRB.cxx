/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright(c) 1997-2009, NCFS, All Rights Reserved.                          *
 *                                                                             *
 * Author: The Netherlands Center for Fundamental Studies (NCFS).              *
 *         http://sites.google.com/site/ncfsweb ncfs.nl@gmail.com              *
 *                                                                             *
 * Contributors are mentioned in the code where appropriate.                   *
 *                                                                             * 
 * No part of this software may be used, copied, modified or distributed       *
 * by any means nor transmitted or translated into machine language without    *
 * written permission by the NCFS.                                             *
 * Permission to use the documentation strictly for non-commercial purposes    *
 * is hereby granted without fee, provided that the above copyright notice     *
 * appears in all copies and that both the copyright notice and this           *
 * permission notice appear in the supporting documentation.                   *
 * This software is provided "as is" without express or implied warranty.      *
 * The authors make no claims that this software is free of error, or is       *
 * consistent with any particular standard of merchantability, or that it      *
 * will meet your requirements for any particular application, other than      *
 * indicated in the corresponding documentation.                               *
 * This software should not be relied on for solving a problem whose           *
 * incorrect solution could result in injury to a person or loss of property.  *
 * If you do use this software in such a manner, it is at your own risk.       *
 * The authors disclaim all liability for direct or consequential damage       *
 * resulting from your use of this software.                                   *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

///////////////////////////////////////////////////////////////////////////
// Class NcGRB
// NcAstrolab (and TTask) derived class to perform coincidence studies of
// (transient) cosmic phenomena.
//
// Example :
// =========
//
//
//--- Author: Nick van Eijndhoven 17-feb-2006 Utrecht University
//- Modified: Nick van Eijndhoven February 27, 2019 12:36 IIHE-VUB, Brussel
///////////////////////////////////////////////////////////////////////////

#include "NcGRB.h"
#include "Riostream.h"
 
ClassImp(NcGRB) // Class implementation to enable ROOT I/O
 
NcGRB::NcGRB(const char* name,const char* title) : NcAstrolab(name,title)
{
// Constructor and initialisation of default parameters.

 fNmax=-1;         // Maximal number of GRBs to be accepted for analysis (<0 : no limitation)
 fDeclmin=-90;     // Minimal declination (in degrees) for GRB position acceptance
 fDeclmax=90;      // Maximal declination (in degrees) for GRB position acceptance
 fT90min=1e-6;     // Minimal duration (t90 in sec) for GRB acceptance
 fT90max=1e6;      // Maximal duration (t90 in sec) for GRB acceptance
 fZmin=-1e-6;      // Minimal redshift for GRB acceptance (<0 : [|fZmin|,fZmax] with random z when redshift is unknown)
 fZmax=9999;       // Maximal redshift for GRB acceptance
 fSigmagrb=-2.5;   // Angular uncertainty (sigma in degrees) on GRB position (<0 : determine from observations)
 fMaxsigma=999;    // Maximal combined GRB position and track angular uncertainty (sigma in degrees) for acceptance
 fGrbnu=-0.05;     // Maximum number of detectable neutrinos per GRB in (<0 : no stat. fluct.)
 fAvgrbz=-1;       // Average GRB redshift (<0 : determine from observations)
 fAvgrbt90=-1;     // Average GRB duration (T90) in seconds (<0 : determine from observations)
 fInburst=0;       // Flag to indicate that neutrinos are produced coupled (1) or not (0) to the gamma flare duration
 fDtnu=-60;        // Mean time diff. (in sec) between gamma and nu production (decoupled) or in T90 units w.r.t. trigger (coupled)
 fDtnus=-0.5;      // Sigma of time difference (in sec) between gamma and nu production (<0 is in T90 units)
 fAngres=0.5;      // Detector angular resolution (degrees)
 fTimres=1e-5;     // Detector time resolution (sec)
 fBkgrate=0.003;   // Mean rate (in Hz) of upgoing bkg muons
 fNevtdt=2;        // Number of events within a dt cell for which the inter-muon dt statistics will be performed 
 fDtwin=7200.;     // Total search time window (in sec) centered at GRB trigger
 fDawin=5;         // Angular search circle (<0 is decl. band) in degrees or sigma around (above/below) GRB position
 fDatype=0;        // Type of angular window specification (0=in degrees 1=in units of combined GRB/track sigma) 
 fNbkg=0.5;        // Mean number of counts per bin for auto-binning
 fTbint90=1;       // Time bin size in units of average T90 (0 : Time bin size determined by fTbin) 
 fTbin=1;          // Time bin size in seconds (0=variable bins  <0 will result in a mean fNbkg counts/bin)
 fVarTbin=10;      // Size (in sec) of the first time bin in case of variable time bins
 fAbin=1;          // Angular bin size in degrees (<0 will result in a mean fNbkg counts/bin)
 fFreq=0;          // Use frequentist's approximation (1) or exact Bayesian expression (0)
 fNpsi=0;          // Number of psi entries for bkg psi-value distributions (<0 : time shuffling)
 fUsetott=1;       // Use the observed tott number of entries in case of time shuffling 
 fGrbpos=1;        // Use the original burst locations (1) or random ones (0) for bkg studies
 fNrandom=-1e7;    // (Maximum) number of randomised configurations for direct psi P-value determination (0 means n=1e19)
 fNcut=10;         // Number of obtained psi>psi0 entries at which randomisations will be terminated (0= no early termination)

 fNgrbs=0;
 fMaxtotsigma=-1;
 fMup=0;
 fNmupday=0;
}
///////////////////////////////////////////////////////////////////////////
NcGRB::~NcGRB()
{
// Default destructor.

 // Remove the subtasks from the internal TTask list without deleting them
 if (fTasks) fTasks->Clear();

}
///////////////////////////////////////////////////////////////////////////
NcGRB::NcGRB(const NcGRB& q) : NcAstrolab(q)
{
// Copy constructor

 fNmax=q.fNmax;
 fDeclmin=q.fDeclmin;
 fDeclmax=q.fDeclmax;
 fT90min=q.fT90min;
 fT90max=q.fT90max;
 fZmin=q.fZmin;
 fZmax=q.fZmax;
 fSigmagrb=q.fSigmagrb;
 fMaxsigma=q.fMaxsigma;
 fGrbnu=q.fGrbnu;
 fAvgrbz=q.fAvgrbz;
 fAvgrbt90=q.fAvgrbt90;
 fInburst=q.fInburst;
 fDtnu=q.fDtnu;
 fDtnus=q.fDtnus;
 fAngres=q.fAngres;
 fTimres=q.fTimres;
 fBkgrate=q.fBkgrate;
 fNevtdt=q.fNevtdt;
 fDtwin=q.fDtwin;
 fDawin=q.fDawin;
 fDatype=q.fDatype;
 fNbkg=q.fNbkg;
 fTbint90=q.fTbint90;
 fTbin=q.fTbin;
 fVarTbin=q.fVarTbin;
 fAbin=q.fAbin;
 fFreq=q.fFreq;
 fNpsi=q.fNpsi;
 fUsetott=q.fUsetott;
 fGrbpos=q.fGrbpos;
 fNrandom=q.fNrandom;
 fNcut=q.fNcut;

 fNgrbs=q.fNgrbs;
 fMaxtotsigma=q.fMaxtotsigma;
 fMup=q.fMup;
 fNmupday=q.fNmupday;

 TObjArray arr=q.fHistos;
 Int_t n=arr.GetEntries();
 for (Int_t i=0; i<n; i++)
 {
  TH1* h=(TH1*)arr.At(i);
  if (h) fHistos.Add(h->Clone());
 }
}
///////////////////////////////////////////////////////////////////////////
Int_t NcGRB::SetParameter(TString name,Double_t value)
{
// Specification of a certain parameter setting.
// For default values please refer to the constructor of this class.
//
// Input arguments :
// -----------------
// name  : Name of the parameter to be set
// value : The parameter value to be set
//
// The available parameter names are :
//
// Nmax;     // Maximal number of GRBs to be accepted for analysis (<0 : no limitation)
// Declmin   // Minimal declination (in degrees) for GRB position acceptance
// Declmax   // Maximal declination (in degrees) for GRB position acceptance
// T90min    // Minimal duration (t90 in sec) for GRB acceptance
// T90max    // Maximal duration (t90 in sec) for GRB acceptance
// Zmin      // Minimal redshift for GRB acceptance (<0 : [|fZmin|,fZmax] with random z when redshift is unknown)
// Zmax      // Maximal redshift for GRB acceptance
// Sigmagrb  // Angular uncertainty (sigma in degrees) on GRB position (<0 : determine from observations)
// Maxsigma  // Maximal combined GRB position and track angular uncertainty (sigma in degrees) for acceptance
// Grbnu     // Maximum number of detectable neutrinos per GRB (<0 : no stat. fluct.)
// Avgrbz    // Average GRB redshift (<0 : determine from observations)
// Avgrbt90  // Average GRB duration (T90) in seconds (<0 : determine from observations)
// Inburst   // Flag to indicate that neutrinos are produced coupled (1) or not (0) to the gamma flare duration
// Dtnu      // Mean time diff. (in sec) between gamma and nu production (decoupled) or in T90 units w.r.t. trigger (coupled)
// Dtnus     // Sigma of time difference (in sec) between gamma and nu production (<0 is in T90 units)
// Angres    // Detector angular resolution (degrees)
// Timres    // Detector time resolution (sec)
// Bkgrate   // Mean rate (in Hz) of upgoing bkg muons
// Nevtdt    // Number of events within a dt cell for which the inter-muon dt statistics will be performed 
// Dtwin     // Total search time window (in sec) centered at GRB trigger
// Dawin     // Angular search circle (<0 is decl. band) in degrees or sigma around (above/below) GRB position
// Datype    // Type of angular window specification (0=in degrees 1=in units of combined GRB/track sigma) 
// Nbkg      // Mean number of counts per bin for auto-binning
// Tbint90   // Time bin size in units of average T90 (0 : Time bin size determined by Tbin) 
// Tbin      // Time bin size in seconds (0=variable bins  <0 will result in a mean Nbkg counts/bin)
// VarTbin   // Size (in sec) of the first time bin in case of variable time bins
// Abin      // Angular bin size in degrees (<0 will result in a mean fNbkg counts/bin)
// Freq      // Use frequentist's approximation (1) or exact Bayesian expression (0)
// Npsi      // Number of psi entries for bkg psi-value distributions (<0 : time shuffling)
// Usetott   // Use the observed tott number of entries in case of time shuffling 
// Grbpos    // Use the original burst locations (1) or random ones (0) for bkg studies
// Nrandom   // (Maximum) number of randomised configurations for direct psi P-value determination (0 means n=1e19)
// Ncut      // Number of obtained psi>psi0 entries at which randomisations will be terminated (0= no early termination)
//
// The meaning of the return value is the following :
// 1 : Parameter has been set
// 0 : Parameter has not been set

 Int_t ret=0;

 if (name=="Nmax")
 {
  fNmax=value;
  ret=1;
 }
 if (name=="Declmin")
 {
  fDeclmin=value;
  ret=1;
 }
 if (name=="Declmax")
 {
  fDeclmax=value;
  ret=1;
 }
 if (name=="T90min")
 {
  fT90min=value;
  ret=1;
 }
 if (name=="T90max")
 {
  fT90max=value;
  ret=1;
 }
 if (name=="Zmin")
 {
  fZmin=value;
  ret=1;
 }
 if (name=="Zmax")
 {
  fZmax=value;
  ret=1;
 }
 if (name=="Sigmagrb")
 {
  fSigmagrb=value;
  ret=1;
 }
 if (name=="Maxsigma")
 {
  fMaxsigma=value;
  ret=1;
 }
 if (name=="Grbnu")
 {
  fGrbnu=value;
  ret=1;
 }
 if (name=="Avgrbz")
 {
  fAvgrbz=value;
  ret=1;
 }
 if (name=="Avgrbt90")
 {
  fAvgrbt90=value;
  ret=1;
 }
 if (name=="Inburst")
 {
  fInburst=value;
  ret=1;
 }
 if (name=="Dtnu")
 {
  fDtnu=value;
  ret=1;
 }
 if (name=="Dtnus")
 {
  fDtnus=value;
  ret=1;
 }
 if (name=="Angres")
 {
  fAngres=value;
  ret=1;
 }
 if (name=="Timres")
 {
  fTimres=value;
  ret=1;
 }
 if (name=="Bkgrate")
 {
  fBkgrate=value;
  ret=1;
 }
 if (name=="Nevtdt")
 {
  fNevtdt=value;
  ret=1;
 }
 if (name=="Dtwin")
 {
  fDtwin=value;
  ret=1;
 }
 if (name=="Dawin")
 {
  fDawin=value;
  ret=1;
 }
 if (name=="Datype")
 {
  fDatype=value;
  ret=1;
 }
 if (name=="Nbkg")
 {
  fNbkg=value;
  ret=1;
 }
 if (name=="Tbint90")
 {
  fTbint90=value;
  ret=1;
 }
 if (name=="Tbin")
 {
  fTbin=value;
  ret=1;
 }
 if (name=="VarTbin")
 {
  fVarTbin=value;
  ret=1;
 }
 if (name=="Abin")
 {
  fAbin=value;
  ret=1;
 }
 if (name=="Freq")
 {
  fFreq=value;
  ret=1;
 }
 if (name=="Npsi")
 {
  fNpsi=value;
  ret=1;
 }
 if (name=="Usetott")
 {
  fUsetott=value;
  ret=1;
 }
 if (name=="Grbpos")
 {
  fGrbpos=value;
  ret=1;
 }
 if (name=="Nrandom")
 {
  fNrandom=value;
  ret=1;
 }
 if (name=="Ncut")
 {
  fNcut=value;
  ret=1;
 }

 return ret;
}
///////////////////////////////////////////////////////////////////////////
void NcGRB::Exec(Option_t* opt)
{
// Perform the analysis

/******
@@@ To be activated when real data treatment is implemented
 TString jobname=opt;
 NcJob* parent=(NcJob*)(gROOT->GetListOfTasks()->FindObject(jobname.Data()));

 if (!parent) return;

 fEvt=(NcEvent*)parent->GetObject("NcEvent");
 if (!fEvt) return;

 // Only process accepted events
 NcDevice* seldev=(NcDevice*)fEvt->GetDevice("NcEventSelector");
 if (seldev)
 {
  if (seldev->GetSignal("Select") < 0.1) return;
 }
******/

 if (fT90min<0) fT90min=0;

 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 Float_t t90,z;
 TString name;

 ////////////////////////////////////////////////
 // Some Burst statistics from the loaded data //
 ////////////////////////////////////////////////

 fNgrbs=GetNsignals(0);

 // Creation of the GRB position uncertainty histo
 TH1F* hsigmagrb=new TH1F("hsigmagrb","GRB position uncertainty",450,0,90);
 fHistos.Add(hsigmagrb);
 hsigmagrb->GetXaxis()->SetTitle("GRB position uncertainty (sigma in degrees)");
 hsigmagrb->GetYaxis()->SetTitle("Counts");

 // Creation of the combined GRB position and track resolution uncertainty histo
 TH1F* htotsigma=new TH1F("htotsigma","Combined GRB position and track resolution uncertainty",450,0,90);
 fHistos.Add(htotsigma);
 htotsigma->GetXaxis()->SetTitle("Combined GRB position and track resolution uncertainty (sigma in degrees)");
 htotsigma->GetYaxis()->SetTitle("Counts");

 NcSignal* sx=0;
 NcSample zsample;
 zsample.SetStoreMode();
 NcSample t90sample;
 t90sample.SetStoreMode();
 Int_t nsig=GetNsignals(0,1);
 for (Int_t i=1; i<=nsig; i++)
 {
  sx=GetSignal(i,0);

  if (!sx) continue;

  hsigmagrb->Fill(sx->GetSignal("sigmagrb"));
  htotsigma->Fill(sx->GetSignal("totsigma"));

  if (fAvgrbz<0) zsample.Enter(sx->GetSignal("z"));
  if (fAvgrbt90<0) t90sample.Enter(sx->GetSignal("t90"));
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
 fMup=fBkgrate*fDtwin;

 // Mean number of upgoing bkg muons per day
 fNmupday=fBkgrate*24.*3600.;

 Float_t danglow=0;    // Lower value (in degrees) of angular difference histo
 Float_t dangup=fDawin; // Upper value (in degrees) of angular difference histo
 if (fDatype) dangup=fDawin*fabs(fMaxtotsigma);
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
   if (fTbint90) fTbin=fTbint90*fabs(fAvgrbt90); 
   ntbins=int(fDtwin/fTbin);
  }
  else // Automatic time binning to get the specified maximal bkg counts per bin
  {
   ntbins=int(fMup*float(fNgrbs)/fNbkg);
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
 if (fAbin<0) nabins=int(((dangup-danglow)/180.)*fMup*float(fNgrbs)/fNbkg);

 // Binning for the cos(opening angle) histo
 Float_t upcos=cos(danglow*pi/180.);
 Float_t lowcos=cos(dangup*pi/180.);
 Int_t nabins2=int((upcos-lowcos)/(1.-cos(fAbin*pi/180.)));
 if (fAbin<0) nabins2=int(((upcos-lowcos)/2.)*fMup*float(fNgrbs)/fNbkg);
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
 fHistos.Add(bkgtfine);

 title="t of all mu-up in twindow";
 title+=";Upgoing #mu arrival time (in sec) w.r.t. GRB #gamma trigger;Counts per bin of size %-10.3g";
 s=title.Format(title.Data(),tbinfine);
 TH1F* tottfine=new TH1F("tottfine",s.Data(),ntbinsfine,-fDtwin/2.,fDtwin/2.);
 fHistos.Add(tottfine);

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
 fHistos.Add(bkgt);
 fHistos.Add(tott);
 fHistos.Add(bkg2);
 fHistos.Add(tot2);

 // The opening angle histo
 TH1F* bkga=new TH1F("bkga","dang of bkg mu-up in twindow",nabins,danglow,dangup);
 TH1F* tota=new TH1F("tota","dang of all mu-up in twindow",nabins,danglow,dangup);
 fHistos.Add(bkga);
 fHistos.Add(tota);

 // The cos(opening angle) histo
 TH1F* bkgcosa=new TH1F("bkgcosa","cos(dang) of bkg mu-up in twindow",nabins2,lowcos,upcos);
 TH1F* totcosa=new TH1F("totcosa","cos(dang) of all mu-up in twindow",nabins2,lowcos,upcos);
 fHistos.Add(bkgcosa);
 fHistos.Add(totcosa);

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

 Float_t t90grb=0;
 Double_t zgrb=0;
//@@@ TString grbname;
 Float_t sigmagrb=0;
 Float_t totsigma=0;
 NcPosition rgrb;
 Int_t nmup;
 Double_t thetagrb,phigrb;
 Float_t thetamu,phimu,cost;
 Float_t dt=0;
 NcPosition rgrb2; // Actual GRB position from which the neutrinos/muons arrive
 NcPosition rmu;
 Float_t dang;
 Float_t ranlow,ranup;
 Float_t thlow,thup;
 Int_t nmugrb=0;
 NcTimestamp* tx=0;
 Float_t solidangle=0; // Total stacked solid angle 

 // Dummy invokation to ensure proper initialisation of the random number generator
 RandomPosition(rmu,90,180,0,360);

 // Obtain the (fictative) GRB space-time positions in the declination acceptance
 for (Int_t igrb=0; igrb<fNgrbs; igrb++)
 {
  sx=GetSignal(igrb+1);

  if (!sx) continue;

  tx=sx->GetTimestamp();
  t90grb=sx->GetSignal("t90");
  sigmagrb=sx->GetSignal("sigmagrb");
  totsigma=sx->GetSignal("totsigma");
  GetSignal(zgrb,thetagrb,"deg",phigrb,"deg","loc",tx,igrb+1);
  rgrb.SetPosition(zgrb,thetagrb,phigrb,"sph","deg");

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

  solidangle+=GetSolidAngle(thlow,thup,"deg",0,360,"deg");    

  // Generate the upgoing bkg muons in the search time window 
  // for both this GRB angular cone and the corresponding "opposite RA" bkg patch
  for (Int_t bkgpatch=0; bkgpatch<=1; bkgpatch++)
  {
   nmup=int(fRan->Poisson(fMup));
   for (Int_t imup=0; imup<nmup; imup++)
   {
    ranlow=-fDtwin/2.;
    ranup=fDtwin/2.;
    dt=fRan->Uniform(ranlow,ranup);
    // Smear the time difference with the Gaussian time resolution
    if (fTimres>0) dt=fRan->Gauss(dt,fTimres); //@@@ Is this needed ?
    RandomPosition(rmu,90,180,0,360);
    // Smear the direction of the upgoing bkg muon according to  the detector resolution
    SmearPosition(rmu,fAngres); //@@@ Is this needed
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
  SmearPosition(rgrb2,sigmagrb);

  nmup=int(fabs(fGrbnu));
  if (!nmup && fRan->Uniform()<fabs(fGrbnu)) nmup=1;
  for (Int_t imup2=0; imup2<nmup; imup2++)
  {
   nmugrb++;
   if (!fInburst) // Neutrino and gamma production decoupled
   {
    if (fDtnus<0) // Sigma in units of T90
    {
     dt=fRan->Gauss(fDtnu,fabs(fDtnus)*t90grb);
    }
    else
    {
     dt=fRan->Gauss(fDtnu,fDtnus);
    }
    dt=dt*(zgrb+1.);
   }
   else // Coupled neutrino and gamma production
   {
    if (fDtnus<0) // Sigma in units of T90
    {
     dt=fRan->Gauss(fDtnu*t90grb,fabs(fDtnus)*t90grb);
    }
    else
    {
     dt=fRan->Gauss(fDtnu*t90grb,fDtnus);
    }
   }
   if (fTimres>0) dt=fRan->Gauss(dt,fTimres);
   rmu.Load(rgrb2);
   // Smear the direction of the upgoing GRB muon according to the detector resolution
   if (fAngres>0) SmearPosition(rmu,fAngres);

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
 if (fGrbnu<0) Compensate(nmugrb);


 /////////////////////////////////////////////////////////////////
 // Creation of the inter-muon time histogram                   //
 // to reflect possible concentration(s) of muon arrival times. //
 // The fine time histos are used for the basic input.          //
 /////////////////////////////////////////////////////////////////

 NcMath math;

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

 ////////////////////////////////////////////////////////////////////////////////
 // Statistical evaluation of the generated signal and background observations //
 //                                                                            //
 // Determination of the Bayesian psi value for the time and angular histos    //
 // under the assumption that there is no GRB signal.                          //
 // This corresponds to searching out the Bernoulli class B_m                  //
 // with m=nbins of the histogram.                                             //
 // An orthodox chi-squared analysis is also performed.                        //
 ////////////////////////////////////////////////////////////////////////////////

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

/*****************
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
     sx=GetSignal(igrb+1);

     if (!sx) continue;

     tx=sx->GetTimestamp();
     totsigma=sx->GetSignal("totsigma");
     GetSignal(zgrb,thetagrb,"deg",phigrb,"deg","loc",tx,igrb+1);
     rgrb.SetPosition(zgrb,thetagrb,phigrb,"sph","deg");
    }
    else
    {
     thlow=fDeclmin+90.;
     thup=fDeclmax+90.;
     if (thup>180) thup=180;
     zgrb=hz->GetRandom();
     rgrb.SetPosition(zgrb,0.,0.,"sph","deg");
     RandomPosition(rgrb,thlow,thup,0,360);
     thetagrb=rgrb.GetX(2,"sph","deg");
     phigrb=rgrb.GetX(3,"sph","deg");

     // Determine the combined GRB position and track resolution uncertainty.
     totsigma=sigmagrb*sigmagrb+fAngres*fAngres;
     totsigma=sqrt(totsigma);
    }

    // Generate the upgoing bkg muons in the search window for this GRB
    nmup=int(fRan->Poisson(fMup));
    for (Int_t jmup=0; jmup<nmup; jmup++)
    {
     ranlow=-fDtwin/2.;
     ranup=fDtwin/2.;
     dt=fRan->Uniform(ranlow,ranup);
     // Smear the time difference with the Gaussian time resolution
     if (fTimres>0) dt=fRan->Gauss(dt,fTimres);
     RandomPosition(rmu,90,180,0,360);
     // Smear the direction of the upgoing bkg muon according to the detector resolution
     SmearPosition(rmu,fAngres);
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
     dt=fRan->Uniform(ranlow,ranup);
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
     dt=fRan->Uniform(ranlow,ranup);
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

 if (fNpsi)
 {
  cout << " Background studies : hpsibkgt hpsibkga hpsit hpsia" << endl;
  cout << " and all hpsibkgt and hpsibkga histos in the TObjArrays bkgthists and bkgahists" << endl;
 }
****************/
}
///////////////////////////////////////////////////////////////////////////
void NcGRB::PrintSettings() const
{
// Print all the parameter settings

 cout << " ========================= NcGRB user provided settings ===============================" << endl;
 if (fNmax<0)
 {
  cout << " No limitation has been put on the number of bursts to be accepted for analysis." << endl;
 }
 else
 {
  cout << " Maximal number of bursts to be accepted for analysis : " << fNmax << endl;
 }
 cout << " Declination interval (in degrees) for GRB position acceptance : [" << fDeclmin << "," << fDeclmax << "]" << endl;
 cout << " Duration interval (t90 in sec) for GRB acceptance : [" << fT90min << "," << fT90max << "]" << endl;
 cout << " Redshift interval for GRB acceptance : [" << fabs(fZmin) << "," << fZmax << "]" << endl;
 if (fZmin<0) cout << " Random redshift values taken from z-distribution in case of unknown redshift" << endl;
 if (fSigmagrb>=0) cout << " Fixed GRB position uncertainty (sigma in degrees) : " << fSigmagrb << endl;
 cout << " Maximal combined GRB position and track angular uncertainty (sigma in degrees) for acceptance : " << fMaxsigma << endl;
 if (fAvgrbz>=0) cout << " User defined average GRB redshift : " << fAvgrbz << endl;
 if (fAvgrbt90>=0) cout << " User defined average GRB T90 duration : " << fAvgrbt90 << endl;
 if (!fInburst)
 {
  cout << " Neutrino production was assumed to be NOT coupled to the gamma flare duration" << endl;
  cout << " Mean decoupled time difference (in sec) between GRB gammas and nus : " << fDtnu << endl;
 }
 else
 {
  cout << " Neutrino production was assumed to be coupled to the gamma flare duration" << endl;
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
 if (fGrbnu<0)
 {
  cout << " Number of generated neutrinos per GRB : " << fabs(fGrbnu) << " without statistical fluctuations" << endl; 
 }
 else
 {
  cout << " Maximum number of generated neutrinos per GRB : " << fGrbnu << endl;
  cout << " The actual number of neutrinos may be less due to statistical fluctuations" << endl;
 }
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
 if (!fTbin) cout << " Variable time binning with as size (in sec) for the first time : " << fVarTbin << endl;
 if (fTbin>0)
 {
  if (fTbint90)
  {
   cout << " Time bin size in average T90 units : " << fTbint90;
   if (fAvgrbt90>0 || fNgrbs>0) cout << " (=" << fTbin << " sec)";
   cout << endl;
  }
  else
  {
   cout << " Time bin size in seconds : " << fTbin << endl;
  }
 }
 if (fAbin<0)
 {
  cout << " Automatic angular binning with as mean number of bkg counts per bin : " << fNbkg << endl;
 }
 else
 {
  cout << " Angular bin size in degrees : " << fAbin << endl;
 }
 if (fFreq) cout << " Frequentist's approximation used for psi determination" << endl;
 if (fNpsi)
 {
  cout << " Number of psi entries for bkg psi-value distributions : " << abs(fNpsi);
  if (fNpsi<0)
  {
   cout << " by means of time shuffling";
   if (fUsetott) cout << " of the observed tott number of entries";
  }
  cout << endl;
 }
 if (fGrbpos)
 {
  cout << " Original burst locations are used for bkg studies" << endl;
 }
 else
 {
  cout << " Random burst locations are used for bkg studies" << endl;
 }
 if (fNrandom>0)
 {
  if (!fNcut)
  {
   cout << " Number of randomised configurations for direct psi P-value determination : " << fNrandom << endl;
  }
  else
  {
   cout << " Maximum number of randomised configurations for direct psi P-value determination : " << fNrandom << endl;
   cout << " Randomisation will be terminated when " << fNcut << " entries with psi>psi0 have been obtained." << endl;
  }
 }
 cout << " ======================================================================================" << endl;
 cout << endl;
 if (fNgrbs>0)
 {
  cout << " ============================== Derived parameters ====================================" << endl;
  cout << " Mean number of upgoing bkg muons in the time window in the detector : " << fMup << endl;
  cout << " Mean number of upgoing bkg muons per day in the detector : " << fNmupday << endl;
  cout << " Number of bursts accepted for analysis : " << fNgrbs << endl;
  cout << " Median redshift from the data sample : " << fabs(fAvgrbz) << endl;
  cout << " Median T90 duration from the data sample : " << fabs(fAvgrbt90) << endl;
  cout << " ======================================================================================" << endl;
  cout << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcGRB::LoadBurstData(TString file,TString tree,Int_t date1,Int_t date2,Int_t nmax)
{
// Load observed burst data, e.g. GRB data from GCN notices as available from
// https://icecube.wisc.edu/~grbweb_public.
// The input data has to be provided via a ROOT Tree which contains at least
// the specified Branch names (see below) filled with data of the indicated type.
//
//  Name        Type      Description
// ----------------------------------
// "date"       Int_t     Observation date as yyyymmdd
// "ra"         Float_t   Right ascension (J2000) in decimal degrees
// "decl"       Float_t   Declination (J2000) in decimal degrees
// "sigpos      Float_t   The error on the burst position in decimal degrees
// "t90"        Float_t   The T90 burst duration in seconds
// "mjdtrig"    Double_t  The (fractional) MJD of the burst trigger
// "mjdt90start Double_t  The (fractional) MJD of the T90 start time
// "t100"       Float_t   The T100 burst duration in seconds
// "fluence"    Float_t   The observed fluence of the burst in erg/cm^2
// "z"          Float_t   The redshift of the burst
//
// Input arguments :
// -----------------
// file  : Name of the input file containing the ROOT Tree (wildcards are allowed)
// tree  : Name of the Tree containing the data
// date1 : The date (yyyymmdd) of the start of the observation period [date1,date2] (0=No restriction)
// date2 : The date (yyyymmdd) of the end of the observation period [date1,date2] (0=No restriction)
// nmax  : Maximum number of bursts to be accepted from this input file (<0 : no limitation)
//
// The default values are date1=0, date2=0 and nmax=-1.
//
// Note : This memberfunction make be invoked several times to read different files
//        to accumulate data.

 // The TTree containing the burst data
 TChain gcn(tree.Data());
 gcn.Add(file.Data());

 Int_t date;
 Float_t ra,decl,sigmapos,t90,t100,fluence,z;
 Double_t mjdtrig,mjdt90start;

 // The variables from the Tree
 gcn.SetBranchAddress("date",&date);
 gcn.SetBranchAddress("ra",&ra);
 gcn.SetBranchAddress("decl",&decl);
 gcn.SetBranchAddress("sigmapos",&sigmapos);
 gcn.SetBranchAddress("t90",&t90);
 gcn.SetBranchAddress("mjdtrig",&mjdtrig);
 gcn.SetBranchAddress("mjdt90start",&mjdt90start);
 gcn.SetBranchAddress("t100",&t100);
 gcn.SetBranchAddress("fluence",&fluence);
 gcn.SetBranchAddress("z",&z);

 // Get access to a redshift distribution to draw randomly redshifts if needed
 TH1* zdist=0;
 if (fZmin<0)
 {
  zdist=(TH1*)fHistos.FindObject("hz");
  if (!zdist)
  {
   cout << " *" << ClassName() << "::LoadBurstData* Archival observed redshift distribution not found." << endl;
   cout << " A Landau fit from Swift GRB redshift data will be used to provide missing z values." << endl;

   zdist=(TH1*)fHistos.FindObject("hzfit");
   if (!zdist)
   {
    TF1 f("f","59.54*TMath::Landau(x,1.092,0.5203)");
    f.SetRange(0,10);
    f.SetNpx(10000);
    TH1* hf=f.GetHistogram();
    zdist=(TH1*)hf->Clone();
    zdist->SetNameTitle("hzfit","Landau fit for Swift GRB z data");
    zdist->GetXaxis()->SetTitle("GRB redshift");
    zdist->GetYaxis()->SetTitle("Counts");
    fHistos.Add(zdist);
   }
  }
 }

 Float_t sigmagrb=fSigmagrb;
 Float_t totsigma=0;
 fNgrbs=GetNsignals(0);
 Float_t t90grb=0;
 Double_t zgrb=0;
 Int_t idate=0;
 TString grbname;
 NcTimestamp ts;
 NcSignal* sx=0;
 Int_t ngcn=0;
 for (Int_t ient=0; ient<gcn.GetEntries(); ient++)
 {
  if (nmax>=0 && ngcn>=nmax) break;
  if (fNmax>=0 && (fNgrbs+ngcn)>=fNmax) break;

  gcn.GetEntry(ient);

  if (date1 && date<date1) continue;
  if (date2 && date>date2) continue;

  if (mjdtrig<0 || mjdt90start<0) continue;

  // Use the GCN 1 sigma position uncertainty (in degrees) if requested
  if (fSigmagrb<0) sigmagrb=fabs(sigmapos);

  // Determine the combined GRB position and track resolution uncertainty.
  totsigma=sigmagrb*sigmagrb+fAngres*fAngres;
  totsigma=sqrt(totsigma);
   
  if (decl<fDeclmin || decl>fDeclmax || totsigma>fMaxsigma) continue;

  t90grb=t90;
  if (t90grb<=0) t90grb=t100;

  if (t90grb<fT90min || t90grb>fT90max) continue;

  zgrb=z;
  if (fZmin<0 && zgrb<0 && zdist) zgrb=zdist->GetRandom();

  if (zgrb<fabs(fZmin) || zgrb>fZmax) continue;

  idate=date%1000000;
  grbname="GRB";
  grbname+=idate;
  ts.SetMJD(mjdtrig);
  sx=SetSignal(1,ra,"deg",decl,"deg","equ",&ts,-1,"J",grbname);

  if (!sx) continue;

  ngcn++;

  sx->AddNamedSlot("t90");
  sx->SetSignal(t90grb,"t90");
  sx->AddNamedSlot("sigmagrb");
  sx->SetSignal(sigmagrb,"sigmagrb");
  sx->AddNamedSlot("totsigma");
  sx->SetSignal(totsigma,"totsigma");
  sx->AddNamedSlot("fluence");
  sx->SetSignal(fluence,"fluence");
  sx->AddNamedSlot("z");
  sx->SetSignal(zgrb,"z");

  if (totsigma>fMaxtotsigma) fMaxtotsigma=totsigma;
 }

 cout << "*NcGRB::LoadBurstData* " << ngcn << " bursts were stored from Tree:" << tree << " of file(s):" << file << endl;
}
///////////////////////////////////////////////////////////////////////////
void NcGRB::GenBurstData(Int_t n)
{
// Generate fictative burst data for "n" bursts.

 // Get access to a redshift distribution to draw randomly redshifts
 TH1* zdist=(TH1*)fHistos.FindObject("hz");
 if (!zdist)
 {
  cout << " *" << ClassName() << "::GenBurstData* Archival observed redshift distribution not found." << endl;
  cout << " A Landau fit from Swift GRB redshift data will be used to provide random z values." << endl;

  zdist=(TH1*)fHistos.FindObject("hzfit");
  if (!zdist)
  { 
   TF1 fz("fz","59.54*TMath::Landau(x,1.092,0.5203)");
   fz.SetRange(0,10);
   fz.SetNpx(10000);
   TH1* hfz=fz.GetHistogram();
   zdist=(TH1*)hfz->Clone();
   zdist->SetNameTitle("hzfit","Landau fit for Swift GRB z data");
   zdist->GetXaxis()->SetTitle("GRB redshift");
   zdist->GetYaxis()->SetTitle("Counts");
   fHistos.Add(zdist);
  }
 }

 // Get access to a redshift distribution to draw randomly redshifts
 TH1* t90dist=(TH1*)fHistos.FindObject("ht90");
 if (!t90dist) t90dist=(TH1*)fHistos.FindObject("ht90fit");
 if (!t90dist)
 {
  cout << " *" << ClassName() << "::GenBurstData* Observational T90 distribution not found." << endl;
  cout << " A double Gaussian fit from Fermi GRB T90 data will be used to provide random T90 values." << endl;
  cout << endl;

  t90dist=(TH1*)fHistos.FindObject("ht90fit");
  if (!t90dist)
  {
   TF1 ft("ft","44.39*TMath::Gaus(x,-0.131,0.481)+193.8*TMath::Gaus(x,1.447,0.4752)");
   ft.SetRange(-5,5);
   ft.SetNpx(10000);
   TH1* hft=ft.GetHistogram();
   t90dist=(TH1*)hft->Clone();
   t90dist->SetNameTitle("ht90fit","Double Gauss fit for Fermi t90 data");
   t90dist->GetXaxis()->SetTitle("GRB duration ^{10}log(T90) in sec.");
   t90dist->GetYaxis()->SetTitle("Counts");
   fHistos.Add(t90dist);
  } 
 }

 Float_t thlow=fDeclmin+90.;
 Float_t thup=fDeclmax+90.;
 if (thup>180) thup=180;

 NcSignal* sx=0;
 NcPosition rgrb;
 Float_t t90grb=0;
 Double_t zgrb=0;
 TString grbname;
 Float_t sigmagrb=fabs(fSigmagrb);
 Float_t totsigma=0;
 Double_t thetagrb,phigrb;
 Int_t ngen=0;
 fNgrbs=GetNsignals(0);

 for (Int_t igrb=1; igrb<=n; igrb++)
 {
  if (fNmax>=0 && (fNgrbs+ngen)>=fNmax) break;

  zgrb=-1;
  while (zgrb<fabs(fZmin) || zgrb>fZmax)
  {
   zgrb=zdist->GetRandom();
  }
  rgrb.SetPosition(zgrb,0,0,"sph","deg");
  RandomPosition(rgrb,thlow,thup,0,360);
  thetagrb=rgrb.GetX(2,"sph","deg");
  phigrb=rgrb.GetX(3,"sph","deg");

  t90grb=-1;
  while (t90grb<fT90min || t90grb>fT90max)
  {
   t90grb=t90dist->GetRandom();
   t90grb=pow(float(10),t90grb);
  }

  // Determine the combined GRB position and track resolution uncertainty.
  totsigma=sigmagrb*sigmagrb+fAngres*fAngres;
  totsigma=sqrt(totsigma);
   
  if (totsigma>fMaxsigma) continue;

  grbname="Random-GRB";
  grbname+=igrb;
  sx=SetSignal(zgrb,thetagrb,"deg",phigrb,"deg","loc",0,-1,"M",grbname);

  if (!sx) continue;

  ngen++;

  sx->AddNamedSlot("t90");
  sx->SetSignal(t90grb,"t90");
  sx->AddNamedSlot("sigmagrb");
  sx->SetSignal(sigmagrb,"sigmagrb");
  sx->AddNamedSlot("totsigma");
  sx->SetSignal(totsigma,"totsigma");
  sx->AddNamedSlot("z");
  sx->SetSignal(zgrb,"z");

  if (totsigma>fMaxtotsigma) fMaxtotsigma=totsigma;
 }

 cout << "*NcGRB::GenBurstData* " << ngen << " generated bursts were stored." << endl;
}
///////////////////////////////////////////////////////////////////////////
void NcGRB::MakeZdist(TString file,TString tree,TString branch,Int_t nb,Float_t zmin,Float_t zmax)
{
// Read observed archival redshift data and create the corresponding distribution.
// Also the corresponding distribution of the Proper Distance at the time of
// observation (called the Physical Distance) will be created.
// If this memberfunction is invoked before LoadBurstData() or GenBurstData(),
// the plain observed redshift distribution will be used to draw random z values
// (if requested) for the bursts without redshift information.
//
// The input data has to be provided via a ROOT Tree which contains at least
// the specified Branch name (see below) filled with data of type Float_t.
//
// Input arguments :
// -----------------
// file   : Name of the input file(s) containing the ROOT Tree (wildcards are allowed)
// tree   : Name of the Tree containing the data
// branch : Name of the Branch containing the redshift data in Float_t format
// nb     : Number of bins for the z-distribution
// zmin   : Minimal z-value
// zmax   : Maximal z-value
//
// The default values are nb=200, zmin=0 and zmax=20.
//
// Note : This memberfunction may be invoked several times to read different files
//        to accumulate data.

 Float_t z=0;

 // The Tree containing the archival data
 TChain data(tree.Data());
 data.Add(file.Data());

 Int_t nen=data.GetEntries();
 TBranch* bx=data.GetBranch(branch.Data());

 if (!nen || !bx) return;

 // The Tree branch with the redshift data
 data.SetBranchAddress(branch.Data(),&z);

 // Create new distributions in case a redshift distribution is not yet present
 TH1* zdist=(TH1*)fHistos.FindObject("hz");
 if (!zdist)
 {
  // Creation of the archival GRB redshift histogram
  TH1F* hz=new TH1F("hz","Archival data of observed GRB redshifts",nb,zmin,zmax);
  fHistos.Add(hz);
  hz->GetXaxis()->SetTitle("GRB redshift");
  hz->GetYaxis()->SetTitle("Counts");

  // Creation of the corresponding physical distance histo
  Float_t dmin=GetPhysicalDistance(zmin);
  Float_t dmax=GetPhysicalDistance(zmax);
  TH1F* hd=new TH1F("hd","GRB distances derived from the archival redshift data",nb,dmin,dmax);
  fHistos.Add(hd);
  hd->GetXaxis()->SetTitle("GRB physical distance in Mpc");
  hd->GetYaxis()->SetTitle("Counts");
 }

 // Get pointers to the relevant histograms 
 TH1* hz=(TH1*)fHistos.FindObject("hz");
 TH1* hd=(TH1*)fHistos.FindObject("hd");

 Int_t nz=0;
 Float_t d=0;
 for (Int_t ien=0; ien<nen; ien++)
 {
  data.GetEntry(ien);
  if (z<zmin || z>zmax) continue;

  hz->Fill(z);
  nz++;

  d=GetPhysicalDistance(z);
  hd->Fill(d);
 }

 cout << "*NcGRB::MakeZdist* " << nz << " archival z-values have been obtained from Branch:" << branch
      << " of Tree:" << tree << " in file(s):" << file << endl;
}
///////////////////////////////////////////////////////////////////////////
void NcGRB::MakeT90dist(TString file,TString tree,TString branch,Int_t nb,Float_t xmin,Float_t xmax)
{
// Read observed archival T90 data and create a log10(T90) distribution.
// If this memberfunction is invoked before LoadBurstData() or GenBurstData(),
// the resulting log10(T90) distribution will be used to draw random T90 values
// (if requested) for the bursts without T90 information.
//
// The input data has to be provided via a ROOT Tree which contains at least
// the specified Branch name (see below) filled with data of type Float_t.
//
// Input arguments :
// -----------------
// file   : Name of the input file containing the ROOT Tree (wildcards are allowed)
// tree   : Name of the Tree containing the data
// branch : Name of the Branch containing the T90 data in Float_t format
// nb     : Number of bins for the T90 distribution
// xmin   : Minimal value for log10(T90)
// xmax   : Maximal value for log10(T90)
//
// The default values are nb=50, xmin=-5 and xmax=5.
//
// Note : This memberfunction may be invoked several times to read different files
//        to accumulate data.

 Float_t t90;

 // The Tree containing the burst data
 TChain data(tree.Data());
 data.Add(file.Data());

 Int_t nen=data.GetEntries();
 TBranch* bx=data.GetBranch(branch.Data());

 if (!nen || !bx) return;

 // The Tree branch with the T90 data
 data.SetBranchAddress(branch.Data(),&t90);

 // Create a new distribution in case a T90 distribution is not yet present
 TH1* t90dist=(TH1*)fHistos.FindObject("ht90");
 if (!t90dist)
 {
  // Creation of observed GRB t90 duration histo
  TH1F* ht90=new TH1F("ht90","Archival data of observed GRB durations",50,-5,5);
  fHistos.Add(ht90);
  ht90->GetXaxis()->SetTitle("GRB duration ^{10}log(T90) in sec.");
  ht90->GetYaxis()->SetTitle("Counts");
 }

 // Get pointer to the relevant histogram 
 TH1* ht90=(TH1*)fHistos.FindObject("ht90");

 Int_t nt90=0;
 for (Int_t ien=0; ien<nen; ien++)
 {
  data.GetEntry(ien);
  if (t90>0)
  {
   ht90->Fill(log10(t90));
   nt90++;
  }
 }

 cout << "*NcGRB::MakeT90dist* " << nt90 << " archival T90 values have been obtained from Branch:" << branch
      << " of Tree:" << tree << " in file(s):" << file << endl;
}
///////////////////////////////////////////////////////////////////////////
void NcGRB::ListHistograms() const
{
// Provide a list of all the stored histograms

 Int_t nh=fHistos.GetEntries();
 cout << endl;
 cout << " =============== The following " << nh << " histograms have been generated ===============" << endl;
 for (Int_t ih=0; ih<nh; ih++)
 {
  TObject* hx=fHistos.At(ih);
  if (!hx) continue;
  cout << " " << hx->GetName() << " : " << hx->GetTitle() << endl;
 }
 cout << " ===============================================================================" << endl;
}
///////////////////////////////////////////////////////////////////////////
void NcGRB::WriteHistograms(TString filename)
{
// Write all the generated histograms to a ROOT file with the specified filename.

 // The output file for the produced histograms
 TFile fout(filename.Data(),"RECREATE","NcGRB analysis results");

 // Write all the histos to the output file
 Int_t nh=fHistos.GetEntries();
 for (Int_t ih=0; ih<nh; ih++)
 {
  TObject* hx=fHistos.At(ih);
  if (!hx) continue;
  hx->Write();
 }

 fout.Write();

 cout << endl;
 cout << " *" << ClassName() << "::WriteHistograms* All generated histograms have been written to file " << filename << endl;
 ListHistograms();
}
///////////////////////////////////////////////////////////////////////////
void NcGRB::Compensate(Int_t& nmugrb)
{
// Compensate statistical underfluctuation in the number of GRB muons.

 Int_t nmup=int(fabs(fGrbnu)*float(fNgrbs));
 Int_t jgrb=0;
 NcSignal* sx=0;
 NcTimestamp* tx=0;
 Float_t t90grb=0;
 Float_t sigmagrb=0;
 Float_t totsigma=0;
 NcPosition rgrb;
 NcPosition rgrb2;
 Float_t dt=0;
 Double_t zgrb=0;
 Double_t thetagrb=0;
 Double_t phigrb=0;
 Float_t dang=0;
 NcPosition rmu;

 TH1* tottfine=(TH1*)fHistos.FindObject("tottfine");
 TH1* tott=(TH1*)fHistos.FindObject("tott");
 TH1* tota=(TH1*)fHistos.FindObject("tota");
 TH1* totcosa=(TH1*)fHistos.FindObject("totcosa");
 TH2* tot2=(TH2*)fHistos.FindObject("tot2");

 Double_t pi=acos(-1.);

 while (nmugrb<nmup)
 {
  // Pick randomly one of the stored GRBs
   jgrb=int(fRan->Uniform(0.,float(fNgrbs)));
   if (jgrb==0) jgrb=1;
   sx=GetSignal(jgrb);

   if (!sx) continue;

   tx=sx->GetTimestamp();
   t90grb=sx->GetSignal("t90");
   GetSignal(zgrb,thetagrb,"deg",phigrb,"deg","loc",tx,jgrb);
   rgrb.SetPosition(zgrb,thetagrb,phigrb,"sph","deg");
   sigmagrb=sx->GetSignal("sigmagrb");
   totsigma=sx->GetSignal("totsigma");

   // Obtain actual GRB position
   rgrb2.Load(rgrb);
   SmearPosition(rgrb2,sigmagrb); //@@@

   nmugrb++;

   if (!fInburst) // Neutrino and gamma production decoupled
   {
    if (fDtnus<0) // Sigma in units of T90
    {
     dt=fRan->Gauss(fDtnu,fabs(fDtnus)*t90grb);
    }
    else
    {
     dt=fRan->Gauss(fDtnu,fDtnus);
    }
    dt=dt*(zgrb+1.);
   }
   else // Coupled neutrino and gamma production
   {
    if (fDtnus<0) // Sigma in units of T90
    {
     dt=fRan->Gauss(fDtnu*t90grb,fabs(fDtnus)*t90grb);
    }
    else
    {
     dt=fRan->Gauss(fDtnu*t90grb,fDtnus);
    }
   }
   if (fTimres>0) dt=fRan->Gauss(dt,fTimres);
   rmu.Load(rgrb2);
   // Smear the direction of the upgoing GRB muon according to the detector resolution
   if (fAngres>0) SmearPosition(rmu,fAngres);

   // Determine angular difference w.r.t. the presumed GRB position
   dang=rgrb.GetOpeningAngle(rmu,"deg");

   if ((!fDatype && dang>fabs(fDawin)) || (fDatype && dang>fabs(fDawin*totsigma))) continue;

   if (tottfine) tottfine->Fill(dt);
   if (tott) tott->Fill(dt);
   if (tota) tota->Fill(dang);
   if (totcosa) totcosa->Fill(cos(dang*pi/180.));
   if (tot2) tot2->Fill(dang,dt);
 }
}
///////////////////////////////////////////////////////////////////////////
TH1* NcGRB::GetBayesianSignalRate(Double_t p,Double_t& rlow,Double_t& rup,Int_t n)
{
// Provide the Bayesian signal rate and the lower and upper bounds of the
// Bayesian "p%" credible interval [rlow,rup] around the mode of the signal PDF.
//
// Input arguments :
// -----------------
// p    : The percentage of the PDF to be covered by the credible interval around the mode.
//        So for a Gaussian PDF, p=68.3 will result in the [mean-sigma,mean+sigma] 68.3% credible interval.
// rlow : The variable for the return of the lower bound of the credible interval. 
// rup  : The variable for the return of the upper bound of the credible interval. 
// n    : The precision on the result, expressed as 1/n.
//
// By default n=1000 which implies that the accuracy of the result is better than 0.1%.
// Note that very large values of "n" may result in a rather long computation time.
//
// The return argument is the histogram for the Bayesian signal rate PDF.
//
// In case of inconsistent data all returned values are 0.

 rlow=0;
 rup=0;

 TH1* tott=(TH1*)fHistos.FindObject("tott");
 TH1* bkgt=(TH1*)fHistos.FindObject("bkgt");

 if (!tott || !bkgt) return 0;

 Double_t underflow,overflow;
 Int_t nbins=0;
 Double_t nentot=tott->GetEntries();
 nbins=tott->GetNbinsX();
 underflow=tott->GetBinContent(0);
 overflow=tott->GetBinContent(nbins+1);
 nentot=nentot-(underflow+overflow);
 Double_t nenbkg=bkgt->GetEntries();
 nbins=bkgt->GetNbinsX();
 underflow=bkgt->GetBinContent(0);
 overflow=bkgt->GetBinContent(nbins+1);
 nenbkg=nenbkg-(underflow+overflow);

 if (nentot<=0 || nenbkg<=0) return 0;

 // The Bayesian posterior background and signal rate PDFs
 Double_t Non=nentot;
 Double_t Ton=fDtwin*float(fNgrbs);
 Double_t Noff=nenbkg;
 Double_t Toff=Ton;
 TF1 fbkgrpdf=GetBackgroundRatePDF(Noff,Toff);
 TF1 fsigrpdf=GetSignalRatePDF(Non,Ton,Noff,Toff);

 // Determine the "p%" credible interval for the signal rate
 Float_t frac=0;
 frac=GetCredibleInterval(fsigrpdf,p,rlow,rup,n);

 // Provide the signal and background rate PDFs as histograms in the output file
 fbkgrpdf.SetRange(0,3.*Noff/Toff);
 fbkgrpdf.SetNpx(n);
 TH1* hpdfbkgr=(TH1*)fbkgrpdf.GetHistogram()->Clone();
 hpdfbkgr->SetName("hpdfbkgr");
 fHistos.Add(hpdfbkgr);
 fsigrpdf.SetRange(0,3.*Non/Ton);
 fsigrpdf.SetNpx(n);
 TH1* hpdfsigr=(TH1*)fsigrpdf.GetHistogram()->Clone();
 hpdfsigr->SetName("hpdfsigr");
 fHistos.Add(hpdfsigr);

 cout << endl;
 cout << " *" << ClassName() << "::GetBayesianSignalRate* Credible interval [rlow,rup] for p=" << p << "%"
      << " with a precision of 1/" << n << endl;
 cout << " The " << frac << "% credible interval from the Bayesian signal pdf :"
      << " rlow=" << rlow << " rup=" << rup << endl;
 cout << " The following signal and background rate PDF histograms have been generated :" << endl;
 cout << " ... " << hpdfsigr->GetName() << " : " << hpdfsigr->GetTitle() << endl;      
 cout << " ... " << hpdfbkgr->GetName() << " : " << hpdfbkgr->GetTitle() << endl;      

 return hpdfsigr;
}
///////////////////////////////////////////////////////////////////////////
Double_t NcGRB::GetLiMaSignificance() const
{
// Provide the Li-Ma signal significance in terms of the amount of
// standard deviations w.r.t. the "on source" and "off source" observations.
//
// In case of inconsistent data the value 0 is returned.

 Double_t sigma=0;

 TH1* tott=(TH1*)fHistos.FindObject("tott");
 TH1* bkgt=(TH1*)fHistos.FindObject("bkgt");

 if (!tott || !bkgt) return 0;

 Double_t underflow,overflow;
 Int_t nbins=0;
 Double_t nentot=tott->GetEntries();
 nbins=tott->GetNbinsX();
 underflow=tott->GetBinContent(0);
 overflow=tott->GetBinContent(nbins+1);
 nentot=nentot-(underflow+overflow);
 Double_t nenbkg=bkgt->GetEntries();
 nbins=bkgt->GetNbinsX();
 underflow=bkgt->GetBinContent(0);
 overflow=bkgt->GetBinContent(nbins+1);
 nenbkg=nenbkg-(underflow+overflow);

 if (nentot<=0 || nenbkg<=0) return 0;

 // The "on source" and "off source" data
 Int_t Non=int(nentot);
 Double_t Ton=fDtwin*float(fNgrbs);
 Int_t Noff=int(nenbkg);
 Double_t Toff=Ton;

 NcMath m;
 sigma=m.LiMaSignificance(Non,Ton,Noff,Toff);

 cout << endl;
 cout << " *" << ClassName() << "::GetLiMaSignificance* The Li-Ma signal significance is : " << sigma << " sigma." << endl;

 return sigma;
}
///////////////////////////////////////////////////////////////////////////
void NcGRB::GetBayesianPsiStatistics(TString type,Int_t ndt,Double_t nr,Int_t ncut,Int_t freq)
{
// Provide the Bayesian Psi statistics for the (stacked) distributions of the
// observed arrival times and opening angles w.r.t. the corresponding bursts.
//
// Consider a hypothesis B_m representing a counting experiment with m different
// possible outcomes and which is completely defined by the probabilities
// of the various outcomes (and the requirement that the sum of all these
// probabilities equals 1).
// In mathematical terms such a hypothesis belongs to the Bernoulli class B_m. 
//
// The Psi value of n trials of B_m provides (in dB scale) the amount of support
// that the data can maximally give to any Bernoulli class hypothesis different
// from the currently specified B_m.
//
// To be specific : Psi=-10*log[p(D|B_m I)]
//
// where p(D|B_m I) represents the likelihood of the data D under the condition
// that B_m (given some prior information I) are true.
//
// In our current situation, the hypotheses B_m (i.e. Probability Distribution Functions)
// for the various observed distributions are known and "m" just represents the number
// of bins, and "n" represents the number of entries of the corresponding histogram.
// As such, the Psi value (psi0) of the actual observation can be determined.
//
// Further mathematical details can be found in the publication
// N. van Eijndhoven, Astropart. Phys. 28 (2008) 540 (astro-ph/0702029).
//
// This memberfunction may also provide the statistical P-value (i.e. the fraction of
// recorded psi values with psi>=psi0) for the actually observed psi value (psi0)
// based on "nr" repetitions of the counting experiment corresponding to B_m
// with "n" independent random trials.
//
// Input arguments :
// -----------------
// type : "time"  --> Provide statistics for the observed arrival times
//                    This will investigate the deviation from a uniform background time spectrum 
//        "angle" --> Provide statistics for the observed opening angles
//                    This will investigate the deviation from a isotropic background angular spectrum 
//        "cosa"  --> Provide statistics for the cosine of the observed opening angles
//                    This will investigate the deviation from a uniform background cos(angle) spectrum 
//        "dt"    --> Provide statistics for the time intervals between the observed arrival times
//                    This will investigate the deviation from dt spectrum expected from Poisson statistics
// ndt  : Number of events within a dt cell for which the dt statistics will be performed 
// nr   : (Maximum) number of randomised configurations for psi P-value determination.
//        nr<0 implies that no psi P-values will be determined (saves CPU time).
//        nr=0 implies the allowed maximum of 1e19 randomisations.
// ncut : Number of obtained randomised psi entries above the actual observed psi value
//        at which randomisations will be terminated (to save CPU time).
//        ncut=0 implies no early termination.
// freq : Use frequentist's approximation (1) or exact Bayesian expression (0)
//
// The default values are ndt=2, nr=-1, ncut=10 and freq=0.

 NcMath math;

 TString text="none";
 if (type=="time") text="arrival time";
 if (type=="angle") text="opening angle";
 if (type=="cosa") text="cos(opening angle)";
 if (type=="dt") text="arrival time interval";

 cout << endl; 
 if (text=="none")
 {
  cout << " *" << ClassName() << "::GetBayesianPsiStatistics* Unknown statistics type : " << type << endl;
  return;
 }
 else
 {
  cout << " *" << ClassName() << "::GetBayesianPsiStatistics* Analysis of " << text << " statistics" << endl;
 }

 Double_t psitot=-1, psibkg=-1;
 Float_t psidif=0;
 Float_t psimintot=-1, psimaxtot=-1, psifractot=0;
 Float_t psiminbkg=-1, psimaxbkg=-1, psifracbkg=0;
 Double_t nrxtot=-1, nrxbkg=-1;
 Double_t pvaluetot=-1, pvaluebkg=-1;
 TH1F* rtot=0;
 TH1F* rbkg=0;

 ////////////////////////////////////////////
 // Arrival time histo Bayesian statistics //
 ////////////////////////////////////////////
 if (type=="time")
 {
  TH1* tot=(TH1*)fHistos.FindObject("tott");
  TH1* bkg=(TH1*)fHistos.FindObject("bkgt");

  if (!tot || !bkg) return;

  psitot=math.PsiValue(tot,0,0,freq);
  psibkg=math.PsiValue(bkg,0,0,freq);
  psidif=psitot-psibkg;

  // Extreme Psi values for a pure background hypothesis of the recorded arrival time entries
  psimintot=math.PsiExtreme(tot,0,0,-2);
  psimaxtot=math.PsiExtreme(tot,0,0,-1);
  psifractot=(psimaxtot-psitot)/(psimaxtot-psimintot);
  psiminbkg=math.PsiExtreme(bkg,0,0,-2); 
  psimaxbkg=math.PsiExtreme(bkg,0,0,-1);
  psifracbkg=(psimaxbkg-psibkg)/(psimaxbkg-psiminbkg);

  // P-value determination
  if (nr>=0)
  {
   TH1F* hrpsitott=(TH1F*)fHistos.FindObject("hrpsitott");
   TH1F* hrpsibkgt=(TH1F*)fHistos.FindObject("hrpsibkgt");

   if (hrpsitott)
   {
    rtot=(TH1F*)hrpsitott->Clone();
    rtot->Reset();
   }
   else
   {
    rtot=new TH1F("hrpsitott","Random #psi distr. for bkg hypothesis of tott",100,psimintot-1.,psimaxtot+1.);
   }

   if (hrpsibkgt)
   {
    rbkg=(TH1F*)hrpsibkgt->Clone();
    rbkg->Reset();
   }
   else
   {
    rbkg=new TH1F("hrpsibkgt","Random #psi distr. for bkg hypothesis of bkgt",100,psiminbkg-1.,psimaxbkg+1.);
   }

   pvaluetot=math.PsiPvalue(-1,nr,tot,0,0,freq,0,rtot,ncut,&nrxtot);
   pvaluebkg=math.PsiPvalue(-1,nr,bkg,0,0,freq,0,rbkg,ncut,&nrxbkg);
   fHistos.Add(rtot);
   fHistos.Add(rbkg);

   cout << " The following randomised Psi histograms have been generated :" << endl;
   cout << " ... " << rtot->GetName() << " : " << rtot->GetTitle() << endl;      
   cout << " ... " << rbkg->GetName() << " : " << rbkg->GetTitle() << endl;
  }
 }
 
 /////////////////////////////////////////////
 // Opening angle histo Bayesian statistics //
 /////////////////////////////////////////////
 if (type=="angle")
 {
  TH1* tot=(TH1*)fHistos.FindObject("tota");
  TH1* bkg=(TH1*)fHistos.FindObject("bkga");

  if (!tot || !bkg) return;

  TF1 pdfa("pdfa","sin(x*acos(-1.)/180.)");
  psitot=math.PsiValue(tot,0,&pdfa,freq);
  psibkg=math.PsiValue(bkg,0,&pdfa,freq);
  psidif=psitot-psibkg;

  // Extreme Psi values for a pure background hypothesis of the recorded opening angle entries
  psimintot=math.PsiExtreme(tot,0,&pdfa,-2);
  psimaxtot=math.PsiExtreme(tot,0,&pdfa,-1);
  psifractot=(psimaxtot-psitot)/(psimaxtot-psimintot);
  psiminbkg=math.PsiExtreme(bkg,0,&pdfa,-2);
  psimaxbkg=math.PsiExtreme(bkg,0,&pdfa,-1);
  psifracbkg=(psimaxbkg-psibkg)/(psimaxbkg-psiminbkg);

  // P-value determination
  if (nr>=0)
  {
   TH1F* hrpsitota=(TH1F*)fHistos.FindObject("hrpsitota");
   TH1F* hrpsibkga=(TH1F*)fHistos.FindObject("hrpsibkga");

   if (hrpsitota)
   {
    rtot=(TH1F*)hrpsitota->Clone();
    rtot->Reset();
   }
   else
   {
    rtot=new TH1F("hrpsitota","Random #psi distr. for bkg hypothesis of tota",100,psimintot-1.,psimaxtot+1.);
   }

   if (hrpsibkga)
   {
    rbkg=(TH1F*)hrpsibkga->Clone();
    rbkg->Reset();
   }
   else
   {
    rbkg=new TH1F("hrpsibkga","Random #psi distr. for bkg hypothesis of bkga",100,psiminbkg-1.,psimaxbkg+1.);
   }

   pvaluetot=math.PsiPvalue(-1,nr,tot,0,&pdfa,freq,0,rtot,ncut,&nrxtot);
   pvaluebkg=math.PsiPvalue(-1,nr,bkg,0,&pdfa,freq,0,rbkg,ncut,&nrxbkg);
   fHistos.Add(rtot);
   fHistos.Add(rbkg);

   cout << " The following randomised Psi histograms have been generated :" << endl;
   cout << " ... " << rtot->GetName() << " : " << rtot->GetTitle() << endl;      
   cout << " ... " << rbkg->GetName() << " : " << rbkg->GetTitle() << endl;
  }      
 }

 ///////////////////////////////////////////////////////
 // Cosine of opening angle histo Bayesian statistics //
 ///////////////////////////////////////////////////////
 if (type=="cosa")
 {
  TH1* tot=(TH1*)fHistos.FindObject("totcosa");
  TH1* bkg=(TH1*)fHistos.FindObject("bkgcosa");

  if (!tot || !bkg) return;

  psitot=math.PsiValue(tot,0,0,freq);
  psibkg=math.PsiValue(bkg,0,0,freq);
  psidif=psitot-psibkg;

  // Extreme Psi values for a pure background hypothesis of the recorded cos(opening angle) entries
  psimintot=math.PsiExtreme(tot,0,0,-2);
  psimaxtot=math.PsiExtreme(tot,0,0,-1);
  psifractot=(psimaxtot-psitot)/(psimaxtot-psimintot);
  psiminbkg=math.PsiExtreme(bkg,0,0,-2);
  psimaxbkg=math.PsiExtreme(bkg,0,0,-1);
  psifracbkg=(psimaxbkg-psibkg)/(psimaxbkg-psiminbkg);

  // P-value determination
  if (nr>=0)
  {
   TH1F* hrpsitotcosa=(TH1F*)fHistos.FindObject("hrpsitotcosa");
   TH1F* hrpsibkgcosa=(TH1F*)fHistos.FindObject("hrpsibkgcosa");

   if (hrpsitotcosa)
   {
    rtot=(TH1F*)hrpsitotcosa->Clone();
    rtot->Reset();
   }
   else
   {
    rtot=new TH1F("hrpsitotcosa","Random #psi distr. for bkg hypothesis of totcosa",100,psimintot-1.,psimaxtot+1.);
   }

   if (hrpsibkgcosa)
   {
    rbkg=(TH1F*)hrpsibkgcosa->Clone();
    rbkg->Reset();
   }
   else
   {
    rbkg=new TH1F("hrpsibkgcosa","Random #psi distr. for bkg hypothesis of bkgcosa",100,psiminbkg-1.,psimaxbkg+1.);
   }
    pvaluetot=math.PsiPvalue(-1,nr,tot,0,0,freq,0,rtot,ncut,&nrxtot);
    pvaluebkg=math.PsiPvalue(-1,nr,bkg,0,0,freq,0,rbkg,ncut,&nrxbkg);
    fHistos.Add(rtot);
    fHistos.Add(rbkg);

    cout << " The following randomised Psi histograms have been generated :" << endl;
    cout << " ... " << rtot->GetName() << " : " << rtot->GetTitle() << endl;      
    cout << " ... " << rbkg->GetName() << " : " << rbkg->GetTitle() << endl;
  }
 }

 ///////////////////////////////////////////////
 // Arrival time interval Bayesian statistics //
 ///////////////////////////////////////////////
 if (type=="dt")
 {
  TH1* tot=(TH1*)fHistos.FindObject("tottfine");
  TH1* bkg=(TH1*)fHistos.FindObject("bkgtfine");

  if (!tot || !bkg) return;

  // Create the delta t histograms
  TString nametot="htotdt";
  nametot+=ndt;
  TString namebkg="hbkgdt";
  namebkg+=ndt;
  TH1F* htotdt=(TH1F*)fHistos.FindObject(nametot.Data());
  TH1F* hbkgdt=(TH1F*)fHistos.FindObject(namebkg.Data());
  Double_t deltatbin=0, deltatmin=0, deltatmax=0;
  if (!htotdt && !hbkgdt)
  {
   htotdt=(TH1F*)(GetDxHistogram(tot,ndt,-1,0,-1).Clone(nametot.Data()));
   deltatbin=htotdt->GetXaxis()->GetBinWidth(1);
   deltatmin=htotdt->GetXaxis()->GetXmin();
   deltatmax=htotdt->GetXaxis()->GetXmax();
   hbkgdt=(TH1F*)(GetDxHistogram(bkg,ndt,deltatbin,deltatmin,deltatmax).Clone(namebkg.Data()));

   // Create titles and labels for the delta t histograms
   TString title="dt to contain ";
   title+=ndt;
   title+=" mu-up in total twindow";
   title+=";Upgoing #mu dt (in sec);Counts per bin of size %-10.3g";
   TString s=title.Format(title.Data(),deltatbin);
   htotdt->SetTitle(s.Data());

   title="dt to contain ";
   title+=ndt;
   title+=" mu-up in bkg twindow";
   title+=";Upgoing #mu dt (in sec);Counts per bin of size %-10.3g";
   s=title.Format(title.Data(),deltatbin);
   hbkgdt->SetTitle(s.Data());

   fHistos.Add(htotdt);
   fHistos.Add(hbkgdt);

   cout << " The following arrival time interval (dt) histograms have been generated :" << endl;
   cout << " ... " << htotdt->GetName() << " : " << htotdt->GetTitle() << endl;      
   cout << " ... " << hbkgdt->GetName() << " : " << hbkgdt->GetTitle() << endl;
  }

  // Creation of the Poisson based dt PDFs from the observed data for a background only hypothesis 
  Double_t underflow,overflow;
  Int_t nbins=0;
  Double_t nentot=tot->GetEntries();
  nbins=tot->GetNbinsX();
  underflow=tot->GetBinContent(0);
  overflow=tot->GetBinContent(nbins+1);
  nentot=nentot-(underflow+overflow);
  Double_t nenbkg=bkg->GetEntries();
  nbins=bkg->GetNbinsX();
  underflow=bkg->GetBinContent(0);
  overflow=bkg->GetBinContent(nbins+1);
  nenbkg=nenbkg-(underflow+overflow);

  if (nentot<=0 || nenbkg<=0) return;

  Double_t ratetot=nentot/(fDtwin);
  Double_t ratebkg=nenbkg/(fDtwin);

  // Determine the corresponding dt PDFs based on Poisson statistics
  // Only the bkg dt PDF is used, since this may be obtained from off-source measurements.
  // Using the total dt PDF would artificially lower the sensitivity due to possible signal events.

  TF1 fdttot=math.PoissonDtDist(ratetot,ndt); // Only for reference, not used in the analysis
  TF1 fdtbkg=math.PoissonDtDist(ratebkg,ndt);

  nametot="hpdftotdt";
  nametot+=ndt;
  namebkg="hpdfbkgdt";
  namebkg+=ndt;
  TH1* hpdftotdt=(TH1*)fHistos.FindObject(nametot.Data());
  TH1* hpdfbkgdt=(TH1*)fHistos.FindObject(namebkg.Data());

  // Provide the dt PDFs as histograms in the output file
  if (!hpdftotdt && !hpdfbkgdt)
  {
   deltatmax=htotdt->GetXaxis()->GetXmax();
   fdttot.SetRange(0,deltatmax);
   fdttot.SetNpx(10000);
   hpdftotdt=(TH1*)fdttot.GetHistogram()->Clone();
   hpdftotdt->SetName(nametot.Data());
   deltatmax=hbkgdt->GetXaxis()->GetXmax();
   fdtbkg.SetRange(0,deltatmax);
   fdtbkg.SetNpx(10000);
   hpdfbkgdt=(TH1*)fdtbkg.GetHistogram()->Clone();
   hpdfbkgdt->SetName(namebkg.Data());
   fHistos.Add(hpdftotdt);
   fHistos.Add(hpdfbkgdt);

   cout << " The following arrival time interval (dt) PDFs have been generated :" << endl;
   cout << " ... " << hpdftotdt->GetName() << " : " << hpdftotdt->GetTitle() << endl;      
   cout << " ... " << hpdfbkgdt->GetName() << " : " << hpdfbkgdt->GetTitle() << endl;
  }

  psitot=math.PsiValue(htotdt,0,&fdtbkg,freq);
  psibkg=math.PsiValue(hbkgdt,0,&fdtbkg,freq);
  psidif=psitot-psibkg;

  psimintot=math.PsiExtreme(htotdt,0,&fdtbkg,-2);
  psimaxtot=math.PsiExtreme(htotdt,0,&fdtbkg,-1);
  psifractot=(psimaxtot-psitot)/(psimaxtot-psimintot);
  psiminbkg=math.PsiExtreme(hbkgdt,0,&fdtbkg,-2);
  psimaxbkg=math.PsiExtreme(hbkgdt,0,&fdtbkg,-1);
  psifracbkg=(psimaxbkg-psibkg)/(psimaxbkg-psiminbkg);

  // P-value determination
  if (nr>=0)
  {
   nametot="hrpsitotdt";
   nametot+=ndt;
   namebkg="hrpsibkgdt";
   namebkg+=ndt;

   TH1F* hrpsitotdt=(TH1F*)fHistos.FindObject(nametot.Data());
   TH1F* hrpsibkgdt=(TH1F*)fHistos.FindObject(namebkg.Data());

   TString title;
   if (hrpsitotdt)
   {
    rtot->Reset();
   }
   else
   {
    title="Random #psi distr. for bkg hypothesis of tot dt for n=";
    title+=ndt;
    rtot=new TH1F(nametot.Data(),title.Data(),100,psimintot-1.,psimaxtot+1.);
   }

   if (hrpsibkgdt)
   {
    rbkg->Reset();
   }
   else
   {
    title="Random #psi distr. for bkg hypothesis of bkg dt for n=";
    title+=ndt;
    rbkg=new TH1F(namebkg.Data(),title.Data(),100,psiminbkg-1.,psimaxbkg+1.);
   }

   pvaluetot=math.PsiPvalue(-1,nr,htotdt,0,&fdtbkg,freq,0,rtot,ncut,&nrxtot);
   pvaluebkg=math.PsiPvalue(-1,nr,hbkgdt,0,&fdtbkg,freq,0,rbkg,ncut,&nrxbkg);
   fHistos.Add(rtot);
   fHistos.Add(rbkg);

   cout << " The following randomised Psi histograms have been (re)generated :" << endl;
   cout << " ... " << rtot->GetName() << " : " << rtot->GetTitle() << endl;      
   cout << " ... " << rbkg->GetName() << " : " << rbkg->GetTitle() << endl;
  }
 }

 // Listing of the statistics results
 cout << " *** Observed Psi values (in dB) for the hypothesis of no GRB signal ***" << endl;
 cout << " For the \"on source\" stacked patches : psi = " << psitot << endl;
 cout << " For the corresponding \"opposite RA\" stacked \"off source\" patches : psi = " << psibkg << endl;
 cout << " --> Difference between observed \"on source\" and \"off source\" psi values : " << psidif << endl;
 cout << " *** Extreme Psi values for the case of pure background ***" << endl;
 cout << " For \"on source\"  psimin : " << psimintot << " psimax : " << psimaxtot << " (psimax-psi)/range : " << psifractot << endl;
 cout << " For \"off source\" psimin : " << psiminbkg << " psimax : " << psimaxbkg << " (psimax-psi)/range : " << psifracbkg << endl;

 if (nr>=0)
 {
  cout << " *** P-values of the observed \"on source\" and \"off source\" psi values ***" << endl;
  cout << " For the \"on source\"  stacked patches : P-value = " << pvaluetot << " Used number of randomisations : " << nrxtot << endl;
  cout << " For the \"off source\" stacked patches : P-value = " << pvaluebkg << " Used number of randomisations : " << nrxbkg << endl;
 }
}
///////////////////////////////////////////////////////////////////////////
void NcGRB::GetChi2Statistics(TString type,Int_t ndt)
{
// Provide the Chi-squared statistics for the (stacked) distributions of the
// observed arrival times and opening angles w.r.t. the corresponding bursts.
//
// Input arguments :
// -----------------
// type : "time"  --> Provide statistics for the observed arrival times
//                    This will investigate the deviation from a uniform background time spectrum 
//        "angle" --> Provide statistics for the observed opening angles
//                    This will investigate the deviation from a isotropic background angular spectrum 
//        "cosa"  --> Provide statistics for the cosine of the observed opening angles
//                    This will investigate the deviation from a uniform background cos(angle) spectrum 
//        "dt"    --> Provide statistics for the time intervals between the observed arrival times
//                    This will investigate the deviation from dt spectrum expected from Poisson statistics
// ndt  : Number of events within a dt cell for which the dt statistics will be performed 
//
// The default value is ndt=2.

 NcMath math;

 TString text="none";
 if (type=="time") text="arrival time";
 if (type=="angle") text="opening angle";
 if (type=="cosa") text="cos(opening angle)";
 if (type=="dt") text="arrival time interval";

 cout << endl; 
 if (text=="none")
 {
  cout << " *" << ClassName() << "::GetChi2Statistics* Unknown statistics type : " << type << endl;
  return;
 }
 else
 {
  cout << " *" << ClassName() << "::GetChi2Statistics* Analysis of " << text << " statistics" << endl;
 }

 TH1F* rtot=0;
 TH1F* rbkg=0;
 Int_t ndftot=0;
 Int_t ndfbkg=0;
 Float_t chitot=0;
 Float_t chibkg=0;

 ///////////////////////////////////////////////
 // Arrival time histo Chi-squared statistics //
 ///////////////////////////////////////////////
 if (type=="time")
 {
  TH1* tot=(TH1*)fHistos.FindObject("tott");
  TH1* bkg=(TH1*)fHistos.FindObject("bkgt");

  if (!tot || !bkg) return;

  chitot=math.Chi2Value(tot,0,0,&ndftot);
  chibkg=math.Chi2Value(bkg,0,0,&ndfbkg);
 }

 ////////////////////////////////////////////////
 // Opening angle histo Chi-squared statistics //
 ////////////////////////////////////////////////
 if (type=="angle")
 {
  TH1* tot=(TH1*)fHistos.FindObject("tota");
  TH1* bkg=(TH1*)fHistos.FindObject("bkga");

  if (!tot || !bkg) return;

  TF1 pdf("pdf","sin(x*acos(-1.)/180.)");
  chitot=math.Chi2Value(tot,0,&pdf,&ndftot);
  chibkg=math.Chi2Value(bkg,0,&pdf,&ndfbkg);
 }

 //////////////////////////////////////////////////////////
 // Cosine of opening angle histo Chi-squared statistics //
 //////////////////////////////////////////////////////////
 if (type=="cosa")
 {
  TH1* tot=(TH1*)fHistos.FindObject("totcosa");
  TH1* bkg=(TH1*)fHistos.FindObject("bkgcosa");

  if (!tot || !bkg) return;

  chitot=math.Chi2Value(tot,0,0,&ndftot);
  chibkg=math.Chi2Value(bkg,0,0,&ndfbkg);
 }

 ///////////////////////////////////////////////
 // Arrival time interval Bayesian statistics //
 ///////////////////////////////////////////////
 if (type=="dt")
 {
  TH1* tot=(TH1*)fHistos.FindObject("tottfine");
  TH1* bkg=(TH1*)fHistos.FindObject("bkgtfine");

  if (!tot || !bkg) return;
 
  // Create the delta t histograms
  TString nametot="htotdt";
  nametot+=ndt;
  TString namebkg="hbkgdt";
  namebkg+=ndt;
  TH1F* htotdt=(TH1F*)fHistos.FindObject(nametot.Data());
  TH1F* hbkgdt=(TH1F*)fHistos.FindObject(namebkg.Data());
  Double_t deltatbin=0, deltatmin=0, deltatmax=0;
  if (!htotdt && !hbkgdt)
  {
   htotdt=(TH1F*)(GetDxHistogram(tot,ndt,-1,0,-1).Clone(nametot.Data()));
   deltatbin=htotdt->GetXaxis()->GetBinWidth(1);
   deltatmin=htotdt->GetXaxis()->GetXmin();
   deltatmax=htotdt->GetXaxis()->GetXmax();
   hbkgdt=(TH1F*)(GetDxHistogram(bkg,ndt,deltatbin,deltatmin,deltatmax).Clone(namebkg.Data()));

   // Create titles and labels for the delta t histograms
   TString title="dt to contain ";
   title+=ndt;
   title+=" mu-up in total twindow";
   title+=";Upgoing #mu dt (in sec);Counts per bin of size %-10.3g";
   TString s=title.Format(title.Data(),deltatbin);
   htotdt->SetTitle(s.Data());

   title="dt to contain ";
   title+=ndt;
   title+=" mu-up in bkg twindow";
   title+=";Upgoing #mu dt (in sec);Counts per bin of size %-10.3g";
   s=title.Format(title.Data(),deltatbin);
   hbkgdt->SetTitle(s.Data());

   fHistos.Add(htotdt);
   fHistos.Add(hbkgdt);

   cout << " The following arrival time interval (dt) histograms have been generated :" << endl;
   cout << " ... " << htotdt->GetName() << " : " << htotdt->GetTitle() << endl;      
   cout << " ... " << hbkgdt->GetName() << " : " << hbkgdt->GetTitle() << endl;
  }

  // Creation of the Poisson based dt PDFs from the observed data for a background only hypothesis 
  Double_t underflow,overflow;
  Int_t nbins=0;
  Double_t nentot=tot->GetEntries();
  nbins=tot->GetNbinsX();
  underflow=tot->GetBinContent(0);
  overflow=tot->GetBinContent(nbins+1);
  nentot=nentot-(underflow+overflow);
  Double_t nenbkg=bkg->GetEntries();
  nbins=bkg->GetNbinsX();
  underflow=bkg->GetBinContent(0);
  overflow=bkg->GetBinContent(nbins+1);
  nenbkg=nenbkg-(underflow+overflow);

  if (nentot<=0 || nenbkg<=0) return;

  Double_t ratetot=nentot/(fDtwin);
  Double_t ratebkg=nenbkg/(fDtwin);

  TF1 fdttot=math.PoissonDtDist(ratetot,ndt); // Only for reference, not used in the analysis
  TF1 fdtbkg=math.PoissonDtDist(ratebkg,ndt);

  nametot="hpdftotdt";
  nametot+=ndt;
  namebkg="hpdfbkgdt";
  namebkg+=ndt;
  TH1* hpdftotdt=(TH1*)fHistos.FindObject(nametot.Data());
  TH1* hpdfbkgdt=(TH1*)fHistos.FindObject(namebkg.Data());

  // Provide the dt PDFs as histograms in the output file
  if (!hpdftotdt && !hpdfbkgdt)
  {
   deltatmax=htotdt->GetXaxis()->GetXmax();
   fdttot.SetRange(0,deltatmax);
   fdttot.SetNpx(10000);
   hpdftotdt=(TH1*)fdttot.GetHistogram()->Clone();
   hpdftotdt->SetName(nametot.Data());
   deltatmax=hbkgdt->GetXaxis()->GetXmax();
   fdtbkg.SetRange(0,deltatmax);
   fdtbkg.SetNpx(10000);
   hpdfbkgdt=(TH1*)fdtbkg.GetHistogram()->Clone();
   hpdfbkgdt->SetName(namebkg.Data());
   fHistos.Add(hpdftotdt);
   fHistos.Add(hpdfbkgdt);

   cout << " The following arrival time interval (dt) PDFs have been generated :" << endl;
   cout << " ... " << hpdftotdt->GetName() << " : " << hpdftotdt->GetTitle() << endl;      
   cout << " ... " << hpdfbkgdt->GetName() << " : " << hpdfbkgdt->GetTitle() << endl;
  }

  chitot=math.Chi2Value(htotdt,0,&fdttot,&ndftot);
  chibkg=math.Chi2Value(hbkgdt,0,&fdtbkg,&ndfbkg);
 }

 // Listing of the statistics results
 Float_t chidif=chitot-chibkg;
 cout << " *** Observed Chi-squared values for the hypothesis of no GRB signal ***" << endl;
 cout << " For the \"on source\" stacked patches : chi2 = " << chitot << " ndf = " << ndftot << endl;
 cout << " For the corresponding \"opposite RA\" stacked \"off source\" patches : chi2 = " << chibkg << " ndf = " << ndfbkg << endl;
 cout << " --> Difference between observed \"on source\" and \"off source\" chi2 values : " << chidif << endl;

 Float_t ptot=math.Chi2Pvalue(chitot,ndftot);
 Float_t sigmatot=math.Chi2Pvalue(chitot,ndftot,0,1);
 Float_t pbkg=math.Chi2Pvalue(chibkg,ndfbkg);
 Float_t sigmabkg=math.Chi2Pvalue(chibkg,ndfbkg,0,1);

 cout << " *** P-values of the observed \"on source\" and \"off source\" chi2 values ***" << endl;
 cout << " For the \"on source\"  stacked patches : P-value = " << ptot << " (" << sigmatot << " sigma)" << endl;
 cout << " For the \"off source\" stacked patches : P-value = " << pbkg << " (" << sigmabkg << " sigma)" << endl;
}
///////////////////////////////////////////////////////////////////////////
TObject* NcGRB::Clone(const char* name) const
{
// Make a deep copy of the current object and provide the pointer to the copy.
// This memberfunction enables automatic creation of new objects of the
// correct type depending on the object type, a feature which may be very useful
// for containers when adding objects in case the container owns the objects.

 NcGRB* grb=new NcGRB(*this);
 if (name)
 {
  if (strlen(name)) grb->SetName(name);
 }
 return grb;
}
///////////////////////////////////////////////////////////////////////////
