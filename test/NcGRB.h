#ifndef NcGRB_h
#define NcGRB_h
// Copyright(c) 1997-2009, NCFS, All Rights Reserved.
// See cxx source for full Copyright notice.

#include <iostream>
#include <fstream>

#include "TStyle.h"
#include "TFile.h"
#include "TChain.h"

#include "NcJob.h"
#include "NcAstrolab.h"
#include "NcSample.h"

class NcGRB : public NcAstrolab
{
 public:
  NcGRB(const char* name="NcGRB",const char* title="Cosmic coincidence analysis"); // Default constructor
  virtual ~NcGRB();                                             // Default destructor
  NcGRB(const NcGRB& q);                                        // Copy constructor
  Int_t SetParameter(TString name,Double_t value);              // Specification of a certain parameter setting
  void LoadBurstData(TString file,TString tree,Int_t date1=0,Int_t date2=0,Int_t nmax=-1);  // Load observed burst data
  void GenBurstData(Int_t n);                                   // Generate fictative burst data
  void PrintSettings() const;                                   // Print all parameter settings
  void ListHistograms() const;                                  // Provide a list of all the stored histograms
  void WriteHistograms(TString filename);                       // Write all stored histograms to a ROOT output file
  void MakeZdist(TString file,TString tree,TString branch,Int_t nb=200,Float_t zmin=0,Float_t zmax=20); // Make observed redshift distribution
  void MakeT90dist(TString file,TString tree,TString branch,Int_t nb=50,Float_t xmin=-5,Float_t xmax=5); // Make observed T90 distribution
  TH1* GetBayesianSignalRate(Double_t p,Double_t& rlow,Double_t& rup,Int_t n=1000); // Provide Bayesian signal rate and credible interval 
  Double_t GetLiMaSignificance() const;                         // Provide the Li-Ma signal significance
  void GetBayesianPsiStatistics(TString type,Int_t ndt=2,Double_t nr=-1,Int_t ncut=10,Int_t freq=0); // Provide Bayesian Psi statistics
  void GetChi2Statistics(TString type,Int_t ndt=2);             // Provide Chi-squared statistics
  virtual void Exec(Option_t* opt);                             // Perform the analysis
  virtual TObject* Clone(const char* name="") const;            // Make a deep copy and provide its pointer
 
 protected:
  Int_t fNmax;       // Maximal number of GRBs to be accepted for analysis (<0 : no limitation)
  Float_t fDeclmin;  // Minimal declination (in degrees) for GRB position acceptance
  Float_t fDeclmax;  // Maximal declination (in degrees) for GRB position acceptance
  Float_t fT90min;   // Minimal duration (t90 in sec) for GRB acceptance
  Float_t fT90max;   // Maximal duration (t90 in sec) for GRB acceptance
  Float_t fZmin;     // Minimal redshift for GRB acceptance (<0 : [|fZmin|,fZmax] with random z when redshift is unknown)
  Float_t fZmax;     // Maximal redshift for GRB acceptance
  Float_t fSigmagrb; // Angular uncertainty (sigma in degrees) on GRB position (<0 : use GCN data) 
  Float_t fMaxsigma; // Maximal combined GRB position and track angular uncertainty (sigma in degrees) for acceptance
  Float_t fGrbnu;    // Maximum number of detectable neutrinos per GRB (<0 : no statistical fluctuations)
  Float_t fAvgrbz;   // Average GRB redshift (<0 : determine from observations)
  Float_t fAvgrbt90; // Average GRB duration (T90) in seconds (<0 : determine from observations)
  Int_t fInburst;    // Flag to indicate that neutrinos are produced coupled (1) or not (0) to the gamma flare duration
  Float_t fDtnu;     // Mean time diff. (in sec) between gamma and nu production (decoupled) or in T90 units w.r.t. trigger (coupled)
  Float_t fDtnus;    // Sigma of time difference (in sec) between gamma and nu production (<0 is in T90 units)
  Float_t fAngres;   // Detector angular resolution (degrees)
  Float_t fTimres;   // Detector time resolution (sec)
  Float_t fBkgrate;  // Mean rate (in Hz) of upgoing bkg muons
  Int_t fNevtdt;     // Number of events within a dt cell for which the inter-muon dt statistics will be performed 
  Float_t fDtwin;    // Total search time window (in sec) centered at GRB trigger
  Float_t fDawin;    // Angular search circle (<0 is decl. band) in degrees or sigma around (above/below) GRB position
  Int_t fDatype;     // Type of angular window specification (0=in degrees 1=in units of combined GRB/track sigma) 
  Float_t fNbkg;     // Mean number of counts per bin for auto-binning
  Float_t fTbint90;  // Time bin size in units of average T90 (0 : Time bin size determined by fTbin) 
  Float_t fTbin;     // Time bin size in seconds (0=variable bins  <0 will result in a mean fNbkg counts/bin)
  Float_t fVarTbin;  // Size (in sec) of the first time bin in case of variable time bins
  Float_t fAbin;     // Angular bin size in degrees (<0 will result in a mean fNbkg counts/bin)
  Int_t fFreq;       // Use frequentist's approximation (1) or exact Bayesian expression (0)
  Int_t fNpsi;       // Number of psi entries for bkg psi-value distributions (<0 : time shuffling)
  Int_t fUsetott;    // Use the observed tott number of entries in case of time shuffling 
  Int_t fGrbpos;     // Use the original burst locations (1) or random ones (0) for bkg studies
  Float_t fNrandom;  // (Maximum) number of randomised configurations for direct psi P-value determination (0 means n=1e19)
  Int_t fNcut;       // Number of obtained psi>psi0 entries at which randomisations will be terminated (0= no early termination)

  Int_t fNgrbs;         // Number of GRBs accepted for analysis
  Float_t fMaxtotsigma; // Maximum of the encountered totsigma values
  Float_t fMup;         // Resulting number of @@@@
  Float_t fNmupday;     // Resulting number of @@@@

  TObjArray fHistos; // Storage of all the produced histograms
  void Compensate(Int_t& nmugrb); // Compensate for statistical underfluctuation
 
 ClassDef(NcGRB,1) // Perform coincidence studies of (transient) cosmic phenomena.
};
#endif
