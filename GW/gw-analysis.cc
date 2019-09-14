/////////////////////////////////////////////////////////////
// ROOT macro to investigate the neutrino time clustering  //
// around the Gravitational Wave event GW170608.           //
//                                                         //
// To run the macro, just invoke the command ($=prompt)    //
// $root -b -q gw-analysis.cc                              //
//                                                         //
// All produced histograms are available in the produced   //
// output file "output.root".                              //
//                                                         //
// Note: This macro needs the NCFSPack library, which is   //
// available on the IIHE cluster after issuing the command //
// $ source /ice3/software/iihe/ncfs.sh                    // 
//                                                         //
// Nick van Eijndhoven, 17-sep-2018, IIHE-VUB Brussel.     //
/////////////////////////////////////////////////////////////
{
 gSystem->Load("ncfspack");

 // Provide various mathematical tools
 NcMath math;

 // The virtual lab providing various analysis tools
 NcAstrolab lab;

 // The timestamp of the GW170608 event
 NcTimestamp tsgw;
 tsgw.SetUT("2017-06-08","02:01:16.0",1);

 // The internal storage array for the produced histograms
 TObjArray histos;

 // The data file for observed neutrino arrival times
 ifstream nudata;
 nudata.clear();
 nudata.open("GW170608-IceCube.txt");
 if (!nudata.good())
 {
  cout << " *** Data file for observed neutrino arrival times not found ***" << endl;
  return;
 }
 nudata.seekg(0); // Position at begin of file

 // Obtain  the relative Icecube event times
 TString date;
 TString time;
 NcTimestamp tsnu;
 Float_t dt=0;
 Float_t dtmin=1e10;
 NcSample nutimes;
 nutimes.SetStoreMode();
 while (nudata >> date >> time)
 {
  tsnu.SetUT(date,time,1);
  dt=tsgw.GetDifference(tsnu,"s");
  nutimes.Enter(dt);
  if (fabs(dt)<dtmin) dtmin=fabs(dt);
 }

 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 // Creation of the histo of the observed IceCube event times w.r.t. the GW trigger time
 Float_t tbin=dtmin/2.; // Bin size (in sec) for the relative event time histo
 if (tbin>0.1) tbin=0.1;
 Float_t twinlow=-500;  // Lower bound (in sec) of the time window around the GW trigger time
 Float_t twinup=500;    // Upper bound (in sec) of the time window around the GW trigger time
 Float_t tnuwin=twinup-twinlow;
 Int_t ntbins=int(tnuwin/tbin);
 TString title,s;
 title="IceCube event times w.r.t. GW trigger";
 title+=";Event time (in sec) w.r.t. GW trigger;Counts per bin of size %-10.3g";
 s=title.Format(title.Data(),tbin);
 TH1F* htnu=new TH1F("htnu",s.Data(),ntbins,twinlow,twinup);
 histos.Add(htnu);

 // Fill the histo with the relative IceCube event times
 for (Int_t i=1; i<=nutimes.GetN(); i++)
 {
  dt=nutimes.GetEntry(i,1);
  htnu->Fill(dt);
 }
 
 //////////////////////////////////
 // Create the delta t histogram //
 //////////////////////////////////

 Int_t nevtdt=2; // Number of events within a dt cell for which the inter-event dt statistics will be performed 
 TH1F* hdtnu=(TH1F*)(lab.GetDxHistogram(htnu,nevtdt,-1,0,-1).Clone("hdtnu"));
//@@@ TH1F* hdtnu=(TH1F*)(lab.GetDxHistogram(htnu,nevtdt,0,0,-1).Clone("hdtnu"));
 histos.Add(hdtnu);

 Int_t nbdtnu=hdtnu->GetNbinsX();
 Double_t deltatbin=hdtnu->GetXaxis()->GetBinWidth(1);
 Double_t deltatmin=hdtnu->GetXaxis()->GetXmin();
 Double_t deltatmax=hdtnu->GetXaxis()->GetXmax();

 // Create title and labels for the delta t histogram
 title="Delta t distribution to contain exactly ";
 title+=nevtdt;
 title+=" events";
 title+=";#Deltat (in sec);Counts per bin of size %-10.3g";
 s=title.Format(title.Data(),deltatbin);

 hdtnu->SetTitle(s.Data());

 ////////////////////////////////////////////////////////////
 // Specification of the signal and background event rates //
 ////////////////////////////////////////////////////////////

 Int_t nbtnu=htnu->GetNbinsX();
 Float_t nentnu=htnu->GetEntries();
 Float_t underflow=htnu->GetBinContent(0);
 Float_t overflow=htnu->GetBinContent(nbtnu+1);
 nentnu=nentnu-(underflow+overflow);

 // The Bayesian posterior background and signal rate PDFs
 Int_t Non=nentnu;
 Double_t Ton=tnuwin;
 Int_t Noff=545;   // The number of recorded off-time events
 Double_t Toff=28692; // The off-time recording period
//@@@@ For consitency check only
//@@@ Noff=Non;
//@@@ Toff=Ton;
//@@@@
 TF1 fbkgrpdf=lab.GetBackgroundRatePDF(Noff,Toff);
 TF1 fsigrpdf=lab.GetSignalRatePDF(Non,Ton,Noff,Toff);

 // Provide the rate PDFs as histograms in the output file
 fbkgrpdf.SetRange(0,3.*float(Noff)/Toff);
 fbkgrpdf.SetNpx(10000);
 TH1* hpdfbkgr=(TH1*)fbkgrpdf.GetHistogram()->Clone();
 hpdfbkgr->SetName("hpdfbkgr");
 histos.Add(hpdfbkgr);
 fsigrpdf.SetRange(0,3.*float(Non)/Ton);
 fsigrpdf.SetNpx(10000);
 TH1* hpdfsigr=(TH1*)fsigrpdf.GetHistogram()->Clone();
 hpdfsigr->SetName("hpdfsigr");
 histos.Add(hpdfsigr);

 // Creation of the Poisson based dt PDF for a background only hypothesis
 Float_t ratebkg=float(Noff)/Toff;
 TF1 fdtbkg=math.PoissonDtDist(ratebkg,nevtdt);
 fdtbkg.SetRange(0,deltatmax);

 // Provide the dt PDF histogram with the same binning as the Delta t distribution
 fdtbkg.SetNpx(nbdtnu);
 TH1* hpdfbkgdt=(TH1*)fdtbkg.GetHistogram()->Clone();
 hpdfbkgdt->SetName("hpdfbkgdt");
 histos.Add(hpdfbkgdt);

 // Provide the cumulative distributions for the dt PDF and observed distribution
 TH1F hpdfbkgdtCDH=lab.GetCumulHistogram(hpdfbkgdt,"hpdfbkgdtCDH","FN");
 histos.Add(&hpdfbkgdtCDH);
 TH1F hdtnuCDH=lab.GetCumulHistogram(hdtnu,"hdtnuCDH","FN");
 histos.Add(&hdtnuCDH);

 // Provide the dt PDF histogram with fine binning
 fdtbkg.SetNpx(10000);
 TH1* hpdfbkgdtfine=(TH1*)fdtbkg.GetHistogram()->Clone();
 hpdfbkgdtfine->SetName("hpdfbkgdtfine");
 histos.Add(hpdfbkgdtfine);

 ////////////////////////////////////////////////////////////////////////////////
 // Statistical evaluation of the signal and background observations           //
 //                                                                            //
 // Determination of the Bayesian psi value assumes that there is no signal.   //
 // This corresponds to searching out the Bernoulli class B_m                  //
 // with m=nbins of the histogram.                                             //
 ////////////////////////////////////////////////////////////////////////////////

 Double_t nrandom=1e3; // (Maximum) number of randomised configurations for direct psi P-value determination (0 means n=1e19)
 Int_t ncut=10;        // Number of obtained psi>psi0 entries at which randomisations will be terminated (0= no early termination)

 Float_t onrate=float(Non)/Ton;
 cout << endl;
 cout << " *** Statistics of the observed event samples ***" << endl;
 cout << " *On time* Number of entries : " << Non << " time window (sec) : " << Ton << " --> Event rate (Hz) : " << onrate << endl;
 cout << "  (Used bin sizes (sec) for the distributions time : " << tbin << " Delta t : " << deltatbin << ")" << endl;
 cout << " *Off time* Number of entries : " << Noff << " time window (sec) : " << Toff
      <<  " --> Event rate (Hz) : " << ratebkg << endl;

 // The Li-Ma significance for the observation w.r.t. the background
 Float_t sigma=math.LiMaSignificance(Non,Ton,Noff,Toff);
 cout << endl;
 cout << " Li-Ma significance of the observed event rate is " << sigma << " sigma." << endl;

 // The Bayesian 90% credible interval for the signal rate
 Double_t rlow,rup;
 Float_t frac;
 frac=lab.GetCredibleInterval(fsigrpdf,90,rlow,rup);
 frac=frac*100.;
 cout << endl;
 cout << " === The " << frac << "% credible interval (in Hz) from the Bayesian signal rate pdf :"
      << " rlow=" << rlow << " rup=" << rup << endl;

 // IceCube event time Bayesian psi statistics
/***
 Double_t psitnu=math.PsiValue(htnu,0,0);

 Float_t psimintnu=math.PsiExtreme(htnu,0,0,-2);
 Float_t psimaxtnu=math.PsiExtreme(htnu,0,0,-1);

 TH1F* hrpsitnu=new TH1F("hrpsitnu","Random #psi distr. for bkg hypothesis of event times",100,psimintnu-1.,psimaxtnu+1.);
 histos.Add(hrpsitnu);

 Double_t nrxtnu=-1; // Returned number of actually performed randomizations
 Double_t pvaluetnu=math.PsiPvalue(-1,nrandom,htnu,0,0,0,0,hrpsitnu,ncut,&nrxtnu);
 
 cout << endl;
 cout << " *** Psi analysis for the observed IceCube event times ***" << endl;
 cout << " Background event rate (Hz) : " << ratebkg << " --> psi for observed event times : " << psitnu << endl;
 cout << " ==> P-value of the observed psi : " << pvaluetnu << " Used # of randomisations : " << nrxtnu << endl;
 cout << " Extreme Psi values for the case of pure background : psimin=" << psimintnu << " psimax=" << psimaxtnu << endl;
***/

 // The Kolmogorov-Smirnov test of the Delta t distribution w.r.t. the expectation from a Poisson

/***
 cout << endl;
 cout << " *** Kolmogorov-Smirnov test results WITHOUT including normalisation analysis ***" << endl;
 hpdfbkgdt->KolmogorovTest(hdtnu,"XD");

 cout << endl;
 cout << " *** Kolmogorov-Smirnov test results WITH including normalisation analysis ***" << endl;
 hpdfbkgdt->KolmogorovTest(hdtnu,"NXD");
***/

 // My new KS test facility
 TString str="PIN";

 TH1F* hks1=new TH1F("hks1","KS pseudo experiment data",10,-1,1);
 histos.Add(hks1);
 cout << endl;
 lab.KolmogorovTest(str,hdtnu,0,&fdtbkg,nrandom,hks1,ncut);

/****
 TH1F* hks2=new TH1F("hks2","KS pseudo experiment data",10,-1,1);
 histos.Add(hks2);
 cout << endl;
 lab.KolmogorovTest(str,hdtnu,hpdfbkgdt,0,nrandom,hks2,ncut);
****/

 // Bayesian psi value analysis of the Delta t distribution
 Double_t psidtnu=-1;
 psidtnu=math.PsiValue(hdtnu,0,&fdtbkg);
//@@@ psidtnu=math.PsiValue(hdtnu,hpdfbkgdt);

 Float_t psimindtnu=-1;
 Float_t psimaxdtnu=-1;
 psimindtnu=math.PsiExtreme(hdtnu,0,&fdtbkg,-2);
 psimaxdtnu=math.PsiExtreme(hdtnu,0,&fdtbkg,-1);
//@@@ psimindtnu=math.PsiExtreme(hdtnu,hpdfbkgdt,0,-2);
//@@@ psimaxdtnu=math.PsiExtreme(hdtnu,hpdfbkgdt,0,-1);

 TH1F* hrpsidtnu=new TH1F("hrpsidtnu","Random #psi distr. for bkg hypothesis of event dt",1000,psimindtnu-1.,psimaxdtnu+1.);
 histos.Add(hrpsidtnu);

 Double_t nrxdtnu=-1; // Returned number of actually performed randomizations
 Double_t pvaluedtnu=-1;
 pvaluedtnu=math.PsiPvalue(-1,nrandom,hdtnu,0,&fdtbkg,0,0,hrpsidtnu,ncut,&nrxdtnu);
//@@@ pvaluedtnu=math.PsiPvalue(-1,nrandom,hdtnu,hpdfbkgdt,0,0,0,hrpsidtnu,ncut,&nrxdtnu);
 
 cout << endl;
 cout << " *** Psi analysis for the Poisson Dt distributions ***" << endl;
 cout << " Background event rate (Hz) : " << ratebkg << " --> psi for observed neutrino dt distr. : " << psidtnu << endl;
 cout << " ==> P-value of the observed psi : " << pvaluedtnu << " Used # of randomisations : " << nrxdtnu << endl;
 cout << " Extreme Psi values for the case of pure background : psimin=" << psimindtnu << " psimax=" << psimaxdtnu << endl;

 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 // The output file for the produced histograms
 TFile* fout=new TFile("output.root","RECREATE","GW analysis results");

 Int_t nhistos=histos.GetEntries();
 cout << endl;
 cout << " ====== The following " << nhistos << " histograms are generated in the output file ======" << endl;
 for (Int_t ih=0; ih<nhistos; ih++)
 {
  TH1* hx=(TH1*)histos.At(ih);
  if (!hx) continue;
  cout << hx->GetName() << " : " << hx->GetTitle() << endl;
  hx->Write();
 }

 fout->Write();
}
