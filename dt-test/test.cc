/////////////////////////////////////////////////////////////////////
// ROOT macro to test Dt spectra and the corresponding Erlang pdf. //
//                                                                 //
// Nick van Eijndhoven, IIHE-VUB Brussel, May 5, 2020  12:09       //
/////////////////////////////////////////////////////////////////////
{
 gSystem->Load("ncfspack");

 gROOT->LoadMacro("code.cxx+");

 // The virtual lab providing various analysis tools
 NcAstrolab lab;

 NcRandom ran;
 NcMath math;
 NcSample sample;

 // Arrival time window in seconds
 Float_t xmin=-500;
 Float_t xmax=500;
 Int_t nbins=100000;

 // Settings for the Dt histogram
 Double_t deltatbin=-1;
 Double_t deltatmin=-1;
 Double_t deltatmax=-1;
 Double_t fact=1;
 Int_t nevtdt=2; // Number of events within a dt cell for which the inter-event dt statistics will be performed 

 // The arival time window histogram
 TH1F* ht=new TH1F("ht","Arrival times;Arrival time;Counts",nbins,xmin,xmax);
 Double_t tbin=ht.GetXaxis()->GetBinWidth(1);
 Double_t tmin=ht.GetXaxis()->GetXmin();
 Double_t tmax=ht.GetXaxis()->GetXmax();

 // Background and signal rates
 Float_t brate=0.030; // The background event rate in Hz
 Float_t srate=0.005; // The signal event rate in Hz 
 Float_t tsig=0;      // The central signal arrival time
 Float_t dtsig=tbin*5.;  // Spread of the signal arrival time

 Float_t rate=brate+srate;                  // Total event rate
 Int_t nbkg=TMath::Nint(brate*(xmax-xmin)); // The number background of events
 Int_t nsig=TMath::Nint(srate*(xmax-xmin)); // The number of signal events

 cout << " Event rates (Hz) Background:" << brate << " Signal:" << srate << " Total:" << rate << endl;

 Double_t nrandom=1e3; // (Maximum) number of randomised configurations for direct psi P-value determination (0 means n=1e19)
 Int_t ncut=10;        // Number of obtained psi>psi0 entries at which randomisations will be terminated (0= no early termination)

 // Fill the time window histogram with uniform background event arrival times
 // This will mimic a Poissonian background
 // Also a filling with identical time differences is possible

 Float_t step=0;
 if (nbkg) step=(xmax-xmin)/float(nbkg);
 Float_t x=xmin;
 for (Int_t i=0; i<nbkg; i++)
 {
  x=ran.Uniform(xmin,xmax); // Uniform filling
///  x=x+step; // Identical dt filling
  ht->Fill(x);
 }

 // Add some signal(s)
 for (Int_t i=0; i<nsig; i++)
 {
  x=ran.Uniform(tsig-dtsig,tsig+dtsig);
  ht->Fill(x);
 }
 
 TCanvas* c1=new TCanvas("c1","c1");
 ht->Draw();
 
 // Create the delta t histogram
 TH1F hdt=lab.GetDxHistogram(ht,nevtdt,deltatbin,deltatmin,deltatmax);
 TH1F hdt0=GetDxHistogram4(ht,nevtdt,deltatbin,deltatmin,deltatmax,0,fact);
 TH1F hdt1=GetDxHistogram4(ht,nevtdt,deltatbin,deltatmin,deltatmax,1,fact);
 TH1F hdt2=GetDxHistogram4(ht,nevtdt,deltatbin,deltatmin,deltatmax,2,fact);
 TH1F hdt3=GetDxHistogram4(ht,nevtdt,deltatbin,deltatmin,deltatmax,3,fact);
 
 TCanvas* c2=new TCanvas("c2","hdt");
 hdt.Draw();
 
 TCanvas* c20=new TCanvas("c20","hdt0");
 hdt0.Draw();
 
 TCanvas* c21=new TCanvas("c21","hdt1");
 hdt1.Draw();
 
 TCanvas* c22=new TCanvas("c22","hdt2");
 hdt2.Draw();
 
 TCanvas* c23=new TCanvas("c23","hdt3");
 hdt3.Draw();

 TH1F* hdtx=&hdt1;

 Int_t nbdt=hdtx->GetNbinsX();
 deltatbin=hdtx->GetXaxis()->GetBinWidth(1);
 deltatmin=hdtx->GetXaxis()->GetXmin();
 deltatmax=hdtx->GetXaxis()->GetXmax();

 cout << " dtnbins:" << nbdt << " dtbsize:" << deltatbin << " dtmin:" << deltatmin << " dtmax:" << deltatmax << endl;

 // Creation of the Poisson based dt PDF for a background only hypothesis
 TF1 fdtbkg=math.PoissonDtDist(rate,nevtdt);
 fdtbkg.SetRange(deltatmin,deltatmax);

 // Provide the dt PDF histogram with the same binning as the Delta t distribution
 Int_t npx=1000;
 npx=nbdt;
 fdtbkg.SetNpx(npx);
 TH1* hpdfbkgdt=(TH1*)fdtbkg.GetHistogram()->Clone();
 hpdfbkgdt->SetName("hpdfbkgdt");

 Float_t hdtmedian=sample.GetMedian(hdtx);
 Float_t fdtmedian=sample.GetMedian(hpdfbkgdt);

 cout << " hdtx median:" << hdtmedian << " pdfmedian:" << fdtmedian << endl;

 TCanvas* c3=new TCanvas("c3","c3");
 hpdfbkgdt->Draw();

 // The Kolmogorov-Smirnov test of the Delta t distribution w.r.t. the expectation from a Poisson
 TString str="PIN";
 TH1F* hks=new TH1F("hks","KS pseudo experiment data",10,-1,1);
 cout << endl;
 lab.KolmogorovTest(str,hdtx,0,&fdtbkg,nrandom,hks,ncut);
 
 TCanvas* c4=new TCanvas("c4","c4");
 hks->Draw();

 // Bayesian psi value analysis of the Delta t distribution
 Double_t psidt=-1;
 psidt=math.PsiValue(hdtx,0,&fdtbkg);
//@@@ psidt=math.PsiValue(&hdt,hpdfbkgdt);
 
 Float_t psimindt=-1;
 Float_t psimaxdt=-1;
 psimindt=math.PsiExtreme(hdtx,0,&fdtbkg,-2);
 psimaxdt=math.PsiExtreme(hdtx,0,&fdtbkg,-1);
//@@@ psimindt=math.PsiExtreme(hdtx,hpdfbkgdt,0,-2);
//@@@ psimaxdt=math.PsiExtreme(hdtx,hpdfbkgdt,0,-1);

 TH1F* hrpsidt=new TH1F("hrpsidt","Random #psi distr. for bkg hypothesis of event dt",1000,psimindt-1.,psimaxdt+1.);

 Double_t nrxdt=-1; // Returned number of actually performed randomizations
 Double_t pvaluedt=-1;
 pvaluedt=math.PsiPvalue(-1,nrandom,hdtx,0,&fdtbkg,0,0,hrpsidt,ncut,&nrxdt);
//@@@ pvaluedtnu=math.PsiPvalue(-1,nrandom,hdtx,hpdfbkgdt,0,0,0,hrpsidtnu,ncut,&nrxdtnu);
 
 cout << endl;
 cout << " *** Psi analysis for the Poisson Dt distributions ***" << endl;
 cout << " Event rate (Hz) : " << rate << " --> psi for observed dt distr. : " << psidt << endl;
 cout << " ==> P-value of the observed psi : " << pvaluedt << " Used # of randomisations : " << nrxdt << endl;
 cout << " Extreme Psi values for the case of pure background : psimin=" << psimindt << " psimax=" << psimaxdt << endl;
 
 TCanvas* c5=new TCanvas("c5","c5");
 hrpsidt->Draw();
}
