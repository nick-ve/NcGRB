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
 Float_t xmin=-1000;
 Float_t xmax=1000;
 Int_t nbins=200000;

 Float_t brate=0.020; // The background event rate in Hz
 Float_t srate=0.000; // The signal event rate in Hz 
 Float_t tsig=300;    // The central signal arrival time
 Float_t dtsig=0.1;   // Spread of the signal arrival time
 Int_t nevtdt=2;      // Number of events within a dt cell for which the inter-event dt statistics will be performed 

 Float_t rate=brate+srate;          // Total event rate
 Int_t nevt=int(brate*(xmax-xmin)); // The number background of events
 Int_t nsig=int(srate*(xmax-xmin)); // The number of signal events

 cout << " Event rates (Hz) Background:" << rate << " Signal:" << srate << " Total:" << rate << endl;

 Double_t nrandom=1e3; // (Maximum) number of randomised configurations for direct psi P-value determination (0 means n=1e19)
 Int_t ncut=10;        // Number of obtained psi>psi0 entries at which randomisations will be terminated (0= no early termination)

 // Fill a histogram with uniform event arrival times
 // This will mimic a Poissonian background
 // Also a filling with identical time differences is possible

 TH1F* ht=new TH1F("ht","Arrival times;Arrival time;Counts",nbins,xmin,xmax);
 Double_t tbin=ht.GetXaxis()->GetBinWidth(1);
 Double_t tmin=ht.GetXaxis()->GetXmin();
 Double_t tmax=ht.GetXaxis()->GetXmax();

 Float_t step=(xmax-xmin)/float(nevt);
 Float_t x=xmin;
 for (Int_t i=0; i<nevt; i++)
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

/// TH1F hdt=lab.GetDxHistogram(ht,nevtdt,1,0,-1);
//// TH1F hdt=lab.GetDxHistogram(ht,nevtdt,0,0,-1);
///// TH1F hdt=lab.GetDxHistogram(ht,nevtdt,0,tbin,-1);
 TH1F hdt=GetDxHistogram2(ht,nevtdt,0,0,-1,0);

 Int_t nbdt=hdt.GetNbinsX();
 Double_t deltatbin=hdt.GetXaxis()->GetBinWidth(1);
 Double_t deltatmin=hdt.GetXaxis()->GetXmin();
 Double_t deltatmax=hdt.GetXaxis()->GetXmax();

/*********
 // Create title and labels for the delta t histogram
 TString title;
 TString s;
 title="Delta t distribution to contain exactly ";
 title+=nevtdt;
 title+=" events";
 title+=";#Deltat (in sec);Counts per bin of size %-10.3g";
 s=title.Format(title.Data(),deltatbin);

 hdt.SetTitle(s.Data());
***********/
 
 TCanvas* c2=new TCanvas("c2","c2");
 hdt.Draw();

 // Creation of the Poisson based dt PDF for a background only hypothesis
 TF1 fdtbkg=math.PoissonDtDist(rate,nevtdt);

/*********
 // Alternative pdf for hit-empty Deta t-hit
 TF1 fdtbkg("fdtbkg","[0]*[0]*exp(-[0]*x)");
 fdtbkg.SetParameter(0,rate);
**********/

 fdtbkg.SetRange(0,deltatmax);

 // Provide the dt PDF histogram with the same binning as the Delta t distribution
 fdtbkg.SetNpx(1000);
 TH1* hpdfbkgdt=(TH1*)fdtbkg.GetHistogram()->Clone();
 hpdfbkgdt->SetName("hpdfbkgdt");

 Float_t hdtmedian=sample.GetMedian(&hdt);
 Float_t fdtmedian=sample.GetMedian(hpdfbkgdt);

 cout << " hdt median:" << hdtmedian << " pdfmedian:" << fdtmedian << endl;

 TCanvas* c3=new TCanvas("c3","c3");
 hpdfbkgdt->Draw();

 // The Kolmogorov-Smirnov test of the Delta t distribution w.r.t. the expectation from a Poisson
 TString str="PIN";
 TH1F* hks=new TH1F("hks","KS pseudo experiment data",10,-1,1);
 cout << endl;
 lab.KolmogorovTest(str,&hdt,0,&fdtbkg,nrandom,hks,ncut);
 
 TCanvas* c4=new TCanvas("c4","c4");
 hks->Draw();

 // Bayesian psi value analysis of the Delta t distribution
 Double_t psidt=-1;
 psidt=math.PsiValue(&hdt,0,&fdtbkg);
//@@@ psidt=math.PsiValue(&hdt,hpdfbkgdt);

 Float_t psimindt=-1;
 Float_t psimaxdt=-1;
 psimindt=math.PsiExtreme(&hdt,0,&fdtbkg,-2);
 psimaxdt=math.PsiExtreme(&hdt,0,&fdtbkg,-1);
//@@@ psimindt=math.PsiExtreme(&hdt,hpdfbkgdt,0,-2);
//@@@ psimaxdt=math.PsiExtreme(&hdt,hpdfbkgdt,0,-1);

 TH1F* hrpsidt=new TH1F("hrpsidt","Random #psi distr. for bkg hypothesis of event dt",1000,psimindt-1.,psimaxdt+1.);

 Double_t nrxdt=-1; // Returned number of actually performed randomizations
 Double_t pvaluedt=-1;
 pvaluedt=math.PsiPvalue(-1,nrandom,&hdt,0,&fdtbkg,0,0,hrpsidt,ncut,&nrxdt);
//@@@ pvaluedtnu=math.PsiPvalue(-1,nrandom,&hdt,hpdfbkgdt,0,0,0,hrpsidtnu,ncut,&nrxdtnu);
 
 cout << endl;
 cout << " *** Psi analysis for the Poisson Dt distributions ***" << endl;
 cout << " Event rate (Hz) : " << rate << " --> psi for observed dt distr. : " << psidt << endl;
 cout << " ==> P-value of the observed psi : " << pvaluedt << " Used # of randomisations : " << nrxdt << endl;
 cout << " Extreme Psi values for the case of pure background : psimin=" << psimindt << " psimax=" << psimaxdt << endl;
 
 TCanvas* c5=new TCanvas("c5","c5");
 hrpsidt->Draw();

}
