/////////////////////////////////////////////////////////////
// ROOT macro to investigate the performance of the KS and //
// Bayesian psi statistical tests.                         //
//                                                         //
// For this test a certain PDF is chosen and fictative     //
// data entries are generated either following this PDF    //
// or following a completely different distribution.       //
// The two statistical tests are evaluated for different   //
// number of generated data entries to investigate the     //
// dependence of the results with multiplicity.            //
//                                                         //
// To run the macro, just invoke the command ($=prompt)    //
// $root -b -q gw-analysis.cc                              //
//                                                         //
// All produced histograms are available in the produced   //
// output file "simple-test.root".                         //
//                                                         //
// Note: This macro needs the NCFSPack library, which is   //
// available on the IIHE cluster after issuing the command //
// $ source /ice3/software/iihe/ncfs.sh                    // 
//                                                         //
// Nick van Eijndhoven, 23-sep-2018, IIHE-VUB Brussel.     //
/////////////////////////////////////////////////////////////
{
 gSystem->Load("ncfspack");
//// gROOT->LoadMacro("NcAstrolab2.cxx+");

 // Provide various mathematical tools
 NcMath math;

 // The virtual lab providing various analysis tools
 NcAstrolab lab;

 // The internal storage array for the produced histograms
 TObjArray histos;

 Int_t ndata=1000; // Number of fictative data points to be generated
 Float_t xmin=-2;  // Minimum value for data distribution
 Float_t xmax=2;   // Maximum value for data distribution
 Int_t nbins=100;  // Number of bins for the data distribution
 TString title="";

 Double_t nrandom=1e3; // (Maximum) number of randomised configurations for direct psi P-value determination (0 means n=1e19)
 Int_t ncut=0;        // Number of obtained psi>psi0 entries at which randomisations will be terminated (0= no early termination)
 Double_t nrx=-1;      // Returned number of actually performed randomizations
 Double_t pvalue=-1;

 // The PDF to be used for the backgound hypothesis.
 TF1 fbkgpdf("fbkgpdf","TMath::Gaus(x)");
//@@@ TF1 fbkgpdf("fbkgpdf","TMath::Exp(x)");
//@@@ TF1 fbkgpdf("fbkgpdf","TMath::Landau(x)");
 title=fbkgpdf.GetTitle();
 title+=";x;Function value";
 fbkgpdf.SetTitle(title.Data());

 // The PDF to be used for the signal generation
 TF1 fsigpdf("fsigpdf","TMath::Gaus(x)");
 title=fsigpdf.GetTitle();
 title+=";x;Function value";
 fsigpdf.SetTitle(title.Data());

 // Provide the backgound PDF as histogram in the output file
 // Also provide a finely binned version for checking
 fbkgpdf.SetRange(xmin,xmax);
 fbkgpdf.SetNpx(nbins);
 TH1* hpdfbkg=(TH1*)fbkgpdf.GetHistogram()->Clone();
 hpdfbkg->SetName("hpdfbkg");
 histos.Add(hpdfbkg);
 fbkgpdf.SetNpx(10000);
 TH1* hpdfbkgfine=(TH1*)fbkgpdf.GetHistogram()->Clone();
 hpdfbkgfine->SetName("hpdfbkgfine");
 histos.Add(hpdfbkgfine);

 // Provide the signal PDF as histogram in the output file
 // Also provide a finely binned version for checking
 fsigpdf.SetRange(xmin,xmax);
 fsigpdf.SetNpx(nbins);
 TH1* hpdfsig=(TH1*)fsigpdf.GetHistogram()->Clone();
 hpdfsig->SetName("hpdfsig");
 histos.Add(hpdfsig);
 fsigpdf.SetNpx(10000);
 TH1* hpdfsigfine=(TH1*)fsigpdf.GetHistogram()->Clone();
 hpdfsigfine->SetName("hpdfsigfine");
 histos.Add(hpdfsigfine);

 // Obtain  the fictative data entries from the signal pdf
 TH1F* hsig=new TH1F("hsig","Fictative data entries;x;Counts",nbins,xmin,xmax);
 hsig->FillRandom("fsigpdf",ndata);
 // Create some underflow and overflow entries
 hsig->Fill(xmin-1.);
 hsig->Fill(xmin-2.);
 hsig->Fill(xmin-3.);
 hsig->Fill(xmin-4.);
 hsig->Fill(xmin-5.);
 hsig->Fill(xmax+1.);
 hsig->Fill(xmax+2.);
 hsig->Fill(xmax+3.);
 histos.Add(hsig);

 // Provide the cumulative distributions as histograms in the output file
 TH1F hpdfbkgfineCDH=lab.GetCumulHistogram(hpdfbkgfine,"hpdfbkgfineCDH","NB");
 histos.Add(&hpdfbkgfineCDH);
 TH1F hpdfbkgCDH=lab.GetCumulHistogram(&fbkgpdf,"hpdfbkgCDH",nbins,xmin,xmax,"FN");
 histos.Add(&hpdfbkgCDH);
 TH1F hsigCDH=lab.GetCumulHistogram(hsig,"hsigCDH","FN");
 histos.Add(&hsigCDH);
 
 ////////////////////////////////////////////////////////////////////////////////
 // Statistical evaluation of the signal and background observations           //
 //                                                                            //
 // Determination of the Bayesian psi value assumes that there is no signal.   //
 // This corresponds to searching out the Bernoulli class B_m                  //
 // with m=nbins of the histogram.                                             //
 ////////////////////////////////////////////////////////////////////////////////

 // Perform the KS statistics test on the fictative data and the background PDF
 // Provide results for all options supported by ROOT

 cout << endl;
 cout << " *** Kolmogorov-Smirnov test results WITHOUT under/overflow ***" << endl;
 hpdfbkg->KolmogorovTest(hsig,"XDN");

 cout << endl;
 cout << " *** Kolmogorov-Smirnov test results WITH under/overflow ***" << endl;
 hpdfbkg->KolmogorovTest(hsig,"UOXDN");

 // Perform my own KS test
 TString str="PIN";

 TH1F* hks1=new TH1F("hks1","KS pseudo experiment data",10,-1,1);
 histos.Add(hks1);
 cout << endl;
 lab.KolmogorovTest(str,hsig,0,&fbkgpdf,nrandom,hks1,ncut,&nrx);

 str+="UO";
 TH1F* hks2=new TH1F("hks2","KS pseudo experiment data",10,-1,1);
 histos.Add(hks2);
 cout << endl;
 lab.KolmogorovTest(str,hsig,0,&fbkgpdf,nrandom,hks2,ncut,&nrx);

 str="KIN";
 cout << endl;
 lab.KolmogorovTest(str,hsig,hpdfbkg);

 cout << endl;
 lab.KolmogorovTest(str,hsig,0,&fbkgpdf);

 str+="UO";
 cout << endl;
 lab.KolmogorovTest(str,hsig,hpdfbkg);

 cout << endl;
 lab.KolmogorovTest(str,hsig,0,&fbkgpdf);

 str="KI";
 cout << endl;
 lab.KolmogorovTest(str,hsig,hpdfbkg);

 cout << endl;
 lab.KolmogorovTest(str,hsig,0,&fbkgpdf);

 str+="UO";
 cout << endl;
 lab.KolmogorovTest(str,hsig,hpdfbkg);

 cout << endl;
 lab.KolmogorovTest(str,hsig,0,&fbkgpdf);

 // Perform the Bayesian Psi analysis on the fictative data and the background PDF 
 Double_t psi=math.PsiValue(hsig,hpdfbkg,0);
// Double_t psi=math.PsiValue(hsig,0,&fbkgpdf);

 Float_t psimin=-1;
 Float_t psimax=-1;
//@@@ psimin=math.PsiExtreme(hsig,0,&fbkgpdf,-2);
//@@@ psimax=math.PsiExtreme(hsig,0,&fbkgpdf,-1);
 psimin=math.PsiExtreme(hsig,hpdfbkg,0,-2);
 psimax=math.PsiExtreme(hsig,hpdfbkg,0,-1);

 TH1F* hrpsi=new TH1F("hrpsi","Random #psi distr. for bkg hypothesis of the fictative data",100,psimin-1.,psimax+1.);
 histos.Add(hrpsi);

//@@@ pvalue=math.PsiPvalue(-1,nrandom,hsig,0,&fbkgpdf,0,0,hrpsi,ncut,&nrx);
 pvalue=math.PsiPvalue(-1,nrandom,hsig,hpdfbkg,0,0,0,hrpsi,ncut,&nrx);
 
 cout << endl;
 cout << " *** Psi analysis for the observed fictative data ***" << endl;
 cout << " Psi w.r.t. background PDF for the observed fictative data : " << psi << endl;
 cout << " ==> P-value of the observed psi : " << pvalue << " Used # of randomisations : " << nrx << endl;
 cout << " Extreme Psi values for the case of pure background : psimin=" << psimin << " psimax=" << psimax << endl;


 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 // The output file for the produced histograms
 TFile* fout=new TFile("simple-test.root","RECREATE","Statistics test");

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
