// Code to test new implementation of Dt analysis
#include "NcAstrolab.h"
///////////////////////////////////////////////////////////////////////////
TH1F GetDxHistogram2(TH1* hx,Int_t nc,Double_t dxbin,Double_t dxmin,Double_t dxmax,Int_t mode)
{
// Provide the interval size (dx) distribution of X-axis intervals between a certain
// fixed amount of consecutive histogram entries of the specified input histogram.
// This facility can for instance be used to investigate the distribution
// of time intervals between observed events.
//
// Input arguments :
// -----------------
// hx    : The input histogram
// nc    : The step count to arrive at the required consecutive entry (see example below)
// dxbin : The bin size of the X-axis for the dx distribution (see also the note below)
//         dxbin=0  ==> Bin size taken the same as the input histogram hx
//         dxbin=-1 ==> Bin size taken to be the minimal encountered dx interval
//                     (or hx bin size in case the minimal encountered dx was 0) 
//         dxbin=-2 ==> Bin size taken to be nc times the bin size of the input histogram hx
// dxmin : The lower bound of the produced histogram
//         dxmin<0 ==> Lower bound taken to be the minimal encountered dx interval           
// dxmax : The upper bound of the produced histogram
//         dxmax<0 ==> Upper bound taken to be the maximal encountered dx interval
//                     increased with one additional bin size to contain the maximal
//                     encountered dx interval value in the output histogram
// mode  : 0 ==> Bin contents are regarded as event counts (see also note 1) below)
//               Multiple entries in the same bin are treated as multiple events
//         1 ==> Bin contents are treated as (weighted) values
//               Each filled bin is treated as a single event
//
// Returned object : The 1-D histogram (TH1F) containing the dx distribution.
//
// Default values : dxbin=-1, dxmin=-1, dxmax=-1 and mode=0. 
//
// Note : A too coarse bin size "dxbin" may lead to binning effects which may affect
//        a statistical interpretation of the resulting histogram, for instance in the
//        light of a comparison with known background distributions.
//        In case the input histogram "hx" reflects an unbinned situation, i.e. not more
//        than one single entry per bin, then "dxbin=0" assures that also the produced
//        "dx histogram" will be free of binning effects.
//        For sparsely populated histograms, the same may hold for "dxbin=-2", which has
//        the advantage of having fewer bins which may reduce computing time when
//        performing statistical analyses (e.g. randomised combinatorics) on it.
//        It is advised to always first inspect the produced "dx histogram" before
//        performing a statistics analysis on it.
//
// Example :
// ---------
// Histogram "hx" contains the counts of event times according to a Poisson distribution.
// Specifying nc=1 will provide the histogram of dt intervals between each consecutive event,
// i.e. the distribution of time intervals between the events (1,2), (2,3), (3,4) etc.
// Specifying nc=2 will provide the histogram of dt intervals between every 2nd consecutive event,
// i.e. the distribution of time intervals between the events (1,3), (2,4), (3,5) etc.
// In this case of a Poissonian input histogram these produced dt distributions are known
// as the Erlang distributions (see also class NcMath), for which each time interval contains
// exactly "nc" events of which the last event occurs right at the end of the interval.
//
// Notes :
// -------
// 1) In case a certain bin of the input histogram hx contains several entries, all these
//    entries are by default treated as separate entities with exactly the same x-value.
//    Consequently, this may result in (several) dx intervals of value 0, resulting in
//    entries at 0 in the output histogram.
//    In case one wants to avoid this effect, the binning of the input histogram should be
//    fine enough to reflect a basically unbinned situation.
//    Another way of avoiding this effect is to specify "mode=1", which will treat every filled
//    bin as a single event, irrespective of the bin content. This "mode=1" also allows the
//    treatment of input histograms for which the bin contents do not represent event counts,
//    e.g. histograms which have been rescaled or filled with weights.
// 2) In case dxmax<0 the actual bin size of the output histogram (slightly) differs from the
//    provided one because of the addition of one bin size to the auto-detected dxmax.
//    If one wants multiple output histograms resulting from multiple invokations of this
//    function to have identical range and bin size, the user is advised to obtain the output
//    histogram parameters from the first produced histogram and pass these to the subsequent invokations. 
// 3) This member function is used recursively.

 TH1F hdx;

 if (!hx) return hdx;
 
 Int_t nenhx=hx->GetEntries();
 if (!nenhx) return hdx;

 if (!nc) return hdx;

 // Create the output histogram if all parameters have been specified or determined automatically.
 // If not, this will be done at a recursive invokation (see below) once (some of) the
 // parameters have been determined automatically from the input histogram.
 if (dxmin>=0 && dxmax>=dxmin && dxbin>0)
 {
  Int_t nbins=int((dxmax-dxmin)/dxbin);
  hdx.SetBins(nbins,dxmin,dxmax);

  // Add histogram and axes titles
  TString s;

  s="Dx interval distribution between ";
  s+=nc+1;
  s+=" consecutive entries (nc=";
  s+=nc;
  s+=")";
  hdx.SetNameTitle("DxHistogram",s.Data());

  Double_t binwidth=hdx.GetXaxis()->GetBinWidth(1);
  s.Form("Counts per %-10.3g",binwidth);
  hdx.GetXaxis()->SetTitle("Dx interval");
  hdx.GetYaxis()->SetTitle(s.Data());
 }

 // Determine the minimum and maximum encountered dx or fill the output histogram
 Double_t x1,x2,deltax;
 Int_t nx1,nx2;
 Double_t deltaxmin=0;
 Double_t deltaxmax=0;
 Bool_t first=kTRUE;
 Int_t ndxcount=0;
 Int_t jstart;

 Int_t nbhx=hx->GetNbinsX();
 Double_t value=0;
 for (Int_t i=1; i<=nbhx; i++)
 {
  deltax=-1;
  ndxcount=0;
  value=hx->GetBinContent(i);
  nx1=0;
  if (value) nx1=1;
  if (!mode) nx1=int(value+1e-3);
  x1=hx->GetBinCenter(i);

  while (nx1)
  {
   // Check for multiple counts (left) in this bin
   jstart=i+1;
   if (nx1>1) jstart=i;

   for (Int_t j=jstart; j<=nbhx; j++)
   {
    value=hx->GetBinContent(j);
    nx2=0;
    if (value) nx2=1;
    if (!mode) nx2=int(value+1e-3);
    x2=hx->GetBinCenter(j);

    if (j==i) nx2=nx1-1; // Counting within the same bin

    // Empty bin
    if (nx2<1) continue;

    ndxcount+=nx2;

    if (ndxcount>=nc)
    {
     deltax=x2-x1;
     if (dxmin>=0 && dxmax>=dxmin && dxbin>0) // Output histogram has been initialised
     {
      hdx.Fill(deltax);
     }
     else // Auto-determination of the output histogram range
     {
      if (first || deltax<deltaxmin) deltaxmin=deltax;
      if (first || deltax>deltaxmax) deltaxmax=deltax;
      first=kFALSE;
     }
     ndxcount=0;
     break;
    }
   }
   nx1--;
  }
 }

 // Check if a recursive call is needed to actually create and fill the output histogram
 Int_t nen=hdx.GetEntries();
 if (!nen)
 {
  // Set the bin size (if needed) for the output histogram
  if (!dxbin) dxbin=hx->GetBinWidth(1);
  if (dxbin==-1)
  {
   dxbin=hx->GetBinWidth(1);
   if (deltaxmin>0) dxbin=deltaxmin;
  }
  if (dxbin==-2)
  {
   dxbin=hx->GetBinWidth(1);
   dxbin=dxbin*float(nc);
  }

  // Set the auto-determined range of the output histogram
  if (dxmin<0) dxmin=deltaxmin;
  if (dxmax<0) dxmax=deltaxmax+dxbin;

  // Invoke the recursive call to create and fill the output histogram
  hdx=GetDxHistogram2(hx,nc,dxbin,dxmin,dxmax,mode);
 }

 return hdx;
}
///////////////////////////////////////////////////////////////////////////
TH1F GetDxHistogram3(TH1* hx,Int_t nc,Double_t dxbin,Double_t dxmin,Double_t dxmax,Int_t mode)
{
// Provide the interval size (dx) distribution of X-axis intervals between a certain
// fixed amount of consecutive histogram entries of the specified input histogram.
// This facility can for instance be used to investigate the distribution
// of time intervals between observed events.
//
// Input arguments :
// -----------------
// hx    : The input histogram
// nc    : The step count to arrive at the required consecutive entry (see example below)
// dxbin : The bin size of the X-axis for the dx distribution (see also the note below)
//         dxbin=0  ==> Bin size taken the same as the input histogram hx
//         dxbin=-1 ==> Bin size taken to be the minimal encountered dx interval
//                     (or hx bin size in case the minimal encountered dx was 0) 
//         dxbin=-2 ==> Bin size taken to be nc times the bin size of the input histogram hx
// dxmin : The lower bound of the produced histogram
//         dxmin<0 ==> Lower bound taken to be the minimal encountered dx interval           
// dxmax : The upper bound of the produced histogram
//         dxmax<0 ==> Upper bound taken to be the maximal encountered dx interval
//                     increased with one additional bin size to contain the maximal
//                     encountered dx interval value in the output histogram
// mode  : 0 ==> Bin contents are regarded as event counts (see also note 1) below)
//               Multiple entries in the same bin are treated as multiple events with exactly
//               the same x-value, possibly resulting in interval values dx=0.
//         1 ==> Same as "mode=0", but when an interval value dx=0 is encountered, it will be
//               converted to a random fraction of the corresponding bin size (see note 1) below).  
//         2 ==> Bin contents are treated as (weighted) values
//               Each filled bin is treated as a single event
//
// Returned object : The 1-D histogram (TH1F) containing the dx distribution.
//
// Default values : dxbin=-1, dxmin=-1, dxmax=-1 and mode=1. 
//
// Note : A too coarse bin size "dxbin" may lead to binning effects which may affect
//        a statistical interpretation of the resulting histogram, for instance in the
//        light of a comparison with known background distributions.
//        In case the input histogram "hx" reflects an unbinned situation, i.e. not more
//        than one single entry per bin, then "dxbin=0" assures that also the produced
//        "dx histogram" will be free of binning effects.
//        For sparsely populated histograms, the same may hold for "dxbin=-2", which has
//        the advantage of having fewer bins which may reduce computing time when
//        performing statistical analyses (e.g. randomised combinatorics) on it.
//        It is advised to always first inspect the produced "dx histogram" before
//        performing a statistics analysis on it.
//
// Example :
// ---------
// Histogram "hx" contains the counts of event times according to a Poisson distribution.
// Specifying nc=1 will provide the histogram of dt intervals between each consecutive event,
// i.e. the distribution of time intervals between the events (1,2), (2,3), (3,4) etc.
// Specifying nc=2 will provide the histogram of dt intervals between every 2nd consecutive event,
// i.e. the distribution of time intervals between the events (1,3), (2,4), (3,5) etc.
// In this case of a Poissonian input histogram these produced dt distributions are known
// as the Erlang distributions (see also class NcMath), for which each time interval contains
// exactly "nc" events of which the last event occurs right at the end of the interval.
//
// Notes :
// -------
// 1) In case a certain bin of the input histogram hx contains several entries, all these
//    entries are regarded as separate entities with exactly the same x-value.
//    Consequently, this may result in (several) dx intervals of value 0, resulting in
//    entries at 0 in the output histogram.
//    In case one wants to avoid this effect, the binning of the input histogram could be
//    chosen fine enough to reflect a basically unbinned situation.
//    Another way to avoid this effect is to specify "mode=1", for which the corresponding
//    dx=0 values will be converted into random fractions of the corresponding bin size.
//    Yet another way of avoiding this effect is to specify "mode=2", which will treat every filled
//    bin as a single event, irrespective of the bin content. This "mode=2" allows the
//    treatment of input histograms for which the bin contents do not represent event counts,
//    e.g. histograms which have been rescaled or filled with weights.
// 2) In case dxmax<0 the actual bin size of the output histogram (slightly) differs from the
//    provided one because of the addition of one bin size to the auto-detected dxmax.
//    If one wants multiple output histograms resulting from multiple invokations of this
//    function to have identical range and bin size, the user is advised to obtain the output
//    histogram parameters from the first produced histogram and pass these to the subsequent invokations. 
// 3) This member function is used recursively.

 TH1F hdx;

 if (mode<0 || mode>2) return hdx;

 if (!hx) return hdx;
 
 Int_t nenhx=hx->GetEntries();
 if (!nenhx) return hdx;

 if (!nc) return hdx;

 // Create the output histogram if all parameters have been specified or determined automatically.
 // If not, this will be done at a recursive invokation (see below) once (some of) the
 // parameters have been determined automatically from the input histogram.
 if (dxmin>=0 && dxmax>=dxmin && dxbin>0)
 {
  Int_t nbins=int((dxmax-dxmin)/dxbin);
  hdx.SetBins(nbins,dxmin,dxmax);

  // Add histogram and axes titles
  TString s;

  s="Dx interval distribution between ";
  s+=nc+1;
  s+=" consecutive entries (nc=";
  s+=nc;
  s+=", mode=";
  s+=mode;
  s+=")";
  hdx.SetNameTitle("DxHistogram",s.Data());

  Double_t binwidth=hdx.GetXaxis()->GetBinWidth(1);
  s.Form("Counts per %-10.3g",binwidth);
  hdx.GetXaxis()->SetTitle("Dx interval");
  hdx.GetYaxis()->SetTitle(s.Data());
 }

 NcRandom ran;

 // Determine the minimum and maximum encountered dx or fill the output histogram
 Double_t x1,x2,deltax;
 Int_t nx1,nx2;
 Double_t deltaxmin=0;
 Double_t deltaxmax=0;
 Bool_t first=kTRUE;
 Int_t ndxcount=0;
 Int_t jstart;

 Int_t nbhx=hx->GetNbinsX();
 Double_t value=0;
 Double_t bsize=0;
 for (Int_t i=1; i<=nbhx; i++)
 {
  deltax=-1;
  ndxcount=0;
  bsize=hx->GetBinWidth(i);
  value=hx->GetBinContent(i);
  nx1=0;
  if (value) nx1=1;
  if (mode<2) nx1=int(value+1e-3);
  x1=hx->GetBinCenter(i);

  while (nx1)
  {
   // Check for multiple counts (left) in this bin
   jstart=i+1;
   if (nx1>1) jstart=i;

   for (Int_t j=jstart; j<=nbhx; j++)
   {
    value=hx->GetBinContent(j);
    nx2=0;
    if (value) nx2=1;
    if (mode<2) nx2=int(value+1e-3);
    x2=hx->GetBinCenter(j);

    if (j==i) nx2=nx1-1; // Counting within the same bin

    // Empty bin
    if (nx2<1) continue;

    ndxcount+=nx2;

    if (ndxcount>=nc)
    {
     deltax=x2-x1;
     if (dxmin>=0 && dxmax>=dxmin && dxbin>0) // Output histogram has been initialised
     {
      while (mode==1 && deltax<=0)
      {
       deltax=ran.Uniform(0.,bsize);
      }
      hdx.Fill(deltax);
     }
     else // Auto-determination of the output histogram range
     {
      if (first || deltax<deltaxmin) deltaxmin=deltax;
      if (first || deltax>deltaxmax) deltaxmax=deltax;
      first=kFALSE;
     }
     ndxcount=0;
     break;
    }
   }
   nx1--;
  }
 }

 // Check if a recursive call is needed to actually create and fill the output histogram
 Int_t nen=hdx.GetEntries();
 if (!nen)
 {
  // Set the bin size (if needed) for the output histogram
  if (!dxbin) dxbin=hx->GetBinWidth(1);
  if (dxbin==-1)
  {
   dxbin=hx->GetBinWidth(1);
   if (deltaxmin>0) dxbin=deltaxmin;
  }
  if (dxbin==-2)
  {
   dxbin=hx->GetBinWidth(1);
   dxbin=dxbin*float(nc);
  }

  // Set the auto-determined range of the output histogram
  if (dxmin<0) dxmin=deltaxmin;
  if (dxmax<0) dxmax=deltaxmax+dxbin;

  // Invoke the recursive call to create and fill the output histogram
  hdx=GetDxHistogram3(hx,nc,dxbin,dxmin,dxmax,mode);
 }

 return hdx;
}
///////////////////////////////////////////////////////////////////////////
TH1F GetDxHistogram4(TH1* hx,Int_t nc,Double_t dxbin,Double_t dxmin,Double_t dxmax,Int_t mode,Double_t fact)
{
// Provide the interval size (dx) distribution of X-axis intervals between a certain
// fixed amount of consecutive histogram entries of the specified input histogram.
// This facility can for instance be used to investigate the distribution
// of time intervals between observed events.
//
// Input arguments :
// -----------------
// hx    : The input histogram.
// nc    : The step count to arrive at the required consecutive entry (see example below).
// dxbin : The bin size of the X-axis for the dx distribution (see also the note below).
//         dxbin=0  ==> Bin size taken the same as the input histogram "hx".
//         dxbin=-1 ==> Bin size taken to be the minimal encountered dx interval,
//                      or "fact" times the bin size of the input histogram "hx", if the latter
//                      is the larger. This provides a protection against too small bin sizes,
//                      and consequently a very large number of bins, in case an extremely small
//                      dx interval value is encountered.
//                      In case "fact<=0" and no dx>0 interval was encountered, the bin size of
//                      the input histogram "hx" is taken. 
//         dxbin=-2 ==> Bin size taken to be "nc" times the bin size of the input histogram "hx".
// dxmin : The lower bound of the produced histogram.
//         dxmin<0 ==> Lower bound taken to be the minimal encountered dx interval.
// dxmax : The upper bound of the produced histogram.
//         dxmax<0 ==> Upper bound taken to be the maximal encountered dx interval
//                     increased with one additional bin size to contain the maximal
//                     encountered dx interval value in the output histogram.
// mode  : 0 ==> Bin contents of the input histogram "hx" are regarded as event counts
//               (rounded to the nearest integer) with as x-value the center of the bin.
//               Multiple entries in the same bin are treated as multiple events with exactly
//               the same x-value, possibly resulting in interval values dx=0.
//               See also note 1) below.
//         1 ==> Same as "mode=0", but the x-value of each entry of the input histogram "hx" will be
//               assigned to a random value within the corresponding bin. See note 1) below.  
//         2 ==> Bin contents of the input histogram "hx" are treated as (weighted) values.
//               Each filled bin is treated as a single event with as x-value the center of the bin.
//         3 ==> Same as "mode=2", but the x-value of each entry of the input histogram "hx" will be
//               assigned to a random value within the corresponding bin. See also note 1) below.
// fact  : Multiplication factor applied to the bin size of the input histogram "hx" to provide
//         a lower limit to the bin size of the produced Dx histogram when "dxbin=-1". See note 1) below.
//
// Returned object : The 1-D histogram (TH1F) containing the dx distribution.
//
// Default values : dxbin=-1, dxmin=-1, dxmax=-1, mode=1 and fact=1.
//
// Note : A too coarse bin size "dxbin" may lead to binning effects which may affect
//        a statistical interpretation of the resulting histogram, for instance in the
//        light of a comparison with known background distributions.
//        In case the input histogram "hx" reflects an unbinned situation, i.e. not more
//        than one single entry per bin, then "dxbin=0" assures that also the produced
//        "dx histogram" will be free of binning effects.
//        For sparsely populated histograms, the same may hold for "dxbin=-2", which has
//        the advantage of having fewer bins which may reduce computing time when
//        performing statistical analyses (e.g. randomised combinatorics) on it.
//        It is advised to always first inspect the produced "dx histogram" before
//        performing a statistics analysis on it.
//
// Example :
// ---------
// Histogram "hx" contains the counts of event times according to a Poisson distribution.
// Specifying "nc=1" will provide the histogram of dt intervals between each consecutive event,
// i.e. the distribution of time intervals between the events (1,2), (2,3), (3,4) etc.
// Specifying "nc=2" will provide the histogram of dt intervals between every 2nd consecutive event,
// i.e. the distribution of time intervals between the events (1,3), (2,4), (3,5) etc.
// In this case of a Poissonian input histogram these produced dt distributions are known
// as the Erlang distributions (see also class NcMath), for which each time interval contains
// exactly "nc" events of which the last event occurs right at the end of the interval.
//
// Notes :
// -------
// 1) In case a certain bin of the input histogram "hx" contains several entries, all these
//    entries are regarded as separate entities with exactly the same x-value.
//    Consequently, this may result in (several) dx intervals of value 0, resulting in
//    entries at 0 in the output histogram.
//    In case one wants to avoid this effect, the binning of the input histogram "hx" could be
//    chosen fine enough to reflect a basically unbinned situation.
//    Another way to avoid this effect is to specify "mode=1", for which all entries in the
//    input histogram "hx" will be assigned a random value within the corresponding bin, instead of
//    the usual center of bin value. The "sensitivity" to this randomisation effect can be tuned
//    by the input argument "fact" which controls the minimal bin size of the produced Dx histogram.
//    Yet another way of avoiding this effect is to specify "mode=2" or "mode=3", which will treat
//    every filled bin of the input histogram as a single event, irrespective of the bin content.
//    These modes 2 or 3 allow the treatment of input histograms for which the bin contents do
//    not represent event counts, e.g. histograms which have been rescaled or filled with weights.
// 2) In case "dxmax<0" the actual bin size of the output histogram (slightly) differs from the
//    provided one because of the addition of one or two bin size(s) to the auto-detected dxmax
//    in order to prevent entries to fall outside the histogram range. 
//    If one wants multiple output histograms resulting from multiple invokations of this function
//    to have identical range and bin size, the user is advised to obtain the output histogram
//    parameters from the first produced histogram and pass these to the subsequent invokations. 
// 3) This member function is used recursively.

 TH1F hdx;

 if (mode<0 || mode>3) return hdx;

 if (!hx) return hdx;

 if (nc<1) return hdx;
 
 Int_t nenhx=hx->GetEntries();
 if (nenhx<=nc) return hdx;

 Int_t idxbin=TMath::Nint(dxbin);
 if (idxbin<-2) return hdx;

 // Create the output histogram if all parameters have been specified or determined automatically.
 // If not, this will be done at a recursive invokation (see below) once (some of) the
 // parameters have been determined automatically from the input histogram.
 if (dxmin>=0 && dxmax>=dxmin && dxbin>0)
 {
  Int_t nbins=1;
  Double_t range=dxmax-dxmin;
  if (range>dxbin) nbins=TMath::Nint(range/dxbin);
  hdx.SetBins(nbins,dxmin,dxmax);

  // Add histogram and axes titles
  TString s;

  s="Dx interval distribution between ";
  s+=nc+1;
  s+=" consecutive entries (nc=";
  s+=nc;
  s+=", mode=";
  s+=mode;
  s+=")";
  hdx.SetNameTitle("DxHistogram",s.Data());

  Double_t binwidth=hdx.GetXaxis()->GetBinWidth(1);
  s.Form("Counts per bin of size %-10.3g",binwidth);
  hdx.GetXaxis()->SetTitle("Dx interval");
  hdx.GetYaxis()->SetTitle(s.Data());
 }

 NcRandom ran;

 // Determine the minimum and maximum encountered dx or fill the output histogram
 Double_t x1,x2,deltax;
 Int_t nx1,nx2;
 Double_t deltaxmin=0;
 Double_t deltaxmax=0;
 Bool_t found=kFALSE;
 Int_t ndxcount=0;
 Int_t jstart;

 Int_t nbhx=hx->GetNbinsX();
 Double_t value=0;
 Double_t xlow=0;
 Double_t xup=0;
 Double_t bsize=0;
 for (Int_t i=1; i<=nbhx; i++)
 {
  deltax=-1;
  ndxcount=0;
  xlow=hx->GetBinLowEdge(i);
  bsize=hx->GetBinWidth(i);
  xup=xlow+bsize;
  x1=hx->GetBinCenter(i);
  if (mode==1 || mode==3) x1=ran.Uniform(xlow,xup);
  value=hx->GetBinContent(i);
  nx1=0;
  if (value) nx1=1;
  if (mode<2) nx1=TMath::Nint(value);

  while (nx1>0)
  {
   // Check for multiple counts (left) in this bin
   jstart=i+1;
   if (nx1>1) jstart=i;

   for (Int_t j=jstart; j<=nbhx; j++)
   {
    xlow=hx->GetBinLowEdge(j);
    bsize=hx->GetBinWidth(j);
    xup=xlow+bsize;
    x2=hx->GetBinCenter(j);
    if (mode==1 || mode==3) x2=ran.Uniform(xlow,xup);
    value=hx->GetBinContent(j);
    nx2=0;
    if (value) nx2=1;
    if (mode<2) nx2=TMath::Nint(value);

    if (j==i) nx2=nx1-1; // Counting within the same bin

    // Empty bin
    if (nx2<1) continue;

    ndxcount+=nx2;

    if (ndxcount>=nc)
    {
     deltax=fabs(x2-x1);
     if (dxmin>=0 && dxmax>=dxmin && dxbin>0) // Output histogram has been initialised
     {
      hdx.Fill(deltax);
     }
     else // Auto-determination of the output histogram range
     {
       if (!found || deltax<deltaxmin) deltaxmin=deltax;
       if (!found || deltax>deltaxmax) deltaxmax=deltax;
     }
     ndxcount=0;
     found=kTRUE;
     break;
    }
   }
   nx1--;
  }
 }

 // Check if a suitable configuration of entries was encountered
 if (!found) return hdx;

 // Check if a recursive call is needed to actually create and fill the output histogram
 Int_t nen=hdx.GetEntries();
 if (!nen)
 {
  // Set the bin size (if needed) for the output histogram
  if (!idxbin) dxbin=hx->GetBinWidth(1);
  if (idxbin==-1)
  {
   dxbin=(hx->GetBinWidth(1))*fact;
   if (deltaxmin>0 && deltaxmin>dxbin) dxbin=deltaxmin;
   if (dxbin<=0) dxbin=hx->GetBinWidth(1);
  }
  if (idxbin==-2)
  {
   dxbin=hx->GetBinWidth(1);
   dxbin=dxbin*float(nc);
  }

  // Set the auto-determined range of the output histogram
  if (dxmin<0)
  {
   dxmin=deltaxmin;
   // Compensate for randomized x-values within the input histogram bins
   if (mode==1 || mode==3)
   {
    bsize=hx->GetBinWidth(1);
    dxmin=dxmin-2.*bsize;
    if (dxmin<0) dxmin=0;
   }
  }
  if (dxmax<0)
  {
   dxmax=deltaxmax+dxbin;
   // Compensate for randomized x-values within the input histogram bins
   if (mode==1 || mode==3)
   {
    bsize=hx->GetBinWidth(1);
    dxmax=dxmax+2.*bsize;
   }
  }

  // Invoke the recursive call to create and fill the output histogram
  hdx=GetDxHistogram4(hx,nc,dxbin,dxmin,dxmax,mode,fact);
 }

 return hdx;
}
///////////////////////////////////////////////////////////////////////////
