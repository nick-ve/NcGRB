///////////////////////////////////////////////////////////////////////
// ROOT macro to convert Fermi t90 data into plain TTree format.     //
//                                                                   //
// The produced ROOT output file will have the same name as the      //
// .txt input file, but with the extension .txt replaced by .root    //
//                                                                   //
// To execute the user macro, just invoke ($ = prompt)               //
// $ root -b -q GRB-t90-Fermi.cc                                     //
//                                                                   //
//- Nick van Eijndhoven February 19, 2019  15:34 IIHE-VUB, Brussels. //
///////////////////////////////////////////////////////////////////////
{
 // The name of the T90 input txt file
 TString ifname="GRB-t90-Fermi.txt";

 // Construct the name of the produced ROOT output file
 TString ofname=ifname;
 ofname.ReplaceAll(".txt",".root");

 cout << endl;
 cout << " Input  filename : " << ifname << endl;
 cout << " Output filename : " << ofname << endl;

 // Input variables
 TString grbname;
 Float_t t90;

 // The input data file for observed GRB t90 durations
 ifstream ifile;
 ifile.clear();
 ifile.open(ifname.Data());
 if (!ifile.good())
 {
  cout << " *** Data file for observed GRB T90 durations not found ***" << endl;
  return;
 }
 ifile.seekg(0); // Position at begin of file

 // The produced output structure
 TFile* ofile=new TFile(ofname.Data(),"RECREATE","Fermi GRB T90");
 TTree* otree=new TTree("T","Fermi GRB T90 data");
 otree->Branch("t90",&t90,"t90/F");

 // Read the title lines in the input txt file
 string line;
 for (Int_t i=0; i<3; i++)
 {
  getline(ifile,line);
 }

 // Read the data and fill the output Tree
 Int_t n=0;
 while (ifile >> grbname >> t90)
 {
  otree->Fill();
  n++;
 }

 cout << endl;
 cout << " Number of data entries in the output Tree : " << n << endl;

 // Write the produced structure to the output file
 ofile->Write();
 ofile->Close();
}
