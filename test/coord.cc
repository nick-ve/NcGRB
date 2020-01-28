///////////////////////////////////////////////////////////////////
// ROOT macro to test IceCube coordinate system.                 //
//                                                               //
//--- Nick van Eijndhoven 03-sep-2013 IIHE-VUB, Brussels         //
///////////////////////////////////////////////////////////////////
{
 gSystem->Load("ncfspack");

 // Load and compile the user code
 gROOT->LoadMacro("NcAstrolab2.cxx+");

 NcVersion version;
 version.Data();

 NcAstrolab2 lab;
 lab.SetExperiment("IceCube");

 Double_t offset=lab.GetLabTimeOffset();
 lab.Data();

 Double_t alpha,delta,az,zen;
 NcTimestamp tx;

 // IceCube Alert 
 tx.SetUT("20-01-2020","18:49:08",0);
 tx.Date(-1,offset);
 alpha=67.7381;
 delta=-14.5870;
 lab.SetSignal(1,alpha,"deg",delta,"deg","equ",&tx,-1,"J","Alert",1);
 // Add the Sun and the Moon
 lab.GetSignal("Sun");
 lab.GetSignal("Moon");
 lab.PrintSignal("equ","J",&tx,5,1,"T",1); cout << endl;
 lab.PrintSignal("equ","J",&tx,5,1,"T",0); cout << endl;
 lab.PrintSignal("equ","J",&tx,5,2,"T",0); cout << endl;

 // An IceCube event with RA=208.1197 Dec=48.1278 Az=264.6177 Zen=138.0735
 tx.SetMJD(55694.422336338430);
 tx.Date(-1,offset);
 alpha=208.1197;
 delta=48.1278;
 az=264.6177;
 zen=138.0735;
 lab.SetSignal(1,zen,"deg",az,"deg","loc",&tx,-1,"T","IC1",1);
 lab.PrintSignal("equ","J",&tx,5,2,"T",1); cout << endl;

 // IceCube GFU event run=129525 event=15195783
 tx.SetMJD(57891.571848993372);
 tx.Date(-1,offset);
 alpha=17.1372;
 delta=-72.6816;
 az=155.0528;
 zen=17.4066;
 lab.SetSignal(1,zen,"deg",az,"deg","loc",&tx,-1,"T","GFU",1);
 lab.PrintSignal("equ","J",&tx,5,3,"T",1); cout << endl;

 tx.SetMJD(55864.27427937708125682548);
 tx.Date(-1,offset);
 alpha=242.2363;
 delta=53.58844;
 az=344.6789;
 zen=143.5657;
 lab.SetSignal(1,zen,"deg",az,"deg","loc",&tx,-1,"T","IC2",1);
 lab.PrintSignal("equ","J",&tx,5,4,"T",1); cout << endl;

 // Enter the arrival direction of Kloppo as observed event
 tx.SetMJD(56819.20444852863); // The event time of Kloppo
 tx.Date(-1,offset);
 alpha=110.34;
 delta=11.48;
 zen=101.47741;
 az=312.72068;
 lab.SetSignal(1,zen,"deg",az,"deg","loc",&tx,-1,"M","Kloppo",1);
 lab.PrintSignal("equ","J",&tx,5,5,"T",1); cout << endl;

 // Print the signal info for the timestamp at the observation
 cout << endl;
 lab.ListSignals("loc","T",5);
 cout << endl;
 cout << " ----------------------" << endl;
 lab.ListSignals("equ","J",5);
 cout << endl;
 cout << " ----------------------" << endl;
 lab.ListSignals("equ","M",5);
 cout << endl;
 cout << " ----------------------" << endl;
 lab.ListSignals("equ","T",5);
 cout << endl;
 cout << " ----------------------" << endl;
 lab.ListSignals("hor","T",5);
 cout << endl;
 cout << " ----------------------" << endl;
 lab.ListSignals("equ","B",5);
 cout << endl;
 cout << " ----------------------" << endl;
}
