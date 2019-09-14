///////////////////////////////////////////////////////////////////
// ROOT macro to process IceCube data in plain TTree format.     //
// The actual C++ code to be processed is contained in a file    //
// called grbmodel.cxx with "void user()" as main function.      //
//                                                               //
// To execute the user code, just invoke ($ = prompt)            //
// $ root -b -q grbmodel-full.cc                                 //
//                                                               //
//--- Nick van Eijndhoven 03-sep-2013 IIHE-VUB, Brussels         //
///////////////////////////////////////////////////////////////////
{
 gSystem->Load("ncfspack");
 gSystem->Load("icepack");

 // Load and compile the user code
 gROOT->LoadMacro("NcGRB-full.cxx+");

 // Execute the user code
 user();
}
