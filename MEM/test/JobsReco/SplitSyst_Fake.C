
void CopyTree(){


}

void SplitSyst_Fake(string FileName){

  TFile* f = TFile::Open(FileName.c_str());

  //TTree* t_JESplus = (TTree*) f->Get("JES__plus");
  TTree* t_Fakeminus = (TTree*) f->Get("Fakes__minus");

  string str=".root";
  size_t found = FileName.find(str);

  string strinit = "FCNCNTuple";
  size_t foundinit = FileName.find(strinit);

  //string FileName_tmp = FileName;
  //FileName_tmp.replace(found, str.length(), "_JECplus.root");
  //string FileName_JESplus = FileName_tmp.substr(foundinit, FileName_tmp.length()-foundinit);
  //FileName_JESplus = "/tmp/chanon/" + FileName_JESplus;

  FileName_tmp = FileName;
  FileName_tmp.replace(found, str.length(), "_minus.root");
  string FileName_Fakeminus = FileName_tmp.substr(foundinit, FileName_tmp.length()-foundinit);
  FileName_Fakeminus = "/tmp/chanon/" + FileName_Fakeminus;

  //cout << FileName_JESplus.c_str() << endl;
  
  //TFile* f_JECplus = new TFile(FileName_JESplus.c_str(),"RECREATE");
  //f_JECplus->cd();
  //TTree* t_JESplus2 = t_JESplus->CopyTree("");
  //t_JESplus2->SetName("Tree");
  //f_JECplus->Write(); 

  cout << FileName_Fakeminus.c_str() << endl;

  TFile* f_Fakeminus = new TFile(FileName_Fakeminus.c_str(),"RECREATE");
  f_Fakeminus->cd();
  TTree* t_Fakeminus2 = t_Fakeminus->CopyTree("");
  t_Fakeminus2->SetName("Tree");
  f_Fakeminus->Write();


  return;
}
