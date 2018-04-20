
void CopyTree(){


}

void SplitSyst(string FileName){

  TFile* f = TFile::Open(FileName.c_str());

  string str=".root";
  size_t found = FileName.find(str);
  string strinit = "FCNCNTuple";
  size_t foundinit = FileName.find(strinit);
  string FileName_tmp;


  TTree* t_JESplus = (TTree*) f->Get("JES__plus");
  FileName_tmp = FileName;
  FileName_tmp.replace(found, str.length(), "_JECplus.root");
  string FileName_JESplus = FileName_tmp.substr(foundinit, FileName_tmp.length()-foundinit);
  FileName_JESplus = "/tmp/chanon/" + FileName_JESplus;
  cout << FileName_JESplus.c_str() << endl;
  TFile* f_JECplus = new TFile(FileName_JESplus.c_str(),"RECREATE");
  f_JECplus->cd();
  TTree* t_JESplus2 = t_JESplus->CopyTree("");
  t_JESplus2->SetName("Tree");
  f_JECplus->Write(); 

  TTree* t_JESminus = (TTree*) f->Get("JES__minus");
  FileName_tmp = FileName;
  FileName_tmp.replace(found, str.length(), "_JECminus.root");
  string FileName_JESminus = FileName_tmp.substr(foundinit, FileName_tmp.length()-foundinit);
  FileName_JESminus = "/tmp/chanon/" + FileName_JESminus;
  cout << FileName_JESminus.c_str() << endl;
  TFile* f_JECminus = new TFile(FileName_JESminus.c_str(),"RECREATE");
  f_JECminus->cd();
  TTree* t_JESminus2 = t_JESminus->CopyTree("");
  t_JESminus2->SetName("Tree");
  f_JECminus->Write();


  TTree* t_JERplus = (TTree*) f->Get("JER__plus");
  FileName_tmp = FileName;
  FileName_tmp.replace(found, str.length(), "_JERplus.root");
  string FileName_JERplus = FileName_tmp.substr(foundinit, FileName_tmp.length()-foundinit);
  FileName_JERplus = "/tmp/chanon/" + FileName_JERplus;
  cout << FileName_JERplus.c_str() << endl;
  TFile* f_JERplus = new TFile(FileName_JERplus.c_str(),"RECREATE");
  f_JERplus->cd();
  TTree* t_JERplus2 = t_JERplus->CopyTree("");
  t_JERplus2->SetName("Tree");
  f_JERplus->Write(); 

  TTree* t_JERminus = (TTree*) f->Get("JER__minus");
  FileName_tmp = FileName;
  FileName_tmp.replace(found, str.length(), "_JERminus.root");
  string FileName_JERminus = FileName_tmp.substr(foundinit, FileName_tmp.length()-foundinit);
  FileName_JERminus = "/tmp/chanon/" + FileName_JERminus;
  cout << FileName_JERminus.c_str() << endl;
  TFile* f_JERminus = new TFile(FileName_JERminus.c_str(),"RECREATE");
  f_JERminus->cd();
  TTree* t_JERminus2 = t_JERminus->CopyTree("");
  t_JERminus2->SetName("Tree");
  f_JERminus->Write();

  return;
}
