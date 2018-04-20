
void CopyTree(){


}

void SystMerger(string FileName){

  TFile* f = TFile::Open(FileName.c_str());

  string str=".root";
  size_t found = FileName.find(str);
  string strinit = "FCNCNTuple";
  size_t foundinit = FileName.find(strinit);

  string FileName_syst[7];
  string FileName_VarSyst[7];
  TFile* f_VarSyst[7];
  TTree* t_VarSyst[7];

  string extension[7];
  string treeName[7];
  TTree* t[7];
  extension[0] = ".root";
  treeName[0] = "Tree";
  int next = 5;

  FileName_syst[0] = FileName;
  FileName_syst[0].replace(found, str.length(), ".root");
  FileName_VarSyst[0] = FileName_syst[0].substr(foundinit, FileName_syst[0].length()-foundinit);
  string cmpstring = "FCNCNTuple_tZq.root";
  if (FileName_VarSyst[0].compare(cmpstring)==0) next = 7;
  cmpstring = "FCNCNTuple_FakesNew2.root";
  if (FileName_VarSyst[0].compare(cmpstring)==0) next = 2;
  FileName_VarSyst[0] = "/tmp/chanon/" + FileName_VarSyst[0];
  cout << FileName_VarSyst[0].c_str() << " next="<< next<<endl;
  f_VarSyst[0] = new TFile(FileName_VarSyst[0].c_str(), "READ");
  t_VarSyst[0] = (TTree*) f_VarSyst[0]->Get("Tree");
  t_VarSyst[0]->SetBranchStatus("*",1);
  t_VarSyst[0]->SetBranchStatus("multilepton_h0_P4",0);
  t_VarSyst[0]->SetBranchStatus("multilepton_t1_P4",0);
  t_VarSyst[0]->SetBranchStatus("multilepton_t2_P4",0);

  if (next==2){
    extension[1] = "_minus.root";
    treeName[1] = "Fakes__minus";
  }
  if (next>=5) {
    extension[1] = "_JECplus.root";
    extension[2] = "_JECminus.root";
    extension[3] = "_JERplus.root";
    extension[4] = "_JERminus.root";

    treeName[1] = "JES__plus";
    treeName[2] = "JES__minus";
    treeName[3] = "JER__plus";
    treeName[4] = "JER__minus";
  }
  if (next==7){
    extension[5] = "Qdw.root";
    extension[6] = "Qup.root";

    treeName[5] = "Qdw";
    treeName[6] = "Qup";
  }

  for (int i=1; i<next; i++){
    FileName_syst[i] = FileName;
    FileName_syst[i].replace(found, str.length(), extension[i].c_str());
    FileName_VarSyst[i] = FileName_syst[i].substr(foundinit, FileName_syst[i].length()-foundinit);
    FileName_VarSyst[i] = "/tmp/chanon/" + FileName_VarSyst[i];
    f_VarSyst[i] = new TFile(FileName_VarSyst[i].c_str(), "READ");
    t_VarSyst[i] = (TTree*) f_VarSyst[i]->Get("Tree");
    t_VarSyst[i]->SetBranchStatus("*",1);
    t_VarSyst[i]->SetBranchStatus("multilepton_h0_P4",0);
    t_VarSyst[i]->SetBranchStatus("multilepton_t1_P4",0);
    t_VarSyst[i]->SetBranchStatus("multilepton_t2_P4",0);
  }

  string FileNameOutput = FileName;
  FileNameOutput.replace(found, str.length(), "_withSyst.root");
  string FileNameOutputWithDir = FileNameOutput.substr(foundinit, FileNameOutput.length()-foundinit);
  FileNameOutputWithDir = "/tmp/chanon/" + FileNameOutputWithDir;
  TFile* fOutput = new TFile(FileNameOutputWithDir.c_str(), "RECREATE");
  fOutput->cd();
  for (int i=0; i<next; i++){
    t[i] = t_VarSyst[i]->CopyTree("");
    t[i]->SetName(treeName[i].c_str());
  }
  fOutput->Write();


  return;
/*
  FileName_syst[1] = FileName;
  FileName_syst[1].replace(found, str.length(), "_JECplus.root");
  FileName_VarSyst[1] = FileName_syst[1].substr(foundinit, FileName_syst[1].length()-foundinit);
  FileName_VarSyst[1] = "/tmp/chanon/" + FileName_VarSyst[1];
  f_VarSyst[1] = new TFile(FileName_VarSyst[1].c_str(), "READ");
  t_VarSyst[1] = (TTree*) f_VarSyst[1]->Get("Tree");
  t_VarSyst[1]->SetBranchStatus("*",1);
  t_VarSyst[1]->SetBranchStatus("multilepton_h0_P4",0);
  t_VarSyst[1]->SetBranchStatus("multilepton_t1_P4",0);
  t_VarSyst[1]->SetBranchStatus("multilepton_t2_P4",0);

  FileName_syst[2] = FileName;
  FileName_syst[2].replace(found, str.length(), "_JECminus.root");
  FileName_VarSyst[2] = FileName_syst[2].substr(foundinit, FileName_syst[2].length()-foundinit);
  FileName_VarSyst[2] = "/tmp/chanon/" + FileName_VarSyst[2];
  f_VarSyst[2] = new TFile(FileName_VarSyst[2].c_str(), "READ");
  t_VarSyst[2] = (TTree*) f_VarSyst[2]->Get("Tree");
  t_VarSyst[2]->SetBranchStatus("*",1);
  t_VarSyst[2]->SetBranchStatus("multilepton_h0_P4",0);
  t_VarSyst[2]->SetBranchStatus("multilepton_t1_P4",0);
  t_VarSyst[2]->SetBranchStatus("multilepton_t2_P4",0);

  FileName_syst[3] = FileName;
  FileName_syst[3].replace(found, str.length(), "_JERplus.root");
  FileName_VarSyst[3] = FileName_syst[3].substr(foundinit, FileName_syst[3].length()-foundinit);
  FileName_VarSyst[3] = "/tmp/chanon/" + FileName_VarSyst[3];
  f_VarSyst[3] = new TFile(FileName_VarSyst[3].c_str(), "READ");
  t_VarSyst[3] = (TTree*) f_VarSyst[3]->Get("Tree");
  t_VarSyst[3]->SetBranchStatus("*",1);
  t_VarSyst[3]->SetBranchStatus("multilepton_h0_P4",0);
  t_VarSyst[3]->SetBranchStatus("multilepton_t1_P4",0);
  t_VarSyst[3]->SetBranchStatus("multilepton_t2_P4",0);

  FileName_syst[4] = FileName;
  FileName_syst[4].replace(found, str.length(), "_JERminus.root");
  FileName_VarSyst[4] = FileName_syst[4].substr(foundinit, FileName_syst[4].length()-foundinit);
  FileName_VarSyst[4] = "/tmp/chanon/" + FileName_VarSyst[4];
  f_VarSyst[4] = new TFile(FileName_VarSyst[4].c_str(), "READ");
  t_VarSyst[4] = (TTree*) f_VarSyst[4]->Get("Tree");
  t_VarSyst[4]->SetBranchStatus("*",1);
  t_VarSyst[4]->SetBranchStatus("multilepton_h0_P4",0);
  t_VarSyst[4]->SetBranchStatus("multilepton_t1_P4",0);
  t_VarSyst[4]->SetBranchStatus("multilepton_t2_P4",0);


  string FileNameOutput = FileName;
  FileNameOutput.replace(found, str.length(), "_withSyst.root");
  string FileNameOutputWithDir = FileNameOutput.substr(foundinit, FileNameOutput.length()-foundinit);
  FileNameOutputWithDir = "/tmp/chanon/" + FileNameOutputWithDir;
  TFile* fOutput = new TFile(FileNameOutputWithDir.c_str(), "RECREATE");
  fOutput->cd();
  TTree* t[5];
  t[0] = t_VarSyst[0]->CopyTree("");
  t[1] = t_VarSyst[1]->CopyTree("");
  t[1]->SetName("JES__plus");
  t[2] = t_VarSyst[2]->CopyTree("");
  t[2]->SetName("JES__minus");
  t[3] = t_VarSyst[3]->CopyTree("");
  t[3]->SetName("JER__plus");
  t[4] = t_VarSyst[4]->CopyTree("");
  t[4]->SetName("JER__minus");
  fOutput->Write();
*/
  return;

/*
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
*/

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
