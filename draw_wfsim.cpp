Int_t MAX(Int_t a,Int_t b){if(a>b){return a;}else{return b;}}
Double_t fitfunc(Double_t *x, Double_t *par);
Double_t DivisionError(Double_t &value1, Double_t &value1err,const Double_t &value2, Double_t &value2err);
Double_t MultipleError(Double_t &value1, Double_t &value1err, Double_t &value2, Double_t &value2err);

void draw_wfsim(void) {
   TCanvas *cgr0 =new TCanvas("cgr0", "cgr0",600,600);
   TCanvas *cgr1 =new TCanvas("cgr1", "cgr1",600,600);
   Int_t maxlist = 5; // this value must be at least the max size of noiselist
   TClonesArray* cagrQNvar = new TClonesArray("TGraphErrors",maxlist);
   TClonesArray* cafpol1    = new TClonesArray("TF1",maxlist);
   for (int i = 0; i < maxlist; i++) {
      new ((TGraphErrors*)(*cagrQNvar)[i]) TGraphErrors();
      new ((TF1*)(*cafpol1)[i])   TF1(Form("fpol1%d",i),"[0]*x+[1]");
   }

   TFile* fin = new TFile("fout.root");
   TTree* tin = (TTree*)fin->Get("tout");
   Int_t noiseNum,truemaxNum;
   Double_t *charge = new Double_t[maxlist];
   Double_t *noisevar = new Double_t[maxlist];
   Double_t *noiselevel = new Double_t[maxlist];
   tin->SetBranchAddress("charge", charge);
   tin->SetBranchAddress("noisevar", noisevar);
   tin->SetBranchAddress("noiselevel", noiselevel);
   tin->SetBranchAddress("noiselist", &noiseNum);

   for(int ient = 0; ient < tin->GetEntries(); ++ient) {
      tin->GetEntry(ient);
      truemaxNum = MAX(noiseNum,truemaxNum);
      for (int i = 0; i < noiseNum; i++){
         ((TGraphErrors*)(*cagrQNvar)[i])->SetPoint(ient,charge[i],noisevar[i]);
         ((TGraphErrors*)(*cagrQNvar)[i])->SetPointError(ient,0,noisevar[i]*0.2);
      }
   }

   TF1* func= new TF1("fitfunc",fitfunc,-10,5000,3);

   cgr0->cd();
   for (int i = 0; i < truemaxNum; i++) {
      Double_t par[3];
      Double_t err[3];
      TF1* pol1= new TF1("fitfunc",fitfunc,-100,5000,3);
      ((TGraphErrors*)(*cagrQNvar)[i])->Fit("fitfunc","","");
      for (int i = 0; i < 3; i++) {
         par[i]=pol1->GetParameter(i);
         err[i]=pol1->GetParError(i);
      }
      Double_t MesGain=par[1]/(1-par[1]*par[2]);
      Double_t pmul=par[1]*par[2];
      Double_t pmulerr=MultipleError(par[1],err[1],par[2],err[2]);
      Double_t MesGainerr=DivisionError(par[1],err[1],(double)1-pmul,pmulerr);
      std::cout<<"Gain: "<<MesGain<<"+-"<<MesGainerr<<std::endl;
      std::cout<<"Pmul: "<<pmul<<"+-"<<pmulerr<<std::endl;
      ((TGraphErrors*)(*cagrQNvar)[i])->SetTitle("Q_{dint}vs Q^{2}_{drms};Q_{dint};Q^{2}_{drms}");
      ((TGraph*)(*cagrQNvar)[i])->SetMaximum(500);
      ((TGraphErrors*)(*cagrQNvar)[i])->SetMinimum(0);
      ((TGraphErrors*)(*cagrQNvar)[i])->SetMarkerStyle(20);
      ((TGraphErrors*)(*cagrQNvar)[i])->SetMarkerColor(2+i);
      if (i==0) {
         ((TGraphErrors*)(*cagrQNvar)[i])->Draw("ap");
      }else{
         ((TGraphErrors*)(*cagrQNvar)[i])->Draw("p same");
      }
      delete pol1;
   }

}

Double_t fitfunc(Double_t *x, Double_t *par){
   //par[0]: offset;
   //par[1]:
   Double_t  xx        = x[0];
   // Double_t PMultiple = TMath::Sqrt(TMath::Power(par[2]*xx,2)/1+TMath::Power(par[2]*xx,2));
   // Double_t PMultiple = par[2];
   // Double_t Gain = 1;
   return par[0]+par[1]*(xx+par[2]*xx*xx);
}

Double_t DivisionError(Double_t &value1, Double_t &value1err,const Double_t &value2, Double_t &value2err){
   /*calculate error of value 1/ value 2*/
   return TMath::Sqrt(TMath::Power(value1err/value2,2)+TMath::Power(value1*value2err/(value2*value2),2));
}

Double_t MultipleError(Double_t &value1, Double_t &value1err, Double_t &value2, Double_t &value2err){
   /*calculate error of value 1* value 2*/
   return TMath::Sqrt(TMath::Power(value1err* value2,2)+TMath::Power(value1*value2err,2));
}

