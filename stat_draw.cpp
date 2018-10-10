Int_t MAX(Int_t a,Int_t b){if(a>b){return a;}else{return b;}}
Double_t fitfunc(Double_t *x, Double_t *par);
Double_t DivisionError(Double_t &value1, Double_t &value1err,const Double_t &value2, Double_t &value2err);
Double_t MultipleError(Double_t &value1, Double_t &value1err, Double_t &value2, Double_t &value2err);

void stat_draw(void) {
   Double_t min = 0.;
   Double_t max = 1000.;
   TCanvas *cgr0 =new TCanvas("cgr0", "cgr0",600,600);
   TCanvas *cgr1 =new TCanvas("cgr1", "cgr1",600,600);
//   TCanvas *cgr1 =new TCanvas("cgr1", "cgr1",600,600);
   Int_t maxlist = 5; // this value must be at least the max size of noiselist
   TClonesArray* cagrQNvar = new TClonesArray("TGraphErrors",maxlist);
   TClonesArray* cafpol1    = new TClonesArray("TF1",maxlist);
   for (int i = 0; i < maxlist; i++) {
      new ((TGraphErrors*)(*cagrQNvar)[i]) TGraphErrors();
      new ((TF1*)(*cafpol1)[i])   TF1(Form("fpol1%d",i),"[0]*x+[1]");
   }
   
   TGraph *meanvar = new TGraph();
   TFile* fin = new TFile("fout.root");
   TTree* tin = (TTree*)fin->Get("tout");
   Int_t noiseNum,truemaxNum;
   Double_t NphoMean;
   Double_t *charge = new Double_t[maxlist];
   Double_t *noisevar = new Double_t[maxlist];
   Double_t *noiselevel = new Double_t[maxlist];
   Int_t *Npho = new Int_t[maxlist];
   tin->SetBranchAddress("charge", charge);
   tin->SetBranchAddress("noisevar", noisevar);
   tin->SetBranchAddress("noiselevel", noiselevel);
//   tin->SetBranchAddress("noiselist", &noiseNum);
   tin->SetBranchAddress("NphoMean", &NphoMean);
   tin->SetBranchAddress("Npho",Npho);

   tin->GetEntry(0);
   Double_t fix_Npho = NphoMean;
   Int_t index = 0;
   TH1F *hist = new TH1F("hist","mean vs var",300, min, max);
   for(int ient = 0; ient < tin->GetEntries(); ++ient) {
      tin->GetEntry(ient);
      Double_t comp_Npho = NphoMean;
      if (comp_Npho!=fix_Npho){
         fix_Npho = comp_Npho;
         Double_t mean = hist->GetMean();
         Double_t var = hist->GetStdDev();
         var = var*var;
         meanvar->SetPoint(index,mean,var);
         printf("index is %d mean is %lf var is %lf,\n",index,mean,var);
         index++;
         delete hist;
         TH1F *hist = new TH1F("hist","mean vs var",300, min, max);
      }
      hist->Fill(charge[0]);
   }
   Double_t mean = hist->GetMean();
   Double_t var = hist->GetStdDev();
   var = var*var;
   meanvar->SetPoint(index,mean,var);
   meanvar->SetMarkerColor(kRed);
   meanvar->SetMarkerSize(2.);
   meanvar->SetMarkerStyle(22);
    
   cgr0->cd();
   meanvar->Draw("ap");
   cgr1->cd();
   hist->Draw();
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

