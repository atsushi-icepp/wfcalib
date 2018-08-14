#include "phFunc.h"

void phrec(){
	TFile* fin = new TFile("fout.root","read");
	TTree* tin =(TTree*)fin->Get("tout");
	Double_t charge;
	Double_t noiselevel;
	tin->SetBranchAddress("charge",&charge);
	tin->SetBranchAddress("noiselevel",&noiselevel);
	TH1D* hQ=new TH1D("hQ","hQ",100,-0.5,5);
	Int_t Nentry=tin->GetEntries();
	for (Int_t ientry = 0; ientry < Nentry; ientry++) {
		tin->GetEntry(ientry);
		if (noiselevel==0.005) {
			hQ->Fill(charge);
		}
	}
	hQ->Draw();
	Double_t kmax  = 8;
	Double_t alpha = 0.1;
	Double_t mu    = 1;
	Double_t lambda= 0.06;
	Double_t sigma_0=0.08;
	Double_t sigma_1=0.0;
	Double_t beta   =100;
	Double_t scale  =20;
	TF1* fexp= new TF1("fexp",ExpectedSpectrum,-0.5,20,8);
	fexp->SetParameters(kmax,alpha,mu,lambda,sigma_0,sigma_1,beta,scale);
	fexp->SetParNames(
		"kmax",
		"alpha",
		"mu",
		"lamdba",
		"sigma0",
		"sigma1",
		"beta",
		"scale"
	);
	fexp->FixParameter(0,kmax);
	// fexp->FixParameter(2,mu);
	fexp->SetParLimits(1,0,1);
	// fexp->FixParameter(4,0.08);
	fexp->SetParLimits(5,0,20);
	// std::cout<<"debug"<<std::endl;
	// TF1* fborel= new TF1("fborel",Borel);
	// fborel->Draw();
	// fexp->Eval(0);
	// std::cout<<fexp->Eval(0)<<std::endl;
	// fexp->Draw("pl");
	// fexp->SetMarkerStyle(3);
	hQ->Fit("fexp","","",-5,5);
}

