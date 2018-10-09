Double_t ExpectedSpectrum(Double_t *x,Double_t *par);
Double_t GeneralizedPoisson(Double_t *x, Double_t *par);
Int_t    Combination(Int_t n, Int_t r);
Double_t Borel(Double_t *x, Double_t *par);
Double_t GaussPH(Double_t *x,Double_t *par);
Double_t SigmaK(Double_t *x,Double_t *par);
Double_t DiffProb(Double_t *x,Double_t *par);
Double_t SingleAPProb(Double_t *x,Double_t *par);

Double_t ExpectedSpectrum(Double_t *x,Double_t *par){
	Double_t xx    = x[0];
	Double_t kmax  = par[0];
	Double_t alpha = par[1];
	Double_t mu    = par[2];
	Double_t lambda= par[3];
	Double_t sigma_0= par[4];
	Double_t sigma_1= par[5];
	Double_t beta   = par[6];
	Double_t scale  = par[7];
	Double_t gain   =1.;
	Double_t ped    =0.;
	TF1* fgpoi    = new TF1("fgpoi",GeneralizedPoisson,-100,100,2);
	TF1* fgph     =  new TF1("fgph",GaussPH,-100,100,4);
	TF1* fborel   = new TF1("fborel",Borel,-100,100,2);
	TF1* fdpdph_1 = new TF1("fdpdph_1",SingleAPProb,-100,100,5);
	TF1* fdpdph_i = new TF1("fdpdph_i",DiffProb,-100,100,5);
	TF1* fsigmak  = new TF1("fsigmak",SigmaK,0,100,2);
	fgpoi->SetParameters(mu,lambda);
	fgph-> SetParameters(0,sigma_0,ped,gain);
	Double_t vgpoi_ped=fgpoi-> Eval(0);
	Double_t vgph_ped =fgph -> Eval(xx);
	Double_t term_ped = vgpoi_ped*vgph_ped;
	Double_t term_ill = 0;
	fsigmak->SetParameters(sigma_0,sigma_1);
	Double_t maxGP = 0;
	for(Int_t k=1;k<kmax+1;k++){
		Double_t kk=(Double_t)k;
		Double_t sigma_k= fsigmak->Eval(kk);
		fgpoi  -> SetParameters(mu,lambda);
		Double_t GP= fgpoi->Eval(kk);
		if (GP>maxGP) maxGP = GP;
		if (GP<maxGP/TMath::Power(10,5)) continue;
		// i == 0
		fborel -> SetParameters(0,alpha);
		fgph   -> SetParameters(kk,sigma_k,ped,gain);
		Double_t vborel_0 = fborel->Eval(kk);
		Double_t vgrph_0  = fgph  ->Eval(xx);
		Double_t term_0 = vborel_0*vgrph_0;
		// i == 1
		fborel -> SetParameters(1,alpha);
		fdpdph_1 -> SetParameters(kk,beta,sigma_k,ped,gain);
		Double_t vborel_1 = fborel  ->Eval(kk);
		Double_t vdpdph_1 = fdpdph_1->Eval(xx);
		Double_t term_1=vborel_1* vdpdph_1;
		// i > 2
		Double_t term_i=0;
		for(Int_t i=2;i<k+1;i++){
			fborel   -> SetParameters(i,alpha);
			fdpdph_i -> SetParameters(kk,i,beta,ped,gain);
			term_i+=fborel->Eval(kk)*fdpdph_i->Eval(xx);
		}
		term_ill+=GP*(term_0+term_1+term_i);
		// term_ill+=GP*term_0;
		// std::cout<<"term_0"<<term_0<<std::endl;
		// std::cout<<"term_1"<<term_1<<std::endl;
		// std::cout<<"term_i"<<term_i<<std::endl;
		// std::cout<<"GP    "<<GP<<std::endl;
		// std::cout<<"k: "<<k<<" GP: "<<GP<<"term 0: "<<term_0<<" term 1: "<<term_1<<" term_i: "<<term_i<<std::endl;
	}
	// std::cout<<"term_ill"<<term_ill<<std::endl;
	return scale*(term_ped+ term_ill);
}

Double_t GeneralizedPoisson(Double_t *x, Double_t *par){
	Double_t k      = x[0];
	Double_t mu     = par[0];
	Double_t lambda = par[1];
	Double_t dom=mu*TMath::Power(mu+k*lambda,k-1);
	Double_t num=(Double_t)TMath::Factorial((Int_t)k);
	Double_t coeff=TMath::Exp(-mu-k*lambda);
	Double_t rtnval=dom*coeff/num;
	return rtnval;
}

// _nC_r 
Int_t Combination(Int_t n, Int_t r){
	int num = 1;
	for(int i = 1; i <= r; i++){
		num = num * (n - i + 1) / i;
	}
	return num;
}

// _kC_i*alph^i*(1-alpha)^{k-i}
Double_t Borel(Double_t *x, Double_t *par){
	Double_t k    = x[0];
	Double_t i    = par[0];
	Double_t alpha= par[1];
	return (Double_t)Combination((Int_t)k,(Int_t)i)*TMath::Power(alpha,i)*TMath::Power(1-alpha,k-i);
}

Double_t GaussPH(Double_t *x,Double_t *par){
	Double_t PH   = x[0];
	Double_t k    = par[0];
	Double_t sigma= par[1];
	Double_t ped  = par[2];
	Double_t gain = par[3];
	Double_t coeff= 1./TMath::Sqrt(2*TMath::Pi())/sigma;
	Double_t expo = TMath::Exp(-TMath::Power((PH-(ped+k*gain)),2)/(2.*sigma*sigma));
	return coeff*expo;
}

// sigma^2+k*sigma_1^2
Double_t SigmaK(Double_t *x,Double_t *par){
	Double_t k=x[0];
	Double_t sigma0= par[0];
	Double_t sigma1= par[1];
	return TMath::Sqrt(sigma0*sigma0+k*sigma1*sigma1);
}

Double_t DiffProb(Double_t *x,Double_t *par){
	Double_t PH  = x[0];
	Double_t k   = par[0];
	Double_t i   = par[1];
	Double_t beta= par[2];
	Double_t ped = par[3];
	Double_t gain= par[4];
	Double_t PHexp=ped+k*gain;
	if (PH>=PHexp) {
		Double_t dom = TMath::Power(PH-PHexp,i-1);
		Double_t num = TMath::Factorial((Int_t)(i-1))*TMath::Power(beta,i);
		Double_t coeff = TMath::Exp(-(PH-PHexp)/beta);
		return dom*coeff/num;
	}else{
		return 0;
	}
}

Double_t SingleAPProb(Double_t *x,Double_t *par){
	Double_t PH  = x[0];
	Double_t k      = par[0];
	Double_t beta   = par[1];
	Double_t sigma_k= par[2];
	Double_t ped    = par[3];
	Double_t gain   = par[4];
	Double_t PHexp=ped+k*gain;
	Double_t dom = TMath::Exp(-(PH-PHexp)/beta);
	Double_t coeff = TMath::Sqrt(2*TMath::Pi())*sigma_k*beta;
	Double_t inte = TMath::Sqrt(TMath::Pi()/2.)*sigma_k*TMath::Erfc(-(PH-PHexp)/TMath::Sqrt(2)/sigma_k);
	return dom*inte/coeff;
}
