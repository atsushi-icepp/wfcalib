Double_t GeneralizedPoisson(Double_t *x, Double_t *par);
Int_t    Combination(Int_t n, Int_t r);
Double_t Borel(Double_t *x, Double_t *par);
Double_t GaussPH(Double_t *x,Double_t *par);
Double_t SigmaK(Double_t *x,Double_t *par);
Double_t DiffProb(Double_t *x,Double_t *par);
Double_t SingleAPProb(Double_t *x,Double_t *par);

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
	if (PH>PHexp) {
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
