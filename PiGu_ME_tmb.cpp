// 
#include <TMB.hpp>


// function to multiply vector A by matrix B and put into vector C; O. Giminez
template<class Type>

vector<Type> multvecmat(array<Type>  A, matrix<Type>  B) {
  int nrowb = B.rows();
  int ncolb = B.cols(); 
  
  vector<Type> C(ncolb);
  for (int i = 0; i < ncolb; i++)
  {
	    C(i) = Type(0);
      for (int k = 0; k < nrowb; k++){
        C(i) += A(k)*B(k,i);
    }
  }
  return C;
}

template<class Type>
Type objective_function<Type>::operator() () {
  
  PARAMETER_VECTOR(params);
  
  DATA_IMATRIX(CH_MEJS);
  DATA_IVECTOR(first_capt);
  DATA_IVECTOR(init_state);
  
  // OBSERVATIONS
  // 1 = burrow visit
  // 2 = prey visit
  // 3 = not detected
  //   
  // STATES
  // 1 = not entered
  // 2 = egg
  // 3 = nestling
  // 4 = terminated
  //   
  // PARAMETERS
  // phiA: survival probability from egg to chick
  // phiB: survival probability from chick to fledge
  // psiAB: probability of transitioning from egg to chick
  // pA: detection probability of egg burrow
  // pB: detection probability of chick burrow
  // b: conditional on observing a visit to a chick burrow, probability the delivery was a 'burrow visit' 
  // gamma: entry probability
  // alpha: conditional on entry at occasion 1, probability burrow had a chick. 
  //    alpha set to 0 for t>1. So if a burrow is initiated after day one, it must start as an egg burrow.
  
  // logit link for all parameters
  
  int Nind = CH_MEJS.rows();
  int Nocc = CH_MEJS.cols(); 
  int Nobs = 3;
  int Nstate = 4;
  int Npar = params.size();
  vector<Type> par(Npar);
  for (int i = 0; i < Npar; i++) {
    par(i) = Type(1.0) / (Type(1.0) + exp(-params(i)));
  }
  Type phiA = par(0); 
  Type phiB = par(1);
  Type psiAB = par(2);
  Type pA = par(3);
  Type pB = par(4);
  Type b = par(5)
  Type gamma = par(6);
  Type alpha = par(7);
  
  // Create the state and obs matrices in two steps composed of binomial events to use logit link, 
  // if not, must use multinomial logit link 
  
  // prob of obs (cols) cond on states (rows)
  matrix<Type> ObsMat(Nstate, Nobs);
    ObsMat(0,0) = 0;
    ObsMat(0,1) = 0;
    ObsMat(0,2) = 1;

    ObsMat(1,0) = pA[t]*eff[i,t];
    ObsMat(1,1) = 0;
    ObsMat(1,2) = 1-pA[t]*eff[i,t];

    ObsMat(2,0) = b * pB[t]*eff[i,t];
    ObsMat(2,1) = (1 - b) * pB[t]*eff[i,t];
    ObsMat(2,2) = 1 - pB[t]*eff[i,t];

    ObsMat(3,0) = 0;
    ObsMat(3,1) = 0;
    ObsMat(3,2) = 1;

  REPORT(ObsMat);
  
  // first encounter
  // matrix<Type> BE1(3, 3);
  // BE1(0,0) = Type(0.0);
  // BE1(0,1) = Type(1.0);
  // BE1(0,2) = Type(0.0);
  // BE1(1,0) = Type(0.0);
  // BE1(1,1) = Type(0.0);
  // BE1(1,2) = Type(1.0);
  // BE1(2,0) = Type(1.0);
  // BE1(2,1) = Type(0.0);
  // BE1(2,2) = Type(0.0);
  // 
  // matrix<Type> BE2(3, 4);
  // BE2(0,0) = Type(1.0);
  // BE2(0,1) = Type(0.0);
  // BE2(0,2) = Type(0.0);
  // BE2(0,3) = Type(0.0);
  // BE2(1,0) = Type(0.0);
  // BE2(1,1) = deltaNB;
  // BE2(1,2) = Type(0.0);
  // BE2(1,3) = Type(1.0)-deltaNB;
  // BE2(2,0) = Type(0.0);
  // BE2(2,1) = Type(0.0);
  // BE2(2,2) = deltaB;
  // BE2(2,3) = Type(1.0)-deltaB;
  // 
  // matrix<Type> BEE(3,4);
  // BEE = BE1 * BE2;
  // matrix<Type> BE(4,3);
  // BE = BEE.transpose();
  // REPORT(BE);
  
  // state matrix: prob of state at t+1 given state at t
  matrix<Type> StMat(4, 4);
  
  StMat(0,0) = 1-gamma[t];
  StMat(0,1) = gamma[t]*(1-alpha[t]);
  StMat(0,2) = gamma[t]*alpha[t];
  StMat(0,3) = 0;
  
  StMat(1,0) = 0;
  StMat(1,1) = phiA[t]*(1-psiAB[t]);
  StMat(1,2) = phiA[t]*psiAB[t];
  StMat(1,3) = 1-phiA[t];
  
  StMat(2,0) = 0;
  StMat(2,1) = 0;
  StMat(2,2) = phiB[t];
  StMat(2,3) = 1-phiB[t];
  
  StMat(3,0) = 0;
  StMat(3,1) = 0;
  StMat(3,2) = 0;
  StMat(3,3) = 1;
  
  // initial states
  vector<Type> StateInit(3);
  StateInit(0) =   //not entered
  StateInit(1) = alpha;  //egg
  StateInit(2) = Type(1.0)-alpha; //chick
  StateInit(3) = Type(0.0); //term
  REPORT(StateInit);
  
  // likelihood
  Type likelihood = 0;
  Type neglogL;
  array<Type> temp(3);
  for (int t = 0; t < Nocc; t++) {
    int ei = first_capt(i)-1;
    vector<int> evennt = CH_MEJS.col(t);
    temp = StateInit * vector<Type>(BE.row(first_state(t))); // element-wise vector product
    for (int i = ei+1; j < Nobs; i++) {
      temp = multvecmat(temp,StMat) * vector<Type>(ObsMat.row(evennt(i))); // vector matrix product, then element-wise vector product
    }
    likelihood += log(sum(temp));
  }
  neglogL = -likelihood;
  return neglogL;
}

