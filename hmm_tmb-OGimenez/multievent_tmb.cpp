// linear regression
#include <TMB.hpp>

//template<class Type>
//matrix<Type> multmat(array<Type> A, matrix<Type> B) {
//  int nrowa = A.rows();
//  int ncola = A.cols(); 
//  int ncolb = B.cols(); 
//  matrix<Type> C(nrowa,ncolb);
//  for (int i = 0; i < nrowa; i++)
//  {
//    for (int j = 0; j < ncolb; j++)
//    {
//      C(i,j) = Type(0);
//      for (int k = 0; k < ncola; k++)
//        C(i,j) += A(i,k)*B(k,j);
//    }
//  }
//  return C;
//}
//
/* implement the vector - matrix product */
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
  
  // b = parameters
  PARAMETER_VECTOR(b);
  
  // ch = capture-recapture histories (individual format)
  // fc = date of first capture
  // fs = state at first capture
  DATA_IMATRIX(ch);
  DATA_IVECTOR(fc);
  DATA_IVECTOR(fs);
  
  // OBSERVATIONS
  // 0 = non-detected
  // 1 = seen and ascertained as non-breeder
  // 2 = seen and ascertained as breeder
  // 3 = not ascertained
  //   
  // STATES
  // 1 = alive non-breeder
  // 2 = alive breeder
  // 3 = dead
  //   
  // PARAMETERS
  // phiNB  survival prob. of non-breeders
  // phiB  survival prob. of breeders
  // pNB  detection prob. of non-breeders
  // pB  detection prob. of breeders
  // psiNBB transition prob. from non-breeder to breeder
  // psiBNB transition prob. from breeder to non-breeder
  // piNB prob. of being in initial state non-breeder
  // deltaNB prob to ascertain the breeding status of an individual encountered as non-breeder
  // deltaB prob to ascertain the breeding status of an individual encountered as breeder
  //   
  // logit link for all parameters
  // note: below, we decompose the state and obs process in two steps composed of binomial events, 
  // which makes the use of the logit link appealing; 
  // if not, a multinomial (aka generalised) logit link should be used

  int km = ch.rows();
  int nh = ch.cols();  
  int npar = b.size();
  vector<Type> par(npar);
  for (int i = 0; i < npar; i++) {
    par(i) = Type(1.0) / (Type(1.0) + exp(-b(i)));
  }
  Type piNB = par(0); // careful, indexing starts at 0 in Rcpp!
  Type phiNB = par(1);
  Type phiB = par(2);
  Type psiNBB = par(3);
  Type psiBNB = par(4);
  Type pNB = par(5);
  Type pB = par(6);
  Type deltaNB = par(7);
  Type deltaB = par(8);
  
  // prob of obs (rows) cond on states (col)
  matrix<Type> B1(3,3);
  B1(0,0) = Type(1.0)-pNB;
  B1(0,1) = pNB;
  B1(0,2) = Type(0.0);
  B1(1,0) = Type(1.0)-pB;
  B1(1,1) = Type(0.0);
  B1(1,2) = pB;
  B1(2,0) = Type(1.0);
  B1(2,1) = Type(0.0);
  B1(2,2) = Type(0.0);
  
  matrix<Type> B2(3, 4);
  B2(0,0) = Type(1.0);
  B2(0,1) = Type(0.0);
  B2(0,2) = Type(0.0);
  B2(0,3) = Type(0.0);
  B2(1,0) = Type(0.0);
  B2(1,1) = deltaNB;
  B2(1,2) = Type(0.0);
  B2(1,3) = 1-deltaNB;
  B2(2,0) = Type(0.0);
  B2(2,1) = Type(0.0);
  B2(2,2) = deltaB;
  B2(2,3) = Type(1.0)-deltaB;
  
  matrix<Type> BB(4, 3);
  BB = B1 * B2;
  matrix<Type> B(3, 4);
  B = BB.transpose();
  REPORT(B);
  
  // first encounter
  matrix<Type> BE1(3, 3);
  BE1(0,0) = Type(0.0);
  BE1(0,1) = Type(1.0);
  BE1(0,2) = Type(0.0);
  BE1(1,0) = Type(0.0);
  BE1(1,1) = Type(0.0);
  BE1(1,2) = Type(1.0);
  BE1(2,0) = Type(1.0);
  BE1(2,1) = Type(0.0);
  BE1(2,2) = Type(0.0);
  
  matrix<Type> BE2(3, 4);
  BE2(0,0) = Type(1.0);
  BE2(0,1) = Type(0.0);
  BE2(0,2) = Type(0.0);
  BE2(0,3) = Type(0.0);
  BE2(1,0) = Type(0.0);
  BE2(1,1) = deltaNB;
  BE2(1,2) = Type(0.0);
  BE2(1,3) = Type(1.0)-deltaNB;
  BE2(2,0) = Type(0.0);
  BE2(2,1) = Type(0.0);
  BE2(2,2) = deltaB;
  BE2(2,3) = Type(1.0)-deltaB;
  
  matrix<Type> BEE(3,4);
  BEE = BE1 * BE2;
  matrix<Type> BE(4,3);
  BE = BEE.transpose();
  REPORT(BE);
  
  // prob of states at t+1 given states at t
  matrix<Type> A1(3, 3);
  A1(0,0) = phiNB;
  A1(0,1) = Type(0.0);
  A1(0,2) = Type(1.0)-phiNB;
  A1(1,0) = Type(0.0);
  A1(1,1) = phiB;
  A1(1,2) = Type(1.0)-phiB;
  A1(2,0) = Type(0.0);
  A1(2,1) = Type(0.0);
  A1(2,2) = Type(1.0);
  
  matrix<Type> A2(3, 3);
  A2(0,0) = Type(1.0)-psiNBB;
  A2(0,1) = psiNBB;
  A2(0,2) = Type(0.0);
  A2(1,0) = psiBNB;
  A2(1,1) = Type(1.0)-psiBNB;
  A2(1,2) = Type(0.0);
  A2(2,0) = Type(0.0);
  A2(2,1) = Type(0.0);
  A2(2,2) = Type(1.0);

  matrix<Type> A(3, 3);
  A = A1 * A2;
  REPORT(A);
  
  // init states
  vector<Type> PROP(3);
  PROP(0) = piNB;
  PROP(1) = Type(1.0)-piNB;
  PROP(2) = Type(0.0);
  REPORT(PROP);
  
  // likelihood
  Type ll;
  Type nll;
  array<Type> ALPHA(3);
  for (int i = 0; i < nh; i++) {
    int ei = fc(i)-1;
    vector<int> evennt = ch.col(i);
    ALPHA = PROP * vector<Type>(BE.row(fs(i))); // element-wise vector product
    for (int j = ei+1; j < km; j++) {
      ALPHA = multvecmat(ALPHA,A) * vector<Type>(B.row(evennt(j))); // vector matrix product, then element-wise vector product
    }
    ll += log(sum(ALPHA));
  }
  nll = -ll;
  return nll;
}

