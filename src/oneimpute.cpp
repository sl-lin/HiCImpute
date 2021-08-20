// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
  
  
  int matrow(int x, int y);
  
  std::vector<double> uppertriangular(arma::mat mm);
  std::vector<double> vec_nonzero(std::vector<double> v);
  std::vector<double> colmean(arma::mat m);
  
  double a_pdf(double a,double b,double delta,double pii,double 
                 delta_sigma);
  double beta_pdf(double a, double b, double x);
  double dnormal(double x);
  double findMean(arma::mat mm);
  double Gamma_pdf(double alpha, double beta, double x);
  double logitc(double t);
  double mu_pdf(double mu, std::vector<double> mu_k, double alpha, double beta, double mu_sigma);
  double mu_k_pdf(double mu_kk, double y_k, double lam_k, double mu, double mu_sigma, int s);
  double uppersum(arma::mat m);
  double vecmulti(std::vector<double> a);
  double vec_mean(std::vector<double> v);
  double vecsum(arma::vec a);
  double normalCDF(double value);
  double vec_sd(std::vector<double> v);
  
  // [[Rcpp::export]]
  std::vector<double> oneimpute(const int niter, const int burnin, const int n_single, arma::mat m, std::vector<double> lambda, std::vector<double> mu_sigma_vec,
                                std::vector<double> startval, int n, double epsilon1, double epsilon2, int ii){
    
    arma::Row<double> singlerow=arma::rowvec(m.row(ii));
    std::vector<double> cur;
    cur=startval;
    double mu_sig = mu_sigma_vec[ii]+0.001;
    double delta_sigma = 1.5;
    double Bu,alpha,mu_gamma,beta,mu,a,delta,b,pii;
    int s;
    std::vector<double> chain_sum(1+n_single, 0.0);
    std::vector<double> chain_mean(1+n_single, 0.0);
    
    for(int k=0; k<niter; k++){
      Bu = singlerow[n_single+3];
      alpha = cur[0];
      mu_gamma = ::Rf_runif(max(0.0, Bu-epsilon2), Bu+epsilon2);
      beta = mu_gamma/alpha;
      mu = cur[3];
      a = cur[4];
      delta = cur[5];
      b = cur[6];
      pii = cur[7];
      s = cur[8];
      std::vector<double>::const_iterator mu_k_first = cur.begin() + 9;              //get substd::vector
      std::vector<double>::const_iterator mu_k_last = cur.begin() + 9 + n_single;
      std::vector<double> mu_k(mu_k_first, mu_k_last);
      
      int alpha_accept, mu_accept, mu_k_accept, a_accept, delta_accept;
      
      // update alpha
      double current_alpha, proposed_alpha, current_beta, proposed_beta, r_alpha, u_alpha;
      current_alpha =  alpha;
      proposed_alpha = ::Rf_runif(1.0,1000.0);
      current_beta = beta;
      proposed_beta = mu_gamma/proposed_alpha;
      r_alpha = Gamma_pdf(proposed_alpha, proposed_beta, mu)/Gamma_pdf(current_alpha, current_beta, mu);
      u_alpha=::Rf_runif (0, 1);
      if(u_alpha < r_alpha){
        alpha = proposed_alpha;
        beta = proposed_beta;
        alpha_accept = 1;
      }else {
        alpha = current_alpha;
        beta = current_beta;
        alpha_accept = 0;
      }
      
      // update mu
      double current_mu, proposed_mu, r_mu, u_mu;
      current_mu = mu;
      do{
        proposed_mu = ::Rf_rnorm(current_mu, 1);
      } while (proposed_mu <= 0);
      
      r_mu = mu_pdf(proposed_mu, mu_k, alpha, beta, mu_sig)/mu_pdf(current_mu, mu_k, alpha, beta, mu_sig);
      u_mu = ::Rf_runif (0, 1);
      if(u_mu < r_mu){
        mu = proposed_mu;
        mu_accept = 1;
      }else {
        mu = current_mu;
        mu_accept = 0;
      }
      
      // updata mu_k
      double current_mu_k, proposed_mu_k, y_k, lam_k, r_mu_k, u_mu_k, muk;
      for(int q=0; q<n_single; q++){
        current_mu_k = mu_k[q];
        do{
          proposed_mu_k = ::Rf_rnorm(current_mu_k, 1);
        } while (proposed_mu_k <= 0);
        y_k = singlerow[q];
        lam_k = lambda[q];
        r_mu_k = mu_k_pdf(proposed_mu_k, y_k, lam_k, mu, mu_sig, s)/mu_k_pdf(current_mu_k, y_k, lam_k, mu, mu_sig, s);
        u_mu_k = ::Rf_runif(0,1);
        if(u_mu_k < r_mu_k){
          muk = proposed_mu_k;
          mu_k_accept = 1;
        }else {
          muk = current_mu_k;
          mu_k_accept = 0;
        }
        mu_k[q] = muk;
      }
      
      //update a
      double current_a, current_b, proposed_a, proposed_b, r_a, u_a;
      current_a = a;
      current_b = current_a*(1-delta)/delta;
      proposed_a = ::Rf_runif(1,1000);
      proposed_b = proposed_a*(1-delta)/delta;
      r_a = a_pdf(proposed_a, proposed_b, delta, pii, delta_sigma)/a_pdf(current_a, current_b, delta, pii, delta_sigma);
      u_a = ::Rf_runif(0,1);
      if(u_a < r_a){ //isnan(r_a) &&
        a = proposed_a;
        a_accept = 1;
      }else {
        a = current_a;
        a_accept = 0;
      }
      
      //update delta
      double current_delta, proposed_delta, r_delta, u_delta;
      current_delta = delta;
      proposed_delta = ::Rf_runif(max(0.0,singlerow[n_single+1]-epsilon1),min(1.0,singlerow[n_single+1]+epsilon1));
      r_delta = dnormal((logitc(proposed_delta)-logitc(a/(a+b)))/delta_sigma)/dnormal((logitc(current_delta)-logitc(a/(a+b)))/delta_sigma);
      u_delta = ::Rf_runif(0,1);
      if( u_delta < r_delta){
        delta = proposed_delta;
        delta_accept = 1;
      }else {
        delta = current_delta;
        delta_accept = 0;
      }
      
      // updata b
      b = a*(1-delta)/delta;
      
      // update pii
      pii = ::Rf_rbeta(s+a,1-s+b);
      
      // updata s
      double LAM=0, lam;
      for(int w=0; w<n_single; w++){
        if(singlerow[w]==0){
          lam=mu_k[w]*lambda[w];
          LAM += lam;
        }
      }
      s = ::Rf_rbinom(1,pii/(pii+(1-pii)*exp(-1*LAM)));
      
      // Combine data
      cur[0] = alpha;
      cur[1] = mu_gamma;
      cur[2] = beta;
      cur[3] = mu;
      cur[4] = a;
      cur[5] = delta;
      cur[6] = b;
      cur[7] = pii;
      cur[8] = s;
      for(int j=0; j<mu_k.size(); j++){
        cur[9+j] = mu_k[j];
      }
      
      // combine data
      if(k>(burnin)){
        chain_sum[0] += pii;
        for(int cc=0; cc<mu_k.size(); cc++){
          chain_sum[cc+1] += mu_k[cc];
        }
      }
    }// end of for loop of k
    
    //  output
    //std::vector<arma::mat> my_list={chain, acceptance};
    for (int r=0; r<chain_sum.size(); r++){
      chain_mean[r] = chain_sum[r]/(niter-burnin);
    }
    //arma::rowvec output(chain_mean);
    return chain_mean;
  }
  
  //Sum of upper triangular matrix
  double uppersum(arma::mat m){
    int rows = m.n_rows;
    double upper_sum = 0;
    int i,j;
    for (i = 0; i < rows; i++)
      for (j = 0; j < rows; j++) {
        if (i < j) {
          upper_sum += m(i,j);
        }
      }
      return upper_sum;
  }
  
  //Average of matrix
  double findMean(arma::mat mm) {
    int row = mm.n_rows;
    int col = mm.n_cols;
    double sum=0;
    for (int i=0; i<row; i++)
      for (int j=0; j<col; j++)
        sum += mm(i,j);
    return (double)sum/(row*col);
  }
  
  //add up elements in a std::vector
  double vecsum(arma::vec a){
    using std::begin;
    using std::end;
    double s = accumulate(a.begin(),a.end(),0);
    return s;
  }
  
  
  //upper triangular of matrix
  std::vector<double> uppertriangular(arma::mat mm){
    int row = mm.n_rows;
    int column = mm.n_cols;
    std::vector<double> k;
    for (int j = 0; j < column; j++)
    {
      for (int i = 0; i < row; i++)
      {
        if (i < j)
        {
          k.push_back(mm(i,j));
        }
      }
    }
    return k;
  }
  
  //find corresponding row of matrix cell
  int matrow(int x, int y){
    int r = x+1+y*(y-1)/2;
    return r;
  }
  
  // average of a std::vector
  double vec_mean(std::vector<double> v){
    auto n = v.size();
    float average = 0.0f;
    if ( n != 0) {
      average = accumulate( v.begin(), v.end(), 0.0) / n;
    }
    return average;
  }
  
  //standard deviation of a std::vector
  double vec_sd(std::vector<double> v){
    double E=0;
    double ave = vec_mean(v);
    double inverse = 1.0 / static_cast<double>(v.size()-1);
    for(unsigned int i=0;i<v.size();i++)
    {
      E += pow(static_cast<double>(v[i]) - ave, 2);
    }
    return sqrt(inverse * E);
  }
  
  //extract nonzero element of a std::vector
  std::vector<double> vec_nonzero(std::vector<double> v){
    std::vector<double> vv;
    for (int i = 0; i < v.size(); ++i){
      double c = v[i];
      if (c>0.0001) {
        vv.push_back(c);
      }
    }
    return vv;
  }
  
  //Normal pdf
  double dnormal(double x){
    return 1/(sqrt(2*M_PI))*exp(-1*pow(x,2)/2);
  }
  
  //Normal cdf
  double normalCDF(double value){
    return 0.5 * erfc(-value * M_SQRT1_2);
  }
  
  //logit function
  double logitc(double t){
    double tt = log(t/(1-t));
    return tt;
  }
  
  //Beta distribution pdf
  double beta_pdf(double a, double b, double x){
    double pdf = tgamma(a+b)/tgamma(a)/tgamma(b)*pow(x,a-1)*pow(1-x,b-1);
    return pdf;
  }
  
  //Gamma distribution pdf
  double Gamma_pdf(double alpha, double beta, double x){
    double pdf = 1/pow(beta,alpha)/tgamma(alpha)*pow(x,alpha-1)*exp(-x/beta);
    return pdf;
  }
  
  //pdf of a
  double a_pdf(double a,double b,double delta,double pii,double delta_sigma){
    double f=beta_pdf(a,b,pii)*exp(-1/(2*pow(delta_sigma,2))*pow((logitc(a/(a+b))-logitc(delta)),2))*pow(a+b,2)/(a*b);
    return f;
  }
  
  //multiply elements of a std::vector
  double vecmulti(std::vector<double> a){
    double multi = 1;
    for(int i = 0; i < a.size(); i++)
    {
      multi *= a[i];
    }
    return multi;
  }
  
  //pdf of mu
  double mu_pdf(double mu, std::vector<double> mu_k, double alpha, double beta, double mu_sigma){
    std::vector<double> mu_v;
    for(int i=0; i < mu_k.size(); i++){
      double f = dnormal((mu_k[i]-mu)/mu_sigma)/(1-normalCDF(-1*mu/mu_sigma));
      mu_v.push_back(f);
    }
    double v = vecmulti(mu_v)*pow(mu,(alpha-1))*exp(-mu/beta);
    return v;
  }
  
  //pdf of mu_k
  double mu_k_pdf(double mu_kk, double y_k, double lam_k, double mu, double mu_sigma, int s){
    double f;
    if(y_k>0){
      f=pow(mu_kk,y_k)*exp(-1*lam_k*mu_kk)*dnormal((mu_kk-mu)/mu_sigma);
    }else{
      f=exp(-1*lam_k*mu_kk*(1-s))*dnormal((mu_kk-mu)/mu_sigma);
    }
    return f;
  }
  
  //column mean
  std::vector<double> colmean(arma::mat m){
    std::vector<double> output;
    double mean;
    int rows=m.n_rows;
    for(int c=0; c<m.n_cols; c++){
      double sum=0;
      for(int r=0; r<m.n_rows; r++){
        sum += m(r,c);
      }
      mean = sum/rows;
      output.push_back(mean);
    }
    return output;
  }
 