#ifndef LBO_ITERATIVESOLVERS_H
#define LBO_ITERATIVESOLVERS_H

#include "itensor/all.h"
#include "lbo.hpp"
namespace itensor{
template<typename BigMatrixT>
  double
Applymatrix_with_LBO(BigMatrixT const& H, ITensor& v1, ITensor& w, 
		  Index const& id1, Index const& id2,int const& b, MPO& Ham, Args const& args)
{
auto [T1, T2]=applyLbo(v1, id1, id2,b, args);
   auto lbo_id1=findInds(T1, "OM");
   auto lbo_id2=findInds(T2, "OM");
   double new_dim(std::max(static_cast<double>(dim(lbo_id2)), static_cast<double>(dim(lbo_id1))));
auto Ham_temp1=Ham.ref(b);
auto Ham_temp2=Ham.ref(b+1);
v1*=T1;
v1*=T2;


Ham.setA(b, Ham.ref(b)*dag(T1));
Ham.setA(b+1, Ham.ref(b+1)*dag(T2));

        H.product(v1, w);
Ham.setA(b, Ham_temp1);
Ham.setA(b+1, Ham_temp2);

v1*=dag(T1);
v1*=dag(T2);
 return new_dim;
}

  template<typename BigMatrixT, typename ElT>
  std::tuple<double,double>
  Runge_Kutta_4_lbo(BigMatrixT const& H, ITensor& phi,
		 ElT const& tau, Index const& id1, Index const&  id2,int const&  b, MPO& Ham, Args const& args)
    {



    auto debug_level = args.getInt("DebugLevel",-1);

    // Initialize Lanczos vectors
    double lbo_max=0;
    double lbo_sum=0;
    double lbo_max_temp=0;
    ITensor v1 = phi;
    ITensor k1;
    ITensor k2;
    ITensor k3;
    ITensor k4;
    Real nrm = norm(v1);
    
    //generating k1

    lbo_max_temp=Applymatrix_with_LBO(H,v1, k1,id1, id2, b, Ham, args);
    lbo_sum+=lbo_max_temp;   
lbo_max=std::max(lbo_max, lbo_max_temp);
    
 // k2
      auto new_v=v1+0.5*tau*k1;
      lbo_max_temp=Applymatrix_with_LBO(H,new_v, k2,id1, id2, b, Ham, args);
    lbo_sum+=lbo_max_temp;   
    lbo_max=std::max(lbo_max, lbo_max_temp);

new_v=v1+0.5*tau*k2;

 lbo_max_temp=Applymatrix_with_LBO(H,new_v, k3,id1, id2, b, Ham, args);
    lbo_sum+=lbo_max_temp;   
    lbo_max=std::max(lbo_max, lbo_max_temp);

 new_v=v1+tau*k3;

 lbo_max_temp=Applymatrix_with_LBO(H,new_v, k4,id1, id2, b, Ham, args);
    lbo_sum+=lbo_max_temp;   
    lbo_max=std::max(lbo_max, lbo_max_temp);

   // updating the tensor y_{n+1}=y_n+dt/6(k_1+2k_2+2k_3+k4)
   ITensor next_T=phi+(tau/6)*(k1+2*k2+2*k3+k4);
 phi=std::move(next_T);
 // returning maximum and average
    return std::tuple<double,double>(lbo_max,lbo_sum/4);
    //std::tuple<lbo_max,lbo_sum/number_of_iterations>;
    } // End Runge Kutta lbo

   template<typename BigMatrixT, typename ElT>
void
  Runge_Kutta_4(BigMatrixT const& H, ITensor& phi,
		 ElT const&  tau,  Args const& args)
    {



    auto debug_level = args.getInt("DebugLevel",-1);

    // Initialize Lanczos vectors

    ITensor v1 = phi;
    ITensor k1;
    ITensor k2;
    ITensor k3;
    ITensor k4;
    Real nrm = norm(v1);
    
    //generating k1

        H.product(v1, k1);
      auto new_v=v1+0.5*tau*k1;
      H.product(new_v, k2);
      new_v=v1+0.5*tau*k2;
      H.product(new_v, k3);
      new_v=v1+tau*k3;
      H.product(new_v, k4);
   // updating the tensor y_{n+1}=y_n+dt/6(k_1+2k_2+2k_3+k4)
   ITensor next_T=phi+(tau/6)*(k1+2*k2+2*k3+k4);
 phi=std::move(next_T);
 // returning maximum and average
    return;
    //std::tuple<lbo_max,lbo_sum/number_of_iterations>;
    } // End Runge Kutta


template<typename BigMatrixT, typename ElT>
  std::tuple<double,double>
  applyExp_lbo(BigMatrixT const& H, ITensor& phi,
		 ElT tau, Index id1, Index id2,int b, MPO& Ham, Args const& args)
    {

    auto tol = args.getReal("ErrGoal",1E-10);
    auto max_iter = args.getInt("MaxIter",30);
    auto debug_level = args.getInt("DebugLevel",-1);
    auto beta_tol = args.getReal("NormCutoff",1e-7);

    // Initialize Lanczos vectors
    double lbo_max=0;
    double lbo_sum=0;
    int number_of_iterations=0;
    ITensor v1 = phi;
    ITensor v0;
    ITensor w;
    Real nrm = norm(v1);
    v1 /= nrm;


    std::vector<ITensor> lanczos_vectors({v1});
    Matrix bigTmat(max_iter + 2, max_iter + 2);
    std::fill(bigTmat.begin(), bigTmat.begin()+bigTmat.size(), 0.);

    auto nmatvec = 0;

    double beta = 0;
    for (int iter=0; iter < max_iter; ++iter)
        {
        int tmat_size=iter+1;
        // Matrix-vector multiplication
        if(debug_level >= 0)
            nmatvec++;

auto [T1, T2]=applyLbo(v1, id1, id2,b, args);
   auto lbo_id1=findInds(T1, "OM");
   auto lbo_id2=findInds(T2, "OM");
   double new_dim(std::max(static_cast<double>(dim(lbo_id2)), static_cast<double>(dim(lbo_id1))));
lbo_max=std::max(new_dim, lbo_max);
 lbo_sum+=new_dim;
 number_of_iterations+=1;
auto Ham_temp1=Ham.ref(b);
auto Ham_temp2=Ham.ref(b+1);
v1*=T1;
v1*=T2;


Ham.setA(b, Ham.ref(b)*dag(T1));
Ham.setA(b+1, Ham.ref(b+1)*dag(T2));

        H.product(v1, w);
Ham.setA(b, Ham_temp1);
Ham.setA(b+1, Ham_temp2);

v1*=dag(T1);
v1*=dag(T2);
        double avnorm = norm(w);
        double alpha = real(eltC(dag(w) * v1));
        bigTmat(iter, iter) = alpha;
        w -= alpha * v1;
        if (iter > 0)
            w -= beta * v0;
        v0 = v1;
        beta = norm(w);

        // check for Lanczos sequence exhaustion
        if (std::abs(beta) < beta_tol)
            {
            // Assemble the time evolved state
            auto tmat = subMatrix(bigTmat, 0,tmat_size, 0, tmat_size);
            auto tmat_exp = expMatrix(tmat, tau);
            auto linear_comb = column(tmat_exp, 0);
            assembleLanczosVectors(lanczos_vectors, linear_comb, nrm, phi);
            break;
            }

        // update next lanczos vector
        v1 = w;
        v1 /= beta;
        lanczos_vectors.push_back(v1);
        bigTmat(iter+1, iter) = beta;
        bigTmat(iter, iter+1) = beta;

        // Convergence check
        if (iter > 0)
            {
            // Prepare extended T-matrix for exponentiation
            int tmat_ext_size = tmat_size + 2;
            auto tmat_ext = Matrix(tmat_ext_size, tmat_ext_size);
            tmat_ext = subMatrix(bigTmat, 0, tmat_ext_size, 0, tmat_ext_size);

            tmat_ext(tmat_size-1, tmat_size) = 0.;
            tmat_ext(tmat_size+1, tmat_size) = 1.;

            // Exponentiate extended T-matrix
            auto tmat_ext_exp = expMatrix(tmat_ext, tau);

            double phi1 = std::abs( nrm*tmat_ext_exp(tmat_size, 0) );
            double phi2 = std::abs( nrm*tmat_ext_exp(tmat_size + 1, 0) * avnorm );
            double error;
            if (phi1 > 10*phi2) error = phi2;
            else if (phi1 > phi2) error = (phi1*phi2)/(phi1-phi2);
            else error = phi1;
            if(debug_level >= 1)
                println("Iteration: ", iter, ", Error: ", error);
            if ((error < tol) || (iter == max_iter-1))
                {
                if (iter == max_iter-1)
                    printf("warning: applyExp not converged in %d steps\n", max_iter);

                // Assemble the time evolved state
                auto linear_comb = Vec<ElT>(tmat_ext_size);
                linear_comb = column(tmat_ext_exp, 0);
                linear_comb = subVector(linear_comb, 0, tmat_ext_size-1);
                assembleLanczosVectors(lanczos_vectors, linear_comb, nrm, phi);
                if(debug_level >= 0)
                    printf("In applyExp, number of iterations: %d\n", iter);
                break;
                }
            }  // end convergence test

        }  // Lanczos iteratrions

    if(debug_level >= 0)
        println("In applyExp, number of matrix-vector multiplies: ", nmatvec);

    return std::tuple<double,double>(lbo_max,lbo_sum/number_of_iterations);
    //std::tuple<lbo_max,lbo_sum/number_of_iterations>;
    } // End applyExp LBO
    
template<typename BigMatrixT>
double
applyExp_lbo(BigMatrixT const& H, ITensor& phi,
         int tau, Args const& args)
    {
    return applyExp_lbo(H,phi,Real(tau),args);
    
    }

}
#endif /* LBO_ITERATIVESOLVERS_H */
