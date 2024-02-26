#include "parallel_tdvp.h"
#include "itensor/all.h"

using namespace itensor;
int
main(int argc, char* argv[])
    {
    Environment env(argc,argv);

    parallelDebugWait(env);

    int N = 20;

    SpinOne sites;
    MPO H;
    MPS psi;
    Sweeps sweeps;
    if(env.firstNode())
        {
        sites = SpinOne(N); //make a chain of N spin 1's
        auto ampo = AutoMPO(sites);
        for(int j = 1; j < N; ++j)
            {
            ampo += 0.5,"S+",j,"S-",j+1;
            ampo += 0.5,"S-",j,"S+",j+1;
            ampo +=     "Sz",j,"Sz",j+1;
            }
        H = MPO(ampo);
        auto state = InitState(sites);
        for(auto n : range1(N)) state.set(n,n%2==1?"Up":"Dn");
        psi = MPS(state);
        psi.normalize();

        sweeps = Sweeps(1);
        sweeps.maxdim() = 200;
        sweeps.cutoff() = 1E-10;
        sweeps.niter() = 30;
        println(sweeps);
        }
    double dt(-0.1);
  // env.broadcast(sites,H,psi,sweeps);
  //       Partition P;
  //   std::vector<ITensor> Vs;
  //     Args args = Args::global();
  //     // args.add("DoNormalize",false);
  //    splitWavefunction(env,psi,P,Vs,args);
  //     Observer obs;

  //      auto PH = computeHEnvironment(env,P,psi,Vs,H,args);
    
  //     ptdvpWorker(env,P,psi,Vs,PH,sweeps,obs,dt, args);
  //    parallel_tdvp(env,psi,H,sweeps,"Quiet");
    for(int i=0; i<10; i++)
      {
     env.broadcast(sites,H,psi,sweeps);
     
    	   	   

       Partition P;
     std::vector<ITensor> Vs;
     Args args = Args::global();
      args.add("DoNormalize",false);
     splitWavefunction(env,psi,P,Vs,args);
     Observer obs;

      auto PH = computeHEnvironment(env,P,psi,Vs,H,args);
    
      ptdvpWorker(env,P,psi,Vs,PH,sweeps,obs,dt, args);
     //     MPS psi2=psi;
     gather_vector(env,P, psi, Vs);
    	   if(env.firstNode())
    	    	     {
		       // print(psi);
    	   	       psi.position(1);
    	   	       std::cout<< "Norm "<< norm(psi)<<std::endl;
    	   	       psi/=norm(psi);

    	   	       std::cout<< "final energy "<<  real(innerC(psi,H,psi))<<std::endl;	       
    	   	     }
            }

  
    //  Args args={"Quiet","NumCenter",1};
    // for(int i=0; i<3; i++)
    //   {
    // 	 env.broadcast(sites,H,psi,sweeps);
    // parallel_tdvp(env,psi,H,sweeps);
    // parallel_tdvp(env,psi,H,sweeps);
    //   }
    // std::cout<< "final energy "<<  real(innerC(psi,H,psi))<<std::endl;
    // if(env.firstNode())
    //   {
    // 	printfln("There are %d nodes",env.nnodes());
    // 		std::cout<< "final energy "<<  real(innerC(psi,H,psi))<<std::endl;
    // 		}
    return 0;
    }
