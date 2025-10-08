#ifndef __ITENSOR_PARALLEL_TDVP
#define __ITENSOR_PARALLEL_TDVP

#include "itensor/mps/dmrg.h"
#include "itensor/util/parallel.h"
#include "partition.h"
#include "lbo_files/lbo.hpp"
#include"files.hpp"
#include"lbo_files/lbo_solvers.h"
#include"write_data.h"


namespace itensor {


double
gatherMax_val(Environment const& env, double& value)
    { 
    if(env.nnodes() == 1) return value;
    const int root = 0;
    
    if(env.firstNode())
        { 
        for(int i = 1; i < env.nnodes(); i++)
            {
            MailBox mailbox(env,i);
            double tmp;
            mailbox.receive(tmp);
	    if(tmp<0)
	      {std::cout<< "something wrong in gatherMax"<<std::endl;}
	    value=std::max(tmp, value);

            }
        }
    else
        {
        MailBox mailbox(env,root);
        mailbox.send(value);
	value=-1;
        }
    return value;
    }
  
Spectrum
isvd(ITensor T, 
     ITensor & A, 
     ITensor & V, 
     ITensor & B,
     Args const& args = Args::global())
    {
    auto pinv_cut = args.getReal("PInvCut",1E-12);
    auto spec = svd(T,A,V,B,args);
    A *= V;
    B *= V;
    auto pseudoInv = [pinv_cut](Real x)
        {
        return (std::fabs(x) <= pinv_cut) ? 0.0 : 1./x;
        };
    V.apply(pseudoInv);
    V.dag();
    return spec;
    }
  void gather_vector( Environment const& env,  Partition & P, MPS& psi,    std::vector<ITensor> & Vs)
  {

   
        
	   for(int b = 2; b <= P.Nb(); ++b)
                {
		  // for(auto j : range1(P.begin(b),P.end(b)))
		  //{
		   ITensor X;
		   bool edge=true;
		         for(auto j : range1(P.begin(b),P.end(b)))
		  {
	
	    	   if(env.rank()+1==b)
	    	     {

		       MailBox mbox(env, 0);	  
		       	 
	    	  X=psi.A(j);
 if(j==P.begin(b+1)-1)
		    {
		      X*=Vs.at( b);
		      // print(Vs.at(env.rank()));
	
		    }
		  else{
	
		  }
				  	    mbox.send(X);
		  // if(edge)
		  //   {
		  //       mbox.send(Vs.at(env.rank())*X);
		  //     //  std::cout << "V and "<<Vs.size()<<std::endl;
		  //     // print(Vs.at(env.rank()));
		  // edge=false;}
		  // else{
		  //   mbox.send(X);
		  // }
	    	 
	    	  }
		   	
	   
	    	   	   if(env.firstNode())
	    	     {
		   
		         MailBox mbox(env, b-1);
	
	      mbox.receive(X);
	 

		  if(j==P.begin(2))
		    {
		      psi.Aref(j)=Vs.at(b-1)*X;

	
		    }
		  else{
		     psi.Aref(j)=X;
		    //mbox.send(X);
		  }
	       //	      std::cout<< "new "<< psi.Aref(j)<<std::endl;
	      
	     
		   	
	    }
		  }
	}
  }
struct ParallelSweeper
    {
    int j = -1;
    int jl = -1;
    int jr = -1;
    int node = -1; //1-indexed
    int ha = 1; //records 1st or 2nd half sweep
    int sw = 0; //number of full sweeps done

    ParallelSweeper(int jl_,
                    int jr_,
                    int node_)
      : jl(jl_),
        jr(jr_),
        node(node_)
        { 
        if(node <= 0) Error("Node number out of range");
        j = odd() ? jl : (jr-1);
        }

    void
    operator++()
        {
        int inc = 0;
        int end = 0;
        if(odd())
            {
            inc = (ha==1 ? +1 : -1);
            end = (ha==1 ? jr : jl-1);
            }
        else
            {
            inc = (ha==1 ? -1 : +1);
            end = (ha==1 ? jl-1 : jr);
            }
        j += inc;
        if(j == end)
            {
            j -= inc; //went 1 too far, fix
            ++ha;
            }
        }

    void
    newSweep() { ha = 1; }
    bool
    doingHalf() const { return ha == 1; }
    bool
    doingFull() const { return ha < 3; }

    Direction
    dir() const
        {
        if(odd()) return (ha==1 ? Fromleft : Fromright);
        return (ha==1 ? Fromright : Fromleft);
        }

    bool
    atRight() const { return j == jr-1; }
    bool
    atLeft() const { return j == jl; }

    private:

    bool
    odd() const { return node%2==1; }

    };

 template < typename HamT, typename step_type>
Real 
ptdvpWorker(Environment const& env,
            Partition const& P,
            MPS & psi,
            std::vector<ITensor> & Vs,
            HamT & PH,
            Sweeps const& sweeps,
            Observer & obs,
	    step_type dt,
            Args args=Args::global());

LocalMPO 
computeHEnvironment(Environment const& env,
                    Partition const& P,
                    MPS const& psi,
                    std::vector<ITensor> const& Vs,
                    MPO const& H,
                    Args const& args = Args::global())
    {
    std::vector<ITensor> LH;
    std::vector<ITensor> RH;
    if(env.firstNode())
        {
        auto Nnode = env.nnodes();
        LH = std::vector<ITensor>(Nnode+1);
        RH = std::vector<ITensor>(Nnode+1);

        //Make left Hamiltonian environments
        auto L = ITensor(1.);
        LH.at(1) = L;
        for(int b = 1; b < P.Nb(); ++b)
            {
            for(auto j : range1(P.begin(b),P.end(b)))
                {
                L = L*psi.A(j)*H.A(j)*dag(prime(psi.A(j)));
                }
            L *= Vs.at(b);
            L *= dag(prime(Vs.at(b)));
            LH.at(b+1) = L;
            //printfln("LH[%d] = \n%s",b+1,LH[b+1]);
            //printfln("psi.A(%d) = \n%s",P.end(b)+1,psi.A(P.end(b)+1));
            }
        //Make right Hamiltonian environments
        auto R = ITensor(1.);
        RH.at(P.Nb()) = R;
        for(int b = P.Nb(); b > 1; --b)
            {
            for(auto j = P.end(b); j >= P.begin(b); --j)
                {
                R = R*psi.A(j)*H.A(j)*dag(prime(psi.A(j)));
                }
            R *= Vs.at(b-1);
            R *= dag(prime(Vs.at(b-1)));
            RH.at(b-1) = R;
            //printfln("RH[%d] = \n%s",b-1,RH[b-1]);
            }
        }
    env.broadcast(LH,RH);

    auto b = env.rank()+1; //block number of this node
    return LocalMPO(H,LH.at(b),P.begin(b)-1,RH.at(b),P.end(b)+1,args);
    }


LocalMPOSet 
computeHEnvironment(Environment const& env,
                    Partition const& P,
                    MPS const& psi,
                    std::vector<ITensor> const& Vs,
                    std::vector<MPO> const& Hset,
                    Args const& args = Args::global())
    {
    auto nset = Hset.size();
    auto Nnode = env.nnodes();
    std::vector<std::vector<ITensor>> LH(Nnode+1);
    std::vector<std::vector<ITensor>> RH(Nnode+1);
    if(env.firstNode())
        {
        for(auto& lh : LH) lh = std::vector<ITensor>(nset);
        for(auto& rh : RH) rh = std::vector<ITensor>(nset);

        for(auto n : range(nset))
            {
            auto& H = Hset.at(n);

            //Make left Hamiltonian environments
            auto L = ITensor(1.);
            LH.at(1).at(n) = L;
            for(int b = 1; b < P.Nb(); ++b)
                {
                for(auto j : range1(P.begin(b),P.end(b)))
                    {
                    L = L*psi.A(j)*H.A(j)*dag(prime(psi.A(j)));
                    }
                L *= Vs.at(b);
                L *= dag(prime(Vs.at(b)));
                LH.at(b+1).at(n) = L;
                //printfln("LHn[%d] = \n%s",b+1,LHn[b+1]);
                //printfln("psi.A(%d) = \n%s",P.end(b)+1,psi.A(P.end(b)+1));
                }
            //Make right Hamiltonian environments
            auto R = ITensor(1.);
            RH.at(P.Nb()).at(n) = R;
            for(int b = P.Nb(); b > 1; --b)
                {
                for(auto j = P.end(b); j >= P.begin(b); --j)
                    {
                    R = R*psi.A(j)*H.A(j)*dag(prime(psi.A(j)));
                    }
                R *= Vs.at(b-1);
                R *= dag(prime(Vs.at(b-1)));
                RH.at(b-1).at(n) = R;
                //printfln("RHn[%d] = \n%s",b-1,RHn[b-1]);
                }
            }
        }
    for(auto r : range(Nnode+1))
        {
        env.broadcast(LH.at(r),RH.at(r));
        }

    auto b = env.rank()+1; //block number of this node
    return LocalMPOSet(Hset,LH.at(b),P.begin(b)-1,RH.at(b),P.end(b)+1,args);
    }


void
splitWavefunction(Environment const& env,
                  MPS & psi, 
                  Partition & P,
                  std::vector<ITensor> & Vs,
                  Args const& args = Args::global())
    {
    if(env.firstNode()) 
        {
        auto Nnode = env.nnodes();
        if(args.defined("BoundarySize"))
            {
            P = Partition(psi.N(),Nnode,args.getInt("BoundarySize"));
            }
        else
            {
            P = Partition(psi.N(),Nnode);
            }
        println(P);

        Vs = std::vector<ITensor>(Nnode);
        psi.position(1);
        auto c = 1;
        for(int b = 1; b < P.Nb(); ++b)
            {
            auto n = P.end(b);
            //Shift ortho center to one past the end of the b'th block
            while(c < n+1)
                {
                ITensor D;
                svd(psi.A(c)*psi.A(c+1),psi.Aref(c),D,psi.Aref(c+1), args);
                psi.Aref(c+1) *= D;
                c += 1;
                }
            if(c != n+1) Error("c != n+1");
            auto AA = psi.A(n)*psi.A(n+1);
            auto& V = Vs.at(b);
            isvd(AA,psi.Aref(n),V,psi.Aref(n+1), args);
            }
        }
    env.broadcast(P,Vs,psi);
    }


//
// parallel_tdvp with single MPO or IQMPO
// and an observer object
//

Real 
parallel_tdvp(Environment const& env,
              MPS & psi,
              MPO const& H,
              Sweeps const& sweeps,
              Observer & obs,
              Args args = Args::global())
    {
    Partition P;
    std::vector<ITensor> Vs;
    splitWavefunction(env,psi,P,Vs,args);
    //  gather_vector(env,P, psi);
    auto PH = computeHEnvironment(env,P,psi,Vs,H,args);
    double dt=0.1;
      return ptdvpWorker(env,P,psi,Vs,PH,sweeps,obs,dt, args);
    }

//
// parallel_tdvp with single MPO or IQMPO
//

Real 
parallel_tdvp(Environment const& env,
              MPS & psi,
              MPO const& H,
              Sweeps const& sweeps,
              Args const& args = Args::global())
    {
    Observer obs;
    return parallel_tdvp(env,psi,H,sweeps,obs,args);
    }



//
// parallel_tdvp with an (implicit) sum of MPOs or IQMPOs
// and an observer object
//

Real 
parallel_tdvp(Environment const& env,
              MPS & psi,
              std::vector<MPO> const& Hset,
              Sweeps const& sweeps,
              Observer & obs,
              Args args = Args::global())
    {
    Partition P;

     std::vector<ITensor> Vs;
    splitWavefunction(env,psi,P,Vs,args);
    auto PH = computeHEnvironment(env,P,psi,Vs,Hset,args);
    double dt=0.1;
    return ptdvpWorker(env,P,psi,Vs,PH,sweeps,obs,dt, args);
    }

//
// parallel_tdvp with an (implicit) sum of MPOs or IQMPOs
//

Real 
parallel_tdvp(Environment const& env,
              MPS & psi,
              std::vector<MPO> const& Hset,
              Sweeps const& sweeps,
              Args const& args = Args::global())
    {
      
    Observer obs;
    return parallel_tdvp(env,psi,Hset,sweeps,obs,args);
    }

template< typename HType = ITensor>
struct Boundary
    {
    HType HH;
    ITensor A;
    ITensor UU;
    Real energy;

    Boundary() : energy(0) { }

    void
    write(std::ostream& s) const
        {
        itensor::write(s,HH);
        itensor::write(s,A);
        itensor::write(s,UU);
        itensor::write(s,energy);
        }
    void
    read(std::istream& s)
        {
        itensor::read(s,HH);
        itensor::read(s,A);
        itensor::read(s,UU);
        itensor::read(s,energy);
        }
    };

 template < typename HamT, typename step_type>
Real 
ptdvpWorker(Environment const& env,
            Partition const& P,
            MPS & psi,
            std::vector<ITensor> & Vs,
            HamT & PH,
            Sweeps const& sweeps,
            Observer & obs,
	    step_type t,
            Args args)
    {

    using EdgeType = stdx::decay_t<decltype(PH.L())>;
    //Maximum number of Davidson iterations at boundaries
    auto boundary_niter = args.getInt("BoundaryIter",sweeps.niter(1)+6);
    bool DoNormalize=args.getBool("DoNormalize",false);
    auto b = env.rank()+1;
    auto jl = P.begin(b);
    auto jr = P.end(b);

    Real energy = 0.;
  ITensor phi0,phi1;
    auto psw = ParallelSweeper(jl,jr,b);

    MailBox mboxL,mboxR;
    if(not env.firstNode()) mboxL = MailBox(env,env.rank()-1);
    if(not env.lastNode())  mboxR = MailBox(env,env.rank()+1);

    psi.leftLim(jl);
    psi.rightLim(jr);
    psi.position(psw.j);

    //Include MPINode number to use in the observer
    args.add("BlockStart",jl);
    args.add("BlockEnd",jr);
    args.add("MPINode",env.rank()+1);
    int jmid = (jr-jl)/2.;
    const int numCenter = 2;
    //args.getInt("NumCenter",2);
     if(numCenter != 2)
       {std::cout<< "does not work for tdvp1 yet!!"<<std::endl;
	 return 0;}
    // if(numCenter != 1)
    //     args.add("Truncate",args.getBool("Truncate",true));
    // else
    //     args.add("Truncate",args.getBool("Truncate",false));
    //  args.add("DoNormalize",true);
    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
        {
        args.add("Sweep",sw);
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("Noise",sweeps.noise(sw));
        args.add("MaxIter",sweeps.niter(sw));

        if(!PH.doWrite()
           && args.defined("WriteM")
           && sweeps.maxdim(sw) >= args.getInt("WriteM"))
            {
            printfln("\nNode %d turning on write to disk, write_dir = %s",
                     b,args.getString("WriteDir","./"));
            //psi.doWrite(true);
            PH.doWrite(true);
            }

           printfln("Doing sweep %d for node %d (maxm=%d, cutoff=%.0E, mindim=%d)",sw,b,sweeps.maxdim(sw),sweeps.cutoff(sw),sweeps.mindim(sw));

         for(psw.newSweep(); psw.doingFull(); ++psw)
             {
            auto j = psw.j;
            auto dir = psw.dir();
            //printfln("%d j = %d (%d,%d) %s",b,j,jl,jr,dir==Fromleft?"Fromleft":"Fromright");
	    // if(env.firstNode())
	    //   {
	    //	    std::cout<< "node rank " << env.rank()+1<<" act on site  "<< j<< " and "<< j+1<< " and ha = "<<psw.ha <<std::endl;
		//	      }
	      //   if(env.lastNode())
	      // {
	      // 	std::cout<< "node las "<< j <<std::endl;
	      // }
	        PH.numCenter(numCenter);
            PH.position(j,psi);

             phi1 = psi.A(j)*psi.A(j+1);

	     //  energy = davidson(PH,phi1,args);
	      applyExp(PH,phi1,t/2,args);
            //if(env.rank()+1 == 1) printfln("%s j = %d energy = %.10f",dir==Fromleft?"->":"<-",j,energy);
		 if(DoNormalize)
	      
		   {
phi1 /= norm(phi1);
		   }
	     auto energy=0;
            auto spec = psi.svdBond(j,phi1,dir,PH,args);
            
            if(env.rank()+1 == env.nnodes()/2 
            && dir == Fromright 
            && j == jmid)
                {
                printfln("Truncation error for sweep %d (node %d) at site %d is terr=%.12f",
                         sw,b,j,spec.truncerr());
                }

            args.add("AtBond",j);
            args.add("AtBoundary",false);
            args.add("Energy",energy);
            args.add("Direction",dir);
            obs.measure(args);

	     if((psw.ha == 1 && psw.j+numCenter-1 != length(psi) && dir==Fromleft) || (psw.ha == 2 && psw.j != 1 && dir==Fromright))

                {
                auto b1 = (psw.ha == 1 ? psw.j+1 : psw.j);
		//	std::cout<< " XXX   rank "<<env.rank()+1<< " and j " << psw.j<<"  and "<< b1 << std::endl;
		//     if(numCenter == 2)
		//  {
		 phi0 = psi.A(b1);
		    //  }
		 //			 std::cout<< "j "<< j << " b1 "<<b1 <<std::endl;
			 // PH.position(b1,psi);
			 	      PH.numCenter(numCenter-1);
				       PH.position(b1,psi);
 
				       applyExp(PH,phi0,-t/2,args);
				       //					  energy = davidson(PH,phi0,args);
                 if(DoNormalize)
		   {		       	          phi0 /= norm(phi0);
		   }
                // if(numCenter == 2)
                //     {
		        psi.ref(b1) = phi0;
                //     }
		}
	     else if((psw.ha == 1 && psw.j != 1  && dir==Fromright) || (psw.ha == 2 && psw.j+numCenter-1 != length(psi) && dir==Fromleft))
                 {
		   auto b1 = (psw.ha == 1 ?  psw.j : psw.j+1);
		   //  	std::cout<< " YYY  rank "<<env.rank()+1<< " and j " << psw.j<<"  and "<< b1 << std::endl;
		//          if(numCenter == 2)
                //     {
		    phi0 = psi(b1);
                //     }
          
 
		     PH.numCenter(numCenter-1);
		       PH.position(b1,psi);
 
		           applyExp(PH,phi0,-t/2,args);
			   //energy = davidson(PH,phi0,args);
                 if(DoNormalize)
		   {	    phi0 /= norm(phi0);
		   }
                
                // if(numCenter == 2)
                //     {
		        psi.ref(b1) = phi0;
                //     }
  
	       	}
		
            if(psw.atRight() && dir==Fromleft && bool(mboxR))
                {
		  //	  std::cout<< "WWW rank "<<env.rank()+1<< " j "<<j << " n " << j+1<<std::endl;
    //             auto prev_energy = energy;
    //             printfln("Node %d communicating with right, boundary_niter=%d",b,boundary_niter);
     				auto n = j+1;
					std::cout<< " n hete "<< n << " at  "<<env.rank() <<std::endl;
  PH.numCenter(numCenter);
        PH.position(n,psi);

                 Boundary<EdgeType> B;
                 mboxR.receive(B);
                 psi.Aref(n+1) = B.A;
                 PH.R(B.HH);
                 B = Boundary<EdgeType>(); //to save memory

                 auto& V = Vs.at(b);

                auto phi = psi.A(n)*V*psi.A(n+1);
		//  
		//	std::cout<< "Pri 11111"<< V<<std::endl;
		//energy = davidson(PH,phi,{args,"MaxIter=",boundary_niter});
			applyExp(PH,phi,t,args);
		//energy = davidson(PH,phi,args);
			//			std::cout<< "norm s "<< norm(phi)<<std::endl; 
if(DoNormalize)
		  {
				phi /= norm(phi);
		  }
			//	std::cout<< "norm s "<< norm(phi)<<std::endl;
				auto spec = isvd(phi,psi.Aref(n),V,psi.Aref(n+1), args);
				//	std::cout<< "Pri "<< V<<std::endl;

		
                B.HH = PH.L();
                B.UU = psi.A(n)*V;
                B.A = psi.A(n+1);
                B.energy = energy;
                mboxR.send(B);

                psi.Aref(n+1) *= V;
                psi.rightLim(n+1);





		// new part
		phi0 = psi.A(n);
		    //  }
		 //			 std::cout<< "j "<< j << " b1 "<<b1 <<std::endl;
			 // PH.position(b1,psi);
			 	      PH.numCenter(numCenter-1);
				       PH.position(n,psi);
 
				       applyExp(PH,phi0,-t/2,args);
				       //					  energy = davidson(PH,phi0,args);
                if(DoNormalize)
		  {
		       		        	          phi0 /= norm(phi0);
                
		  }
                // if(numCenter == 2)
                //     {
		        psi.ref(n) = phi0;




                // args.add("AtBond",n);
                // args.add("AtBoundary",true);
                // args.add("Energy",energy);
                // obs.measure(args);

		//  printfln("Node %d done with boundary step, energy %.5f -> %.5f",b,prev_energy,energy);
                 }
             else if(psw.atLeft() && dir==Fromright && bool(mboxL))
                 {
		   //	  std::cout<< "LLL rank "<<env.rank()+1<< " j "<<j << " n " << j-1<<std::endl;
			    PH.numCenter(numCenter);
			    //	      PH.position(n,psi);
    //             auto prev_energy = energy;
    //             printfln("Node %d communicating with left, boundary_niter=%d",b,boundary_niter);
			    
                auto n = j-1;
	
  PH.numCenter(numCenter);
                PH.position(n,psi);

                Boundary<EdgeType> B;
                B.A = psi.A(n+1);
                B.HH = PH.R();
                mboxL.send(B);

                mboxL.receive(B);
                PH.L(B.HH);
                PH.shift(n,Fromleft,B.UU);
                psi.Aref(n+1) = B.A;
                energy = B.energy;
		phi0 = psi.A(n+1);
      PH.numCenter(numCenter-1);
      PH.position(n+1,psi);
      
      applyExp(PH,phi0,-t/2,args);
       if(DoNormalize)
		  {
	       	          phi0 /= norm(phi0);
		  }
                // if(numCenter == 2)
                //     {
		        psi.ref(n+1) = phi0;
		// printfln("Node %d done with boundary step, energy %.5f -> %.5f",b,prev_energy,energy);
                 }
             }

    //     if(obs.checkDone(args)) break;
           }

    printfln("Block %d final energy = %.12f",b,energy);

    return energy;

    } // ptdvpWorker

   template < typename HamT, typename step_type>
void 
Act(Environment const& env,
            Partition const& P,
            MPS & psi,
            std::vector<ITensor> & Vs,
            HamT & PH,
            Sweeps const& sweeps,
            Observer & obs,
	    step_type t,
    Args args, bool Testing)
    {

    using EdgeType = stdx::decay_t<decltype(PH.L())>;
    //Maximum number of Davidson iterations at boundaries
    auto boundary_niter = args.getInt("BoundaryIter",sweeps.niter(1)+6);
    bool DoNormalize=args.getBool("DoNormalize",false);
    auto b = env.rank()+1;
    auto jl = P.begin(b);
    auto jr = P.end(b);

    Real energy = 0.;
  ITensor phi0,phi1;
    auto psw = ParallelSweeper(jl,jr,b);

    MailBox mboxL,mboxR;
    if(not env.firstNode()) mboxL = MailBox(env,env.rank()-1);
    if(not env.lastNode())  mboxR = MailBox(env,env.rank()+1);
    if(Testing)
      {
    psi.leftLim(jl);
    psi.rightLim(jr);
    Testing=false;
    psi.position(psw.j);
      }


    //Include MPINode number to use in the observer
    args.add("BlockStart",jl);
    args.add("BlockEnd",jr);
    args.add("MPINode",env.rank()+1);
    int jmid = (jr-jl)/2.;
    const int numCenter = 2;
    //args.getInt("NumCenter",2);
     if(numCenter != 2)
       {std::cout<< "does not work for tdvp1 yet!!"<<std::endl;
	 return;}
    // if(numCenter != 1)
    //     args.add("Truncate",args.getBool("Truncate",true));
    // else
    //     args.add("Truncate",args.getBool("Truncate",false));
    //  args.add("DoNormalize",true);
    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
        {
 args.add("Sweep",sw);
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("Noise",sweeps.noise(sw));
        args.add("MaxIter",sweeps.niter(sw));

        // if(!PH.doWrite()
        //    && args.defined("WriteM")
        //    && sweeps.maxdim(sw) >= args.getInt("WriteM"))
        //     {
        //     printfln("\nNode %d turning on write to disk, write_dir = %s",
        //              b,args.getString("WriteDir","./"));
        //     //psi.doWrite(true);
        //     PH.doWrite(true);
        //     }

	//    printfln("Doing sweep %d for node %d (maxm=%d, cutoff=%.0E, mindim=%d)",sw,b,sweeps.maxdim(sw),sweeps.cutoff(sw),sweeps.mindim(sw));

         for(psw.newSweep(); psw.doingFull(); ++psw)
             {
            auto j = psw.j;
            auto dir = psw.dir();

	        PH.numCenter(numCenter);
            PH.position(j,psi);

             phi1 = psi.A(j)*psi.A(j+1);


	      applyExp(PH,phi1,t/2,args);
      
		 if(DoNormalize)
	      
		   {
phi1 /= norm(phi1);
		   }
	     auto energy=0;
            auto spec = psi.svdBond(j,phi1,dir,PH,args);
            
            if(env.rank()+1 == env.nnodes()/2
            && dir == Fromright
            && j == jmid)
                {
                printfln("Truncation error for sweep %d (node %d) at site %d is terr=%.12f",
                         sw,b,j,spec.truncerr());
                }

            args.add("AtBond",j);
            args.add("AtBoundary",false);
            args.add("Energy",energy);
            args.add("Direction",dir);
            obs.measure(args);

	     if((psw.ha == 1 && psw.j+numCenter-1 != length(psi) && dir==Fromleft) || (psw.ha == 2 && psw.j != 1 && dir==Fromright))

                {
                auto b1 = (psw.ha == 1 ? psw.j+1 : psw.j);
	
		 phi0 = psi.A(b1);

			 	      PH.numCenter(numCenter-1);
				       PH.position(b1,psi);
 
				       applyExp(PH,phi0,-t/2,args);
				
                 if(DoNormalize)
		   {		       	          phi0 /= norm(phi0);
		   }
              
		        psi.ref(b1) = phi0;
              
		}
	     else if((psw.ha == 1 && psw.j != 1  && dir==Fromright) || (psw.ha == 2 && psw.j+numCenter-1 != length(psi) && dir==Fromleft))
                 {
		   auto b1 = (psw.ha == 1 ?  psw.j : psw.j+1);
	
		    phi0 = psi(b1);
        
          
 
		     PH.numCenter(numCenter-1);
		       PH.position(b1,psi);
 
		           applyExp(PH,phi0,-t/2,args);
		
                 if(DoNormalize)
		   {	    phi0 /= norm(phi0);
		   }
                
      
		        psi.ref(b1) = phi0;
  
         
	       	}
		
            if(psw.atRight() && dir==Fromleft && bool(mboxR))
                {
	

     				auto n = j+1;

  PH.numCenter(numCenter);
  //std
          PH.position(n,psi);

                  Boundary<EdgeType> B; 
                  mboxR.receive(B); 
                 psi.Aref(n+1) = B.A; 
                  PH.R(B.HH); 
                  B = Boundary<EdgeType>(); //to save memory 

               auto& V = Vs.at(b); 

                 auto phi = psi.A(n)*V*psi.A(n+1); 
	
 			applyExp(PH,phi,t,args); 
		
/* if(DoNormalize) */
/* 		  { */
/* 				phi /= norm(phi); */
/* 		  } */
		
 				auto spec = isvd(phi,psi.Aref(n),V,psi.Aref(n+1), args); 
			

		
                B.HH = PH.L();
                B.UU = psi.A(n)*V;
                B.A = psi.A(n+1);
                B.energy = energy;
                mboxR.send(B);

                psi.Aref(n+1) *= V;
                psi.rightLim(n+1);





	
		phi0 = psi.A(n);
	
			 	      PH.numCenter(numCenter-1);
				       PH.position(n,psi);
 
				       applyExp(PH,phi0,-t/2,args);
				       //
                if(DoNormalize)
		  {
		       		        	          phi0 /= norm(phi0);
                
		  }
               
		        psi.ref(n) = phi0;
                  } 
              else if(psw.atLeft() && dir==Fromright && bool(mboxL)) 
                  { 


 			    PH.numCenter(numCenter); 
		
			    
                 auto n = j-1; 
	
   PH.numCenter(numCenter); 
                 PH.position(n,psi); 

                Boundary<EdgeType> B;
                B.A = psi.A(n+1);
                B.HH = PH.R();
                mboxL.send(B);

                mboxL.receive(B);
                PH.L(B.HH);
                PH.shift(n,Fromleft,B.UU);
                psi.Aref(n+1) = B.A;
                energy = B.energy;
		phi0 = psi.A(n+1);
      PH.numCenter(numCenter-1);
      PH.position(n+1,psi);
      
      applyExp(PH,phi0,-t/2,args);
       if(DoNormalize)
		  {
	       	          phi0 /= norm(phi0);
		  }
                // if(numCenter == 2)
                //     {
		        psi.ref(n+1) = phi0;
		// printfln("Node %d done with boundary step, energy %.5f -> %.5f",b,prev_energy,energy);
                  } 
              } 

    //     if(obs.checkDone(args)) break;
           }

    //  printfln("Block %d final energy = %.12f",b,energy);

    return;

    } // ptdvpWorker


   template < typename HamT, typename step_type>
double
Act_lbo(Environment const& env,
            Partition const& P,
            MPS & psi,
            std::vector<ITensor> & Vs,
            HamT & PH,
            Sweeps const& sweeps,
            Observer & obs,
	    step_type t,
    Args args, bool Testing)
    {

    using EdgeType = stdx::decay_t<decltype(PH.L())>;
    //Maximum number of Davidson iterations at boundaries
    auto boundary_niter = args.getInt("BoundaryIter",sweeps.niter(1)+6);
    bool DoNormalize=args.getBool("DoNormalize",false);
    auto b = env.rank()+1;
    auto jl = P.begin(b);
    auto jr = P.end(b);
    double lbo_max=0;
    Real energy = 0.;
  ITensor phi0,phi1;
    auto psw = ParallelSweeper(jl,jr,b);

    MailBox mboxL,mboxR;
    if(not env.firstNode()) mboxL = MailBox(env,env.rank()-1);
    if(not env.lastNode())  mboxR = MailBox(env,env.rank()+1);
    if(Testing)
      {
    psi.leftLim(jl);
    psi.rightLim(jr);
    psi.position(psw.j);
    Testing=false;
      }


    //Include MPINode number to use in the observer
    args.add("BlockStart",jl);
    args.add("BlockEnd",jr);
    args.add("MPINode",env.rank()+1);
    int jmid = (jr-jl)/2.;
    const int numCenter = 2;
    //args.getInt("NumCenter",2);
     if(numCenter != 2)
       {std::cout<< "does not work for tdvp1 yet!!"<<std::endl;
	 return 0;}

    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
        {
        args.add("Sweep",sw);
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("Noise",sweeps.noise(sw));
        args.add("MaxIter",sweeps.niter(sw));

         for(psw.newSweep(); psw.doingFull(); ++psw)
             {
            auto j = psw.j;
            auto dir = psw.dir();

	        PH.numCenter(numCenter);
            PH.position(j,psi);

             phi1 = psi.A(j)*psi.A(j+1);


	      applyExp(PH,phi1,t/2,args);
	      	      	auto indx1=findIndex(psi.A(j), "Site");
	      		  auto indx2=findIndex(psi.A(j+1), "Site");
		  	         int old_llim=psi.leftLim();
		  int old_rlim=psi.rightLim();
	      auto [T1, T2]=applyLbo(phi1, indx1, indx2,j, args);
	       auto id1=findInds(T1, "OM");
	       auto id2=findInds(T2, "OM");
	       double new_dim(std::max(static_cast<double>(dim(id2)), static_cast<double>(dim(id1))));
	       lbo_max=std::max(new_dim, lbo_max);

	  phi1*=T1;
	  phi1*=T2;
	  psi.set(j , psi.A(j)*T1);
	  psi.Aref(j+1)*=T2;
	  	  psi.leftLim(old_llim);
	  		  psi.rightLim(old_rlim);
		 if(DoNormalize)
	      
		   {
phi1 /= norm(phi1);
		   }
	     auto energy=0;
            auto spec = psi.svdBond(j,phi1,dir,PH,args);
   old_llim=psi.leftLim();
    old_rlim=psi.rightLim();
            	  psi.Aref(j)*=dag(T1);
   		 psi.Aref(j+1)*=dag(T2);
		 
   		 	  psi.leftLim(old_llim);
   		  psi.rightLim(old_rlim);
            if(env.rank()+1 == env.nnodes()/2 
            && dir == Fromright 
            && j == jmid)
                {
                printfln("Truncation error for sweep %d (node %d) at site %d is terr=%.12f",
                         sw,b,j,spec.truncerr());
                }

            args.add("AtBond",j);
            args.add("AtBoundary",false);
            args.add("Energy",energy);
            args.add("Direction",dir);
            obs.measure(args);

	     if((psw.ha == 1 && psw.j+numCenter-1 != length(psi) && dir==Fromleft) || (psw.ha == 2 && psw.j != 1 && dir==Fromright))

                {
                auto b1 = (psw.ha == 1 ? psw.j+1 : psw.j);
	
		 phi0 = psi.A(b1);

			 	      PH.numCenter(numCenter-1);
				       PH.position(b1,psi);
 
				       applyExp(PH,phi0,-t/2,args);
				
                 if(DoNormalize)
		   {		       	          phi0 /= norm(phi0);
		   }
              
		        psi.ref(b1) = phi0;
              
		}
	     else if((psw.ha == 1 && psw.j != 1  && dir==Fromright) || (psw.ha == 2 && psw.j+numCenter-1 != length(psi) && dir==Fromleft))
                 {
		   auto b1 = (psw.ha == 1 ?  psw.j : psw.j+1);
	
		    phi0 = psi(b1);
        
          
 
		     PH.numCenter(numCenter-1);
		       PH.position(b1,psi);
 
		           applyExp(PH,phi0,-t/2,args);
		
                 if(DoNormalize)
		   {	    phi0 /= norm(phi0);
		   }
                
      
		        psi.ref(b1) = phi0;
  
         
	       	}
		
            if(psw.atRight() && dir==Fromleft && bool(mboxR))
                {
	
     				auto n = j+1;

  PH.numCenter(numCenter);
        PH.position(n,psi);

                 Boundary<EdgeType> B;
                 mboxR.receive(B);
                 psi.Aref(n+1) = B.A;
                 PH.R(B.HH);
                 B = Boundary<EdgeType>(); //to save memory

                 auto& V = Vs.at(b);
		 //////////////////////////////////
		 //step one
		 /////////////////
                auto phi = psi.A(n)*V*psi.A(n+1);
	
			applyExp(PH,phi,t/2,args);
		
if(DoNormalize)
		  {
				phi /= norm(phi);
		  }
  	auto indx1=findIndex(psi.A(n), "Site");
		  auto indx2=findIndex(psi.A(n+1), "Site");
		  
		    			         int old_llim=psi.leftLim();
		  int old_rlim=psi.rightLim();
   auto [T1, T2]=applyLbo(phi, indx1, indx2,n, args);
   auto id1=findInds(T1, "OM");
   auto id2=findInds(T2, "OM");
   double new_dim(std::max(static_cast<double>(dim(id2)), static_cast<double>(dim(id1))));
    	       lbo_max=std::max(new_dim, lbo_max);
   	      phi*=T1;
   	  phi*=T2;
   	  psi.set(n , psi.A(n)*T1);
   	  psi.Aref(n+1)*=T2;
	         	   	  psi.leftLim(old_llim);
	       	  psi.rightLim(old_rlim);

		  // auto spec = isvd(phi,psi.Aref(n),V,psi.Aref(n+1), args);
 auto spec = isvd(phi,psi.Aref(n),V,psi.Aref(n+1), args);
			
     	   old_llim=psi.leftLim();
     old_rlim=psi.rightLim();
    
     		  psi.Aref(n)*=dag(T1);
     		  psi.Aref(n+1)*=dag(T2);
		 
     		  	  psi.leftLim(old_llim);
     		   psi.rightLim(old_rlim);
		   //////////step two //////////
                 phi = psi.A(n)*V*psi.A(n+1);
	
			applyExp(PH,phi,t/2,args);
		
if(DoNormalize)
		  {
				phi /= norm(phi);
		  }
  	 indx1=findIndex(psi.A(n), "Site");
		   indx2=findIndex(psi.A(n+1), "Site");
		  
		    			         old_llim=psi.leftLim();
		   old_rlim=psi.rightLim();
    auto [T1_2, T2_2]=applyLbo(phi, indx1, indx2,n, args);
    id1=findInds(T1_2, "OM");
   id2=findInds(T2_2, "OM");
   new_dim=(std::max(static_cast<double>(dim(id2)), static_cast<double>(dim(id1))));
    	       lbo_max=std::max(new_dim, lbo_max);
   	      phi*=T1_2;
   	  phi*=T2_2;
   	  psi.set(n , psi.A(n)*T1_2);
   	  psi.Aref(n+1)*=T2_2;
	         	   	  psi.leftLim(old_llim);
	       	  psi.rightLim(old_rlim);

		  // auto spec = isvd(phi,psi.Aref(n),V,psi.Aref(n+1), args);
  spec = isvd(phi,psi.Aref(n),V,psi.Aref(n+1), args);
			
     	   old_llim=psi.leftLim();
     old_rlim=psi.rightLim();
    
     		  psi.Aref(n)*=dag(T1_2);
     		  psi.Aref(n+1)*=dag(T2_2);
		 
     		  	  psi.leftLim(old_llim);
     		   psi.rightLim(old_rlim);
		   /////////////////////
                B.HH = PH.L();
                B.UU = psi.A(n)*V;
                B.A = psi.A(n+1);
                B.energy = energy;
                mboxR.send(B);

                psi.Aref(n+1) *= V;
                psi.rightLim(n+1);





	
		phi0 = psi.A(n);
	
			 	      PH.numCenter(numCenter-1);
				       PH.position(n,psi);
 
				       applyExp(PH,phi0,-t/2,args);
				       //				
                if(DoNormalize)
		  {
		       		        	          phi0 /= norm(phi0);
                
		  }
               
		        psi.ref(n) = phi0;
                 }
             else if(psw.atLeft() && dir==Fromright && bool(mboxL))
                 {
		
			    PH.numCenter(numCenter);
		
			    
                auto n = j-1;
	
  PH.numCenter(numCenter);
                PH.position(n,psi);

                Boundary<EdgeType> B;
                B.A = psi.A(n+1);
                B.HH = PH.R();
                mboxL.send(B);

                mboxL.receive(B);
                PH.L(B.HH);
                PH.shift(n,Fromleft,B.UU);
                psi.Aref(n+1) = B.A;
                energy = B.energy;
		phi0 = psi.A(n+1);
      PH.numCenter(numCenter-1);
      PH.position(n+1,psi);
      
      applyExp(PH,phi0,-t/2,args);
       if(DoNormalize)
		  {
	       	          phi0 /= norm(phi0);
		  }
     
		        psi.ref(n+1) = phi0;

                 }
             }

           }



    return lbo_max;

    } // ptdvpWorker

   template < typename HamT, typename step_type>
int 
Act_lboBT(Environment const& env,
            Partition const& P,
            MPS & psi,
            std::vector<ITensor> & Vs,
            HamT & PH,
            Sweeps const& sweeps,
            Observer & obs,
	    step_type t,
    Args args, bool Testing)
    {

    using EdgeType = stdx::decay_t<decltype(PH.L())>;
    //Maximum number of Davidson iterations at boundaries
    auto boundary_niter = args.getInt("BoundaryIter",sweeps.niter(1)+6);
    bool DoNormalize=args.getBool("DoNormalize",false);
    auto b = env.rank()+1;
    auto jl = P.begin(b);
    auto jr = P.end(b);
    double lbo_max=0;
    Real energy = 0.;
  ITensor phi0,phi1;
    auto psw = ParallelSweeper(jl,jr,b);

    MailBox mboxL,mboxR;
    if(not env.firstNode()) mboxL = MailBox(env,env.rank()-1);
    if(not env.lastNode())  mboxR = MailBox(env,env.rank()+1);
    if(Testing)
      {
    psi.leftLim(jl);
    psi.rightLim(jr);
    psi.position(psw.j);
    Testing=false;
      }


    //Include MPINode number to use in the observer
    args.add("BlockStart",jl);
    args.add("BlockEnd",jr);
    args.add("MPINode",env.rank()+1);
    int jmid = (jr-jl)/2.;
    const int numCenter = 2;
    //args.getInt("NumCenter",2);
     if(numCenter != 2)
       {std::cout<< "does not work for tdvp1 yet!!"<<std::endl;
	 return 0;}

    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
        {
        args.add("Sweep",sw);
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("Noise",sweeps.noise(sw));
        args.add("MaxIter",sweeps.niter(sw));

         for(psw.newSweep(); psw.doingFull(); ++psw)
             {
            auto j = psw.j;
            auto dir = psw.dir();

	        PH.numCenter(numCenter);
            PH.position(j,psi);

             phi1 = psi.A(j)*psi.A(j+1);


	      applyExp(PH,phi1,t/2,args);
	      	      	auto indx1=findIndex(psi.A(j), "Site");
	      		  auto indx2=findIndex(psi.A(j+1), "Site");
		  	         int old_llim=psi.leftLim();
		  int old_rlim=psi.rightLim();
	      auto [T1, T2]=applyLbo(phi1, indx1, indx2,j, args);
	       auto id1=findInds(T1, "OM");
	       auto id2=findInds(T2, "OM");
	       double new_dim(std::max(static_cast<double>(dim(id2)), static_cast<double>(dim(id1))));
	       lbo_max=std::max(new_dim, lbo_max);

	  phi1*=T1;
	  phi1*=T2;
	  psi.Aref(j)*=T1;
	  
	  psi.Aref(j+1)*=T2;
	  	  psi.leftLim(old_llim);
	  		  psi.rightLim(old_rlim);
		 if(DoNormalize)
	      
		   {
phi1 /= norm(phi1);
		   }
	     auto energy=0;
            auto spec = psi.svdBond(j,phi1,dir,PH,args);
   old_llim=psi.leftLim();
    old_rlim=psi.rightLim();
            	  psi.Aref(j)*=dag(T1);
   		 psi.Aref(j+1)*=dag(T2);
		 
   		 	  psi.leftLim(old_llim);
   		  psi.rightLim(old_rlim);
            if(env.rank()+1 == env.nnodes()/2 
            && dir == Fromright 
            && j == jmid)
                {
                printfln("Truncation error for sweep %d (node %d) at site %d is terr=%.12f",
                         sw,b,j,spec.truncerr());
                }

            args.add("AtBond",j);
            args.add("AtBoundary",false);
            args.add("Energy",energy);
            args.add("Direction",dir);
            obs.measure(args);

	     if((psw.ha == 1 && psw.j+numCenter-1 != length(psi) && dir==Fromleft) || (psw.ha == 2 && psw.j != 1 && dir==Fromright))

                {
                auto b1 = (psw.ha == 1 ? psw.j+1 : psw.j);
	
		 phi0 = psi.A(b1);

			 	      PH.numCenter(numCenter-1);
				       PH.position(b1,psi);
 
				       applyExp(PH,phi0,-t/2,args);
				
                 if(DoNormalize)
		   {		       	          phi0 /= norm(phi0);
		   }
              
		        psi.ref(b1) = phi0;
              
		}
	     else if((psw.ha == 1 && psw.j != 1  && dir==Fromright) || (psw.ha == 2 && psw.j+numCenter-1 != length(psi) && dir==Fromleft))
                 {
		   auto b1 = (psw.ha == 1 ?  psw.j : psw.j+1);
	
		    phi0 = psi(b1);
        
          
 
		     PH.numCenter(numCenter-1);
		       PH.position(b1,psi);
 
		           applyExp(PH,phi0,-t/2,args);
		
                 if(DoNormalize)
		   {	    phi0 /= norm(phi0);
		   }
                
      
		        psi.ref(b1) = phi0;
  
         
	       	}
		
            if(psw.atRight() && dir==Fromleft && bool(mboxR))
                {
	
     				auto n = j+1;

  PH.numCenter(numCenter);
        PH.position(n,psi);

                 Boundary<EdgeType> B;
                 mboxR.receive(B);
                 psi.Aref(n+1) = B.A;
                 PH.R(B.HH);
                 B = Boundary<EdgeType>(); //to save memory

                 auto& V = Vs.at(b);

                auto phi = psi.A(n)*V*psi.A(n+1);
	
			applyExp(PH,phi,t,args);
		
if(DoNormalize)
		  {
				phi /= norm(phi);
		  }
  	auto indx1=findIndex(psi.A(n), "Site");
		  auto indx2=findIndex(psi.A(n+1), "Site");
		  
		  int old_llim=psi.leftLim();
		  int old_rlim=psi.rightLim();
   auto [T1, T2]=applyLbo(phi, indx1, indx2,n, args);
   auto id1=findInds(T1, "OM");
   auto id2=findInds(T2, "OM");
   double new_dim(std::max(static_cast<double>(dim(id2)), static_cast<double>(dim(id1))));
    	       lbo_max=std::max(new_dim, lbo_max);
   	      phi*=T1;
   	  phi*=T2;
	  psi.Aref(n)*=T1;
   	  psi.Aref(n+1)*=T2;
	  psi.leftLim(old_llim);
	       	  psi.rightLim(old_rlim);

		  // auto spec = isvd(phi,psi.Aref(n),V,psi.Aref(n+1), args);
		  auto args_cp=args;
		  args_cp.add("MaxDim",maxLinkDim(psi));
 auto spec = isvd(phi,psi.Aref(n),V,psi.Aref(n+1), args_cp);
			
     	   old_llim=psi.leftLim();
     old_rlim=psi.rightLim();
    
     		  psi.Aref(n)*=dag(T1);
     		  psi.Aref(n+1)*=dag(T2);
		 
		  psi.leftLim(old_llim);
     		   psi.rightLim(old_rlim);
		  
                B.HH = PH.L();
                B.UU = psi.A(n)*V;
                B.A = psi.A(n+1);
                B.energy = energy;
                mboxR.send(B);

                psi.Aref(n+1) *= V;
                psi.rightLim(n+1);





	
		phi0 = psi.A(n);
	
			 	      PH.numCenter(numCenter-1);
				       PH.position(n,psi);
 
				       applyExp(PH,phi0,-t/2,args);
				       //				
                if(DoNormalize)
		  {
		       		        	          phi0 /= norm(phi0);
                
		  }
               
		        psi.ref(n) = phi0;
                 }
             else if(psw.atLeft() && dir==Fromright && bool(mboxL))
                 {
		
			    PH.numCenter(numCenter);
		
			    
                auto n = j-1;
	
  PH.numCenter(numCenter);
                PH.position(n,psi);

                Boundary<EdgeType> B;
                B.A = psi.A(n+1);
                B.HH = PH.R();
                mboxL.send(B);

                mboxL.receive(B);
                PH.L(B.HH);
                PH.shift(n,Fromleft,B.UU);
                psi.Aref(n+1) = B.A;
                energy = B.energy;
		phi0 = psi.A(n+1);
      PH.numCenter(numCenter-1);
      PH.position(n+1,psi);
      
      applyExp(PH,phi0,-t/2,args);
       if(DoNormalize)
		  {
	       	          phi0 /= norm(phi0);
		  }
     
		        psi.ref(n+1) = phi0;

                 }
             }

           }



    return lbo_max;

    } // ptdvpWorker

   template < typename HamT, typename step_type>
     std::tuple<double,double,double> 
Act_lboBT_secimpl(Environment const& env,
            Partition const& P,
            MPS & psi,
            MPO & H,
            std::vector<ITensor> & Vs,
            HamT & PH,
            Sweeps const& sweeps,
            Observer & obs,
	    step_type t,
    Args args, bool Testing)
    {

    using EdgeType = stdx::decay_t<decltype(PH.L())>;
    //Maximum number of Davidson iterations at boundaries
    auto boundary_niter = args.getInt("BoundaryIter",sweeps.niter(1)+6);
    bool DoNormalize=args.getBool("DoNormalize",false);
    auto b = env.rank()+1;
    auto jl = P.begin(b);
    auto jr = P.end(b);
    double lbo_max_svd=0;
    double lbo_max_solv=0;
    double lbo_mean_solv=0;
    Real energy = 0.;
  ITensor phi0,phi1;
    auto psw = ParallelSweeper(jl,jr,b);

    MailBox mboxL,mboxR;
    if(not env.firstNode()) mboxL = MailBox(env,env.rank()-1);
    if(not env.lastNode())  mboxR = MailBox(env,env.rank()+1);
    if(Testing)
      {
    psi.leftLim(jl);
    psi.rightLim(jr);
    psi.position(psw.j);
    Testing=false;
      }


    //Include MPINode number to use in the observer
    args.add("BlockStart",jl);
    args.add("BlockEnd",jr);
    args.add("MPINode",env.rank()+1);
    int jmid = (jr-jl)/2.;
    const int numCenter = 2;
    //args.getInt("NumCenter",2);
     if(numCenter != 2)
       {std::cout<< "does not work for tdvp1 yet!!"<<std::endl;
	 return std::tuple<double,double,double>(0,0,0);}

    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
        {
        args.add("Sweep",sw);
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("Noise",sweeps.noise(sw));
        args.add("MaxIter",sweeps.niter(sw));

         for(psw.newSweep(); psw.doingFull(); ++psw)
             {
            auto j = psw.j;
            auto dir = psw.dir();

	        PH.numCenter(numCenter);
            PH.position(j,psi);

             phi1 = psi.A(j)*psi.A(j+1);

auto indx1=findIndex(psi.A(j), "Site");
 auto indx2=findIndex(psi.A(j+1), "Site");
	     //	      applyExp(PH,phi1,t/2,args);
 // auto  [lbomax_it,lbo_mean]=applyExp_lbo(PH, phi1, t/2, indx1, indx2, j, H,args);

auto [max_solv,mean_solv]= Runge_Kutta_4_lbo(PH,phi1,t/2,indx1, indx2, j,H,args);;		
lbo_max_solv=std::max(max_solv, lbo_max_solv);
 lbo_mean_solv=std::max(mean_solv, lbo_mean_solv);			
		  	         int old_llim=psi.leftLim();
		  int old_rlim=psi.rightLim();
	      auto [T1, T2]=applyLbo(phi1, indx1, indx2,j, args);
	       auto id1=findInds(T1, "OM");
	       auto id2=findInds(T2, "OM");
	       double new_dim(std::max(static_cast<double>(dim(id2)), static_cast<double>(dim(id1))));
	       lbo_max_svd=std::max(new_dim, lbo_max_svd);

	  phi1*=T1;
	  phi1*=T2;
	  psi.Aref(j)*=T1;
	  
	  psi.Aref(j+1)*=T2;

	  	  psi.leftLim(old_llim);
	  		  psi.rightLim(old_rlim);
		 if(DoNormalize)
	      
		   {
phi1 /= norm(phi1);
		   }
	     auto energy=0;
            auto spec = psi.svdBond(j,phi1,dir,PH,args);
   old_llim=psi.leftLim();
    old_rlim=psi.rightLim();
            	  psi.Aref(j)*=dag(T1);
   		 psi.Aref(j+1)*=dag(T2);
		 
   		 	  psi.leftLim(old_llim);
   		  psi.rightLim(old_rlim);
            if(env.rank()+1 == env.nnodes()/2 
            && dir == Fromright 
            && j == jmid)
                {
                printfln("Truncation error for sweep %d (node %d) at site %d is terr=%.12f",
                         sw,b,j,spec.truncerr());
                }

            args.add("AtBond",j);
            args.add("AtBoundary",false);
            args.add("Energy",energy);
            args.add("Direction",dir);
            obs.measure(args);

	     if((psw.ha == 1 && psw.j+numCenter-1 != length(psi) && dir==Fromleft) || (psw.ha == 2 && psw.j != 1 && dir==Fromright))

                {
                auto b1 = (psw.ha == 1 ? psw.j+1 : psw.j);
	
		 phi0 = psi.A(b1);

			 	      PH.numCenter(numCenter-1);
				       PH.position(b1,psi);
 
				       //				       applyExp(PH,phi0,-t/2,args);
		                Runge_Kutta_4(PH,phi0,-t/2,args);
				
                 if(DoNormalize)
		   {		       	          phi0 /= norm(phi0);
		   }
              
		        psi.ref(b1) = phi0;
              
		}
	     else if((psw.ha == 1 && psw.j != 1  && dir==Fromright) || (psw.ha == 2 && psw.j+numCenter-1 != length(psi) && dir==Fromleft))
                 {
		   auto b1 = (psw.ha == 1 ?  psw.j : psw.j+1);
	
		    phi0 = psi(b1);
        
          
 
		     PH.numCenter(numCenter-1);
		       PH.position(b1,psi);
 
		       //		           applyExp(PH,phi0,-t/2,args);
		       Runge_Kutta_4(PH,phi0,-t/2,args);
		
                 if(DoNormalize)
		   {	    phi0 /= norm(phi0);
		   }
                
      
		        psi.ref(b1) = phi0;
  
         
	       	}
		
            if(psw.atRight() && dir==Fromleft && bool(mboxR))
                {

     				auto n = j+1;

  PH.numCenter(numCenter);
  PH.position(n,psi);

                 Boundary<EdgeType> B;
                 mboxR.receive(B);
                 psi.Aref(n+1) = B.A;
                 PH.R(B.HH);
                 B = Boundary<EdgeType>(); //to save memory
auto indx1=findIndex(psi.A(n), "Site");
 auto indx2=findIndex(psi.A(n+1), "Site");
                 auto& V = Vs.at(b);

                auto phi = psi.A(n)*V*psi.A(n+1);



		//					applyExp(PH,phi,t,args);

		//		auto  [lbomax_it,lbo_mean]=applyExp_lbo(PH, phi, t, indx1, indx2, n, H,args)	    

auto [max_solv,mean_solv]= Runge_Kutta_4_lbo(PH,phi,t,indx1, indx2, n,H,args);;		
lbo_max_solv=std::max(max_solv, lbo_max_solv);
 lbo_mean_solv=std::max(mean_solv, lbo_mean_solv);			
if(DoNormalize)
		  {
				phi /= norm(phi);
		  }
  	
 
 int old_llim=psi.leftLim();
 int old_rlim=psi.rightLim();
   auto [T1, T2]=applyLbo(phi, indx1, indx2,n, args);
   auto id1=findInds(T1, "OM");
   auto id2=findInds(T2, "OM");
   double new_dim(std::max(static_cast<double>(dim(id2)), static_cast<double>(dim(id1))));
    	       lbo_max_svd=std::max(new_dim, lbo_max_svd);
   	      phi*=T1;
   	  phi*=T2;
   	  psi.Aref(n)*=T1;
   	  psi.Aref(n+1)*=T2;
	  psi.leftLim(old_llim);
	       	  psi.rightLim(old_rlim);

		  
		  auto args_cp=args;
		  		  args_cp.add("MaxDim",maxLinkDim(psi));

 auto spec = isvd(phi,psi.Aref(n),V,psi.Aref(n+1), args_cp);
			
     	   old_llim=psi.leftLim();
     old_rlim=psi.rightLim();
    
     psi.Aref(n)*=dag(T1);
     psi.Aref(n+1)*=dag(T2);
		 
     psi.leftLim(old_llim);
     psi.rightLim(old_rlim);
		  
     B.HH = PH.L();
     B.UU = psi.A(n)*V;
     B.A = psi.A(n+1);
     B.energy = energy;
     mboxR.send(B);

     psi.Aref(n+1) *= V;
     psi.rightLim(n+1);





	
		phi0 = psi.A(n);
	
			 	      PH.numCenter(numCenter-1);
				       PH.position(n,psi);
 		                Runge_Kutta_4(PH,phi0,-t/2,args);
				       //			       applyExp(PH,phi0,-t/2,args);
				       //				
                if(DoNormalize)
		  {
		       		        	          phi0 /= norm(phi0);
                
		  }
               
		        psi.ref(n) = phi0;
                 }
             else if(psw.atLeft() && dir==Fromright && bool(mboxL))
                 {
		
			    PH.numCenter(numCenter);
		
			    
                auto n = j-1;
	
  PH.numCenter(numCenter);
                PH.position(n,psi);

                Boundary<EdgeType> B;
                B.A = psi.A(n+1);
                B.HH = PH.R();
                mboxL.send(B);

                mboxL.receive(B);
                PH.L(B.HH);
                PH.shift(n,Fromleft,B.UU);
                psi.Aref(n+1) = B.A;
                energy = B.energy;
		phi0 = psi.A(n+1);
      PH.numCenter(numCenter-1);
      PH.position(n+1,psi);
      
      //applyExp(PH,phi0,-t/2,args);
      Runge_Kutta_4(PH,phi0,-t/2,args);
       if(DoNormalize)
		  {
	       	          phi0 /= norm(phi0);
		  }
     
		        psi.ref(n+1) = phi0;

                 }
             }

           }



    return std::tuple<double,double,double>(lbo_max_svd,lbo_max_solv,lbo_mean_solv);

    } 

   template < typename HamT, typename step_type>
     std::tuple<double,double,double> 
Act_lbo_secimpl(Environment const& env,
            Partition const& P,
            MPS & psi,
            MPO & H,
            std::vector<ITensor> & Vs,
            HamT & PH,
            Sweeps const& sweeps,
            Observer & obs,
	    step_type t,
    Args args, bool Testing)
    {

    using EdgeType = stdx::decay_t<decltype(PH.L())>;
    //Maximum number of Davidson iterations at boundaries
    auto boundary_niter = args.getInt("BoundaryIter",sweeps.niter(1)+6);
    bool DoNormalize=args.getBool("DoNormalize",false);
    auto b = env.rank()+1;
    auto jl = P.begin(b);
    auto jr = P.end(b);
    double lbo_max_svd=0;
    double lbo_max_solv=0;
    double lbo_mean_solv=0;
    Real energy = 0.;
  ITensor phi0,phi1;
    auto psw = ParallelSweeper(jl,jr,b);

    MailBox mboxL,mboxR;
    if(not env.firstNode()) mboxL = MailBox(env,env.rank()-1);
    if(not env.lastNode())  mboxR = MailBox(env,env.rank()+1);
    if(Testing)
      {
    psi.leftLim(jl);
    psi.rightLim(jr);
    psi.position(psw.j);
    Testing=false;
      }


    //Include MPINode number to use in the observer
    args.add("BlockStart",jl);
    args.add("BlockEnd",jr);
    args.add("MPINode",env.rank()+1);
    int jmid = (jr-jl)/2.;
    const int numCenter = 2;
    //args.getInt("NumCenter",2);
     if(numCenter != 2)
       {std::cout<< "does not work for tdvp1 yet!!"<<std::endl;
	 return std::tuple<double,double,double>(0,0,0);}

    for(int sw = 1; sw <= sweeps.nsweep(); ++sw)
        {
        args.add("Sweep",sw);
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("Noise",sweeps.noise(sw));
        args.add("MaxIter",sweeps.niter(sw));

         for(psw.newSweep(); psw.doingFull(); ++psw)
             {
            auto j = psw.j;
            auto dir = psw.dir();

	        PH.numCenter(numCenter);
            PH.position(j,psi);

             phi1 = psi.A(j)*psi.A(j+1);

auto indx1=findIndex(psi.A(j), "Site");
 auto indx2=findIndex(psi.A(j+1), "Site");
	     //	      applyExp(PH,phi1,t/2,args);
 // auto  [lbomax_it,lbo_mean]=applyExp_lbo(PH, phi1, t/2, indx1, indx2, j, H,args);

auto [max_solv,mean_solv]= Runge_Kutta_4_lbo(PH,phi1,t/2,indx1, indx2, j,H,args);;		
lbo_max_solv=std::max(max_solv, lbo_max_solv);
 lbo_mean_solv=std::max(mean_solv, lbo_mean_solv);			
		  	         int old_llim=psi.leftLim();
		  int old_rlim=psi.rightLim();
	      auto [T1, T2]=applyLbo(phi1, indx1, indx2,j, args);
	       auto id1=findInds(T1, "OM");
	       auto id2=findInds(T2, "OM");
	       double new_dim(std::max(static_cast<double>(dim(id2)), static_cast<double>(dim(id1))));
	       lbo_max_svd=std::max(new_dim, lbo_max_svd);

	  phi1*=T1;
	  phi1*=T2;
	  psi.Aref(j)*=T1;
	  
	  psi.Aref(j+1)*=T2;

	  	  psi.leftLim(old_llim);
	  		  psi.rightLim(old_rlim);
		 if(DoNormalize)
	      
		   {
phi1 /= norm(phi1);
		   }
	     auto energy=0;
            auto spec = psi.svdBond(j,phi1,dir,PH,args);
   old_llim=psi.leftLim();
    old_rlim=psi.rightLim();
            	  psi.Aref(j)*=dag(T1);
   		 psi.Aref(j+1)*=dag(T2);
		 
   		 	  psi.leftLim(old_llim);
   		  psi.rightLim(old_rlim);
            if(env.rank()+1 == env.nnodes()/2 
            && dir == Fromright 
            && j == jmid)
                {
                printfln("Truncation error for sweep %d (node %d) at site %d is terr=%.12f",
                         sw,b,j,spec.truncerr());
                }

            args.add("AtBond",j);
            args.add("AtBoundary",false);
            args.add("Energy",energy);
            args.add("Direction",dir);
            obs.measure(args);

	     if((psw.ha == 1 && psw.j+numCenter-1 != length(psi) && dir==Fromleft) || (psw.ha == 2 && psw.j != 1 && dir==Fromright))

                {
                auto b1 = (psw.ha == 1 ? psw.j+1 : psw.j);
	
		 phi0 = psi.A(b1);

			 	      PH.numCenter(numCenter-1);
				       PH.position(b1,psi);
 
				       //				       applyExp(PH,phi0,-t/2,args);
		                Runge_Kutta_4(PH,phi0,-t/2,args);
				
                 if(DoNormalize)
		   {		       	          phi0 /= norm(phi0);
		   }
              
		        psi.ref(b1) = phi0;
              
		}
	     else if((psw.ha == 1 && psw.j != 1  && dir==Fromright) || (psw.ha == 2 && psw.j+numCenter-1 != length(psi) && dir==Fromleft))
                 {
		   auto b1 = (psw.ha == 1 ?  psw.j : psw.j+1);
	
		    phi0 = psi(b1);
        
          
 
		     PH.numCenter(numCenter-1);
		       PH.position(b1,psi);
 
		       //		           applyExp(PH,phi0,-t/2,args);
		       Runge_Kutta_4(PH,phi0,-t/2,args);
		
                 if(DoNormalize)
		   {	    phi0 /= norm(phi0);
		   }
                
      
		        psi.ref(b1) = phi0;
  
         
	       	}
		
            if(psw.atRight() && dir==Fromleft && bool(mboxR))
                {

     				auto n = j+1;

  PH.numCenter(numCenter);
  PH.position(n,psi);

                 Boundary<EdgeType> B;
                 mboxR.receive(B);
                 psi.Aref(n+1) = B.A;
                 PH.R(B.HH);
                 B = Boundary<EdgeType>(); //to save memory
		 ////////////////////////////////////////////////////////
auto indx1=findIndex(psi.A(n), "Site");
 auto indx2=findIndex(psi.A(n+1), "Site");
                 auto& V = Vs.at(b);

                auto phi = psi.A(n)*V*psi.A(n+1);


auto [max_solv,mean_solv]= Runge_Kutta_4_lbo(PH,phi,t/2,indx1, indx2, n,H,args);;		
lbo_max_solv=std::max(max_solv, lbo_max_solv);
 lbo_mean_solv=std::max(mean_solv, lbo_mean_solv);			
if(DoNormalize)
		  {
				phi /= norm(phi);
		  }
  	
 
 int old_llim=psi.leftLim();
 int old_rlim=psi.rightLim();
   auto [T1, T2]=applyLbo(phi, indx1, indx2,n, args);
   auto id1=findInds(T1, "OM");
   auto id2=findInds(T2, "OM");
   double new_dim(std::max(static_cast<double>(dim(id2)), static_cast<double>(dim(id1))));
    	       lbo_max_svd=std::max(new_dim, lbo_max_svd);
   	      phi*=T1;
   	  phi*=T2;
   	  psi.Aref(n)*=T1;
   	  psi.Aref(n+1)*=T2;
	  psi.leftLim(old_llim);
	       	  psi.rightLim(old_rlim);

		  



 auto spec = isvd(phi,psi.Aref(n),V,psi.Aref(n+1), args);
			
     	   old_llim=psi.leftLim();
     old_rlim=psi.rightLim();
    
     psi.Aref(n)*=dag(T1);
     psi.Aref(n+1)*=dag(T2);
		 
     psi.leftLim(old_llim);
     psi.rightLim(old_rlim);
     //////////////////////////////////////////////////////////////////////////////////////////
     // ned to time evolve tensor twice in two steps
     //////////////////////////////////////////////////////////////////////////////////
 indx1=findIndex(psi.A(n), "Site");
  indx2=findIndex(psi.A(n+1), "Site");
                  V = Vs.at(b);

                 phi = psi.A(n)*V*psi.A(n+1);


auto [max_solv_2,mean_solv_2]= Runge_Kutta_4_lbo(PH,phi,t/2,indx1, indx2, n,H,args);;		
lbo_max_solv=std::max(max_solv_2, lbo_max_solv);
 lbo_mean_solv=std::max(mean_solv_2, lbo_mean_solv);			
if(DoNormalize)
		  {
				phi /= norm(phi);
		  }
  	
 
  old_llim=psi.leftLim();
  old_rlim=psi.rightLim();
   auto [T1_2, T2_2]=applyLbo(phi, indx1, indx2,n, args);
    id1=findInds(T1_2, "OM");
    id2=findInds(T2_2, "OM");
   double new_dim_2(std::max(static_cast<double>(dim(id2)), static_cast<double>(dim(id1))));
    	       lbo_max_svd=std::max(new_dim_2, lbo_max_svd);
   	      phi*=T1_2;
   	  phi*=T2_2;
   	  psi.Aref(n)*=T1_2;
   	  psi.Aref(n+1)*=T2_2;
	  psi.leftLim(old_llim);
	       	  psi.rightLim(old_rlim);

		  



  spec = isvd(phi,psi.Aref(n),V,psi.Aref(n+1), args);
			
     	   old_llim=psi.leftLim();
     old_rlim=psi.rightLim();
    
     psi.Aref(n)*=dag(T1_2);
     psi.Aref(n+1)*=dag(T2_2);
		 
     psi.leftLim(old_llim);
     psi.rightLim(old_rlim);
     /////////////////////////////////////////////////////////////////////////////////////////
     B.HH = PH.L();
     B.UU = psi.A(n)*V;
     B.A = psi.A(n+1);
     B.energy = energy;
     mboxR.send(B);

     psi.Aref(n+1) *= V;
     psi.rightLim(n+1);





	
		phi0 = psi.A(n);
	
			 	      PH.numCenter(numCenter-1);
				       PH.position(n,psi);
 		                Runge_Kutta_4(PH,phi0,-t/2,args);
				       //			       applyExp(PH,phi0,-t/2,args);
				       //				
                if(DoNormalize)
		  {
		       		        	          phi0 /= norm(phi0);
                
		  }
               
		        psi.ref(n) = phi0;
                 }
             else if(psw.atLeft() && dir==Fromright && bool(mboxL))
                 {
		
			    PH.numCenter(numCenter);
		
			    
                auto n = j-1;
	
  PH.numCenter(numCenter);
                PH.position(n,psi);

                Boundary<EdgeType> B;
                B.A = psi.A(n+1);
                B.HH = PH.R();
                mboxL.send(B);

                mboxL.receive(B);
                PH.L(B.HH);
                PH.shift(n,Fromleft,B.UU);
                psi.Aref(n+1) = B.A;
                energy = B.energy;
		phi0 = psi.A(n+1);
      PH.numCenter(numCenter-1);
      PH.position(n+1,psi);
      
      //applyExp(PH,phi0,-t/2,args);
      Runge_Kutta_4(PH,phi0,-t/2,args);
       if(DoNormalize)
		  {
	       	          phi0 /= norm(phi0);
		  }
     
		        psi.ref(n+1) = phi0;

                 }
             }

           }



    return std::tuple<double,double,double>(lbo_max_svd,lbo_max_solv,lbo_mean_solv);

    } // ptdvpWorker

} //namespace itensor

#endif
