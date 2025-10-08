
#ifndef LBO_H
#define LBO_H
#include<tuple>
//#include"diag.hpp"
#include "itensor/all.h"
//#include"herm.hpp"
namespace itensor{
  //  const int LboMaxQ=4;
  //const int LboMin=1;
  //const double errLbo=1E-17;
//   ITensor applyLbo_ss2(ITensor& phi, Index i1, size_t n, const Args& args)
// {
//   auto lboCut = args.getReal("CutoffLBO",0.);
 
//   auto lboMax = args.getInt("MaxDimLBO",100);
//   auto lboMin = args.getInt("MinDimLBO",20);

//  auto tags = TagSet("Site, OM");
//  auto inds_site=findInds(phi, "Link");
//      //     std::cout<< "lbo max allowed"<< lboMax<<std::endl; 
//      tags.addTags("n="+str(n));
//      auto [U,S,V_1] = svd(phi,inds_site, {"RightTags", tags,"Cutoff", lboCut,  "MaxDim=", lboMax, "MinDim=",10,"RespectDegenerate" , true,"ShowEigs", true, "Truncate", false});

// 	  auto id=findInds(V_1, "OM");
// 	  if(dim(id)<10)
// 	    {std::cout<< "error "<<std::endl;
// 	      print(V_1);
// 	      }

// 	  return dag(V_1);
// }

  ITensor applyLbo_ss(ITensor phi, Index i1, size_t n, const Args& args)
{
  auto lboCut = args.getReal("CutoffLBO",0.);
  auto lboMax = args.getInt("MaxDimLBO",100);
  auto lboMin = args.getInt("MinDimLBO",20);
  
   auto phidag=dag(phi);
 phidag.prime(i1);
 auto phidag2=dag(phi*phidag);

 auto tags = TagSet("Site, OM");

     tags.addTags("n="+str(n));
     auto [U1, D1] =diagPosSemiDef(phidag2,{"Tags", tags,"Cutoff", lboCut,  "MaxDim=", lboMax, "MinDim=",lboMin,"RespectDegenerate" ,true, "Truncate", true,"ShowEigs", false});
     

  return U1;
}

  std::tuple<ITensor, ITensor> applyLbo(ITensor& phi, Index i1, Index i2, int n, const Args& args)
{
  auto lboCut = args.getReal("CutoffLBO",0.);
 
  auto lboMax = args.getInt("MaxDimLBO",100);
  auto lboMin = args.getInt("MinDimLBO",0);

   auto phidag=dag(phi);

 phidag.prime(i1);

auto phidag2=phi*phidag;

     auto tags = TagSet("Site, OM");
    
     tags.addTags("n="+str(n));
    

     auto [U1, D1] =diagPosSemiDef(phidag2,{"Tags", tags,"Cutoff", lboCut,  "MaxDim=", lboMax, "MinDim=",lboMin,"RespectDegenerate" , true,"Truncate", true,"ShowEigs",false});
 

   auto phidag12=dag(phi);
 phidag12.prime(i2);
 auto phidag22=phi*phidag12;
       tags = TagSet("Site, OM");

     tags.addTags("n="+str(n+1));

     auto [U2, D2] =diagPosSemiDef(phidag22,{ "Tags", tags,"Cutoff", lboCut,"MaxDim=", lboMax,"MinDim=",lboMin,"RespectDegenerate" , true,"Truncate", true, "ShowEigs", false});


 
  std::tuple<ITensor, ITensor> trafos(U1, U2);
  return trafos;
}
}
#endif /* LBO_H */
