#ifndef __ITENSOR_WRITE_DATA
#define __ITENSOR_WRITE_DATA
#include"files.hpp"
#include"vector"

namespace itensor {
void write_data_GS(MPS& psi1x, MPS& psi2x, double& final_lboMaxB, double& final_lboMaxK, double time, MPO Op, Args argsObs, std::string DIR, std::string obsname,Args const& args = Args::global())
{
  std::vector<std::complex<double>> o_vec;
  std::vector<double> t_vec;
  std::vector<double> normB_vec;
  std::vector<double> normK_vec;
  std::vector<double> avBdB_vec;
  std::vector<double> maxBdB_vec;
  std::vector<double> avBdK_vec;
  std::vector<double> maxBdK_vec;
  std::vector<double> maxlboK_vec;
  std::vector<double> maxlboB_vec;

  int N=length(psi1x);
  psi1x.position(int(N/2), args);
  psi2x.position(int(N/2), args);
	 psi2x= applyMPO(Op,psi2x,argsObs); 
	 psi2x.noPrime();
	 std::complex<double> Oval= innerC(psi2x,psi1x);
			  
	 auto NormB=norm(psi2x);
	 auto NormK=norm(psi1x);
	 auto maxBdK=maxLinkDim(psi1x);
	 auto maxBdB=maxLinkDim(psi2x);
	 auto avBdK=averageLinkDim(psi1x);
	 auto avBdB=averageLinkDim(psi2x);
	
 
    std::cout<< time<<" overlap" << Oval<< " mx psi1  "<< maxBdK<<" av psi1  "<< avBdK<<" mx psi2  "<< maxBdB<<" av psi 2  "<< avBdB<<" Norm 1.1= "<< NormK<<" norm 2.1= "<< NormB<<"\n";

    t_vec.push_back(time);
    o_vec.push_back(Oval);
    normK_vec.push_back(NormK);
    normB_vec.push_back(NormB);
 
    avBdK_vec.push_back(avBdK);
    maxBdK_vec.push_back(maxBdK);
    avBdB_vec.push_back(avBdB);
    maxBdB_vec.push_back(maxBdB);
 

    maxlboB_vec.push_back(final_lboMaxB);
    maxlboK_vec.push_back(final_lboMaxK);
    auto vecSize=t_vec.size();
    itensor::ToFile(t_vec, DIR+"/time.dat", vecSize);
    itensor::ToFile(o_vec, DIR+"/"+obsname+".dat",vecSize);
    itensor::ToFile(avBdB_vec, DIR+"/avBdB.dat", vecSize);
    itensor::ToFile(maxBdB_vec, DIR+"/maxBdB.dat", vecSize);
    itensor::ToFile(avBdK_vec, DIR+"/avBdK.dat", vecSize);
    itensor::ToFile(maxBdK_vec, DIR+"/maxBdK.dat", vecSize);
    itensor::ToFile(maxlboK_vec, DIR+"/lbo_svdmxK.dat", vecSize);
    itensor::ToFile(maxlboB_vec, DIR+"/lbo_svdmxB.dat", vecSize);
    itensor::ToFile(normK_vec, DIR+"/normK.dat", vecSize);
    itensor::ToFile(normB_vec, DIR+"/normB.dat", vecSize);
      return;
}
 void write_data_GS_eff(MPS& psi1x, MPS& psi2x, double& final_lboMaxB, double& final_lboMaxK, double time,  Args argsObs, std::string DIR, std::string obsname,double E_0,Args const& args = Args::global())
{
  std::vector<std::complex<double>> o_vec;
  std::vector<double> t_vec;
  std::vector<double> normB_vec;
  std::vector<double> normK_vec;
  std::vector<double> avBdB_vec;
  std::vector<double> maxBdB_vec;
  std::vector<double> avBdK_vec;
  std::vector<double> maxBdK_vec;
  std::vector<double> maxlboK_vec;
  std::vector<double> maxlboB_vec;


  psi1x.position(1, args);
  psi2x.position(1, args);

  std::complex<double> Oval= std::exp(Cplx_i*E_0)*innerC(psi2x,psi1x);
			  
	 auto NormB=norm(psi2x);
	 auto NormK=norm(psi1x);
	 auto maxBdK=maxLinkDim(psi1x);
	 auto maxBdB=maxLinkDim(psi2x);
	 auto avBdK=averageLinkDim(psi1x);
	 auto avBdB=averageLinkDim(psi2x);
	
 
    std::cout<< time<<" overlap" << Oval<< " mx psi1  "<< maxBdK<<" av psi1  "<< avBdK<<" mx psi2  "<< maxBdB<<" av psi 2  "<< avBdB<<" Norm 1.1= "<< NormK<<" norm 2.1= "<< NormB<<"\n";

    t_vec.push_back(time);
    o_vec.push_back(Oval);
    normK_vec.push_back(NormK);
    normB_vec.push_back(NormB);
 
    avBdK_vec.push_back(avBdK);
    maxBdK_vec.push_back(maxBdK);
    avBdB_vec.push_back(avBdB);
    maxBdB_vec.push_back(maxBdB);
 

    maxlboB_vec.push_back(final_lboMaxB);
    maxlboK_vec.push_back(final_lboMaxK);
    auto vecSize=t_vec.size();
    itensor::ToFile(t_vec, DIR+"/time.dat", vecSize);
    itensor::ToFile(o_vec, DIR+"/"+obsname+".dat",vecSize);
    itensor::ToFile(avBdB_vec, DIR+"/avBdB.dat", vecSize);
    itensor::ToFile(maxBdB_vec, DIR+"/maxBdB.dat", vecSize);
    itensor::ToFile(avBdK_vec, DIR+"/avBdK.dat", vecSize);
    itensor::ToFile(maxBdK_vec, DIR+"/maxBdK.dat", vecSize);
    itensor::ToFile(maxlboK_vec, DIR+"/lbo_svdmxK.dat", vecSize);
    itensor::ToFile(maxlboB_vec, DIR+"/lbo_svdmxB.dat", vecSize);
    itensor::ToFile(normK_vec, DIR+"/normK.dat", vecSize);
    itensor::ToFile(normB_vec, DIR+"/normB.dat", vecSize);
      return;

}
void write_data_FT(MPS& psi1x, MPS& psi2x, double& final_lboMaxB, double& final_lboMaxK, double time,  Args argsObs, std::string DIR, std::string obsname,Args const& args = Args::global())
{
  std::vector<std::complex<double>> o_vec;
  std::vector<double> t_vec;
  std::vector<double> normB_vec;
  std::vector<double> normK_vec;
  std::vector<double> avBdB_vec;
  std::vector<double> maxBdB_vec;
  std::vector<double> avBdK_vec;
  std::vector<double> maxBdK_vec;
  std::vector<double> maxlboK_vec;
  std::vector<double> maxlboB_vec;


    psi1x.position(1, args);
   psi2x.position(1, args);

	 std::complex<double> Oval= innerC(psi2x,psi1x);
			  
	 auto NormB=norm(psi2x);
	 auto NormK=norm(psi1x);
	 auto maxBdK=maxLinkDim(psi1x);
	 auto maxBdB=maxLinkDim(psi2x);
	 auto avBdK=averageLinkDim(psi1x);
	 auto avBdB=averageLinkDim(psi2x);
	
 
    std::cout<< time<<" overlap" << Oval<< " mx psi1  "<< maxBdK<<" av psi1  "<< avBdK<<" mx psi2  "<< maxBdB<<" av psi 2  "<< avBdB<<" Norm 1.1= "<< NormK<<" norm 2.1= "<< NormB<<"\n";

    t_vec.push_back(time);
    o_vec.push_back(Oval);
    normK_vec.push_back(NormK);
    normB_vec.push_back(NormB);
 
    avBdK_vec.push_back(avBdK);
    maxBdK_vec.push_back(maxBdK);
    avBdB_vec.push_back(avBdB);
    maxBdB_vec.push_back(maxBdB);
 

    maxlboB_vec.push_back(final_lboMaxB);
    maxlboK_vec.push_back(final_lboMaxK);
    auto vecSize=t_vec.size();
    itensor::ToFile(t_vec, DIR+"/time.dat", vecSize);
    itensor::ToFile(o_vec, DIR+"/"+obsname+".dat",vecSize);
    itensor::ToFile(avBdB_vec, DIR+"/avBdB.dat", vecSize);
    itensor::ToFile(maxBdB_vec, DIR+"/maxBdB.dat", vecSize);
    itensor::ToFile(avBdK_vec, DIR+"/avBdK.dat", vecSize);
    itensor::ToFile(maxBdK_vec, DIR+"/maxBdK.dat", vecSize);
    itensor::ToFile(maxlboK_vec, DIR+"/lbo_svdmxK.dat", vecSize);
    itensor::ToFile(maxlboB_vec, DIR+"/lbo_svdmxB.dat", vecSize);
    itensor::ToFile(normK_vec, DIR+"/normK.dat", vecSize);
    itensor::ToFile(normB_vec, DIR+"/normB.dat", vecSize);
      return;

}
void write_state(MPS& psi1x, MPS& psi2x, std::string& mpsNameK, std::string& mpsNameB, std::string& DIR, std::string stot,bool apply_op_B,bool apply_op_K){
		       std::string sub="MPS";
		       // find right point to inset new file name

		       std::cout<< "data 1 y "<<maxLinkDim(psi1x)<<std::endl;
		       std::cout<< "data 1 y "<<maxLinkDim(psi2x)<<std::endl;
		       std::vector<size_t> position_1;
		       std::vector<size_t> position_2;
		       size_t pos_1=mpsNameK.find(sub,0);
		       size_t pos_2=mpsNameB.find(sub,0);
		       while(pos_1 !=std::string::npos)
			 {
			   position_1.push_back(pos_1);
			   pos_1=mpsNameK.find(sub, pos_1+1);
}

while(pos_2 !=std::string::npos)
			 {
			   position_2.push_back(pos_2);
			   pos_2=mpsNameB.find(sub, pos_2+1);
}
 size_t find1=position_1.back();
 size_t find2=position_2.back();

    		       std::string app_J_B="";
    		       std::string app_J_K="";

    		       if(apply_op_B)
    			 {app_J_B="appliedJ";}
if(apply_op_K)
    			 {app_J_K="appliedJ";}
    		       std::vector<double> dims_K;
    		       std::vector<double> dims_B;
		       //		       print(psi1x);
    		       for(int j=1; j<length(psi1x); j++)
    			 {
    			   auto ind_K=rightLinkIndex(psi1x,j);
    			   auto ind_B=rightLinkIndex(psi2x,j);

    			   dims_K.push_back(static_cast<double>(dim(ind_K)));
    			   dims_B.push_back(static_cast<double>(dim(ind_B)));
    			 }
    auto vecSize=dims_K.size();
    std::cout<<"here" <<std::endl;
    itensor::ToFile(dims_K, DIR+"/bdsKat"+stot+".dat", vecSize);
    itensor::ToFile(dims_B, DIR+"/bdsBat"+stot+".dat", vecSize);
    itensor::writeToFile(DIR+"/MPSKat"+stot+app_J_K+mpsNameK.substr(find1, mpsNameK.size()),psi1x);
    itensor::writeToFile(DIR+"/MPSBat"+stot+app_J_B+mpsNameB.substr(find2, mpsNameB.size()),psi2x);
    		     
  return;
}
 void write_state_bds(MPS& psi1x, MPS& psi2x,  std::string stot,std::string& DIR){

    		       std::vector<double> dims_K;
    		       std::vector<double> dims_B;
		       //		       print(psi1x);
    		       for(int j=1; j<length(psi1x); j++)
    			 {
    			   auto ind_K=rightLinkIndex(psi1x,j);
    			   auto ind_B=rightLinkIndex(psi2x,j);

    			   dims_K.push_back(static_cast<double>(dim(ind_K)));
    			   dims_B.push_back(static_cast<double>(dim(ind_B)));
    			 }
    auto vecSize=dims_K.size();

    itensor::ToFile(dims_K, DIR+"/bdsKat"+stot+".dat", vecSize);
    itensor::ToFile(dims_B, DIR+"/bdsBat"+stot+".dat", vecSize);


  return;
}
}
#endif
