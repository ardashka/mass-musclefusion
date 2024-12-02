#include "header.hpp"
#include "models/drylyotropicF.hpp"
#include "error_msg.hpp"
#include "random.hpp"
#include "lb.hpp"
#include "tools.hpp"
#include <string> 
#include <fstream>
#include <algorithm>
#include <nlohmann/json.hpp>

using namespace std;
using json = nlohmann::json;
namespace opt = boost::program_options;

// from main.cpp:
extern unsigned nthreads, nsubsteps;
extern double time_step;

DryLyotropicF::DryLyotropicF(unsigned LX, unsigned LY, unsigned BC)
  : Model(LX, LY, BC, BC==0 ? GridType::Periodic : GridType::Layer)
{}
DryLyotropicF::DryLyotropicF(unsigned LX, unsigned LY, unsigned BC, GridType Type)
  : Model(LX, LY, BC, Type)
{}

void DryLyotropicF::Initialize()
{
  // initialize variables
  angle = angle_deg*M_PI/180.;

  // allocate memory
  QQxx.SetSize(LX, LY, Type);
  QQyx.SetSize(LX, LY, Type);
  fQQxx.SetSize(LX, LY, Type);
  fQQyx.SetSize(LX, LY, Type);
  QNxx.SetSize(LX, LY, Type);
  QNyx.SetSize(LX, LY, Type);
  fQNxx.SetSize(LX, LY, Type);
  fQNyx.SetSize(LX, LY, Type);
  Px.SetSize(LX, LY, Type);
  Py.SetSize(LX, LY, Type);
  phi.SetSize(LX, LY, Type);
  phi_tmp.SetSize(LX, LY, Type);
  phn.SetSize(LX, LY, Type);
  ux.SetSize(LX, LY, Type);
  uy.SetSize(LX, LY, Type);
  ux_phi.SetSize(LX, LY, Type);
  uy_phi.SetSize(LX, LY, Type);
  HHxx.SetSize(LX, LY, Type);
  HHyx.SetSize(LX, LY, Type);
  MU.SetSize(LX, LY, Type);
  dxQQxx.SetSize(LX, LY, Type);
  dyQQxx.SetSize(LX, LY, Type);
  dxQQyx.SetSize(LX, LY, Type);
  dyQQyx.SetSize(LX, LY, Type);
  del2QQxx.SetSize(LX, LY, Type);
  del2QQyx.SetSize(LX, LY, Type);
  sigmaXX.SetSize(LX, LY, Type);
  sigmaYY.SetSize(LX, LY, Type);
  sigmaYX.SetSize(LX, LY, Type);
  sigmaXY.SetSize(LX, LY, Type);
  outS.SetSize(LX, LY, Type);
  gradS.SetSize(LX, LY, Type);
  maxstress.SetSize(LX, LY, Type); 
  L1v.SetSize(LX, LY, Type);
  L2v.SetSize(LX, LY, Type);
  L3v.SetSize(LX, LY, Type);
  fusionSum.SetSize(LX, LY, Type);
  Mval.SetSize(LX, LY, Type);


  fusion.SetSize(LX, LY, Type);
  fusn.SetSize(LX, LY, Type);
  fusion_tmp.SetSize(LX, LY, Type);

  FFx.SetSize(LX, LY, Type);
  FFy.SetSize(LX, LY, Type);

  if(nsubsteps>1)
    throw error_msg("time stepping not implemented for this model"
                    ", please set nsubsteps=1.");
}

void DryLyotropicF::ConfigureAtNode(unsigned k)
{
  double nematicOrder = 0;
  double theta, ftheta;
  const unsigned x = GetXPosition(k);
  const unsigned y = GetYPosition(k);
  if(init_config=="circle")
  {
    if(pow(diff(x, LX/2), 2) + pow(diff(y, LY/2), 2) <= radius*radius&&pow(diff(x, LX/2), 2) + pow(diff(y, LY/2), 2) >= ((45)*(45)))
      nematicOrder = 1;
  }
  else if(init_config=="square")
  {
    if (diff(LY/2, y) < level/2 && diff(LX/2, x) < level/2)
      nematicOrder = 1;
  }
  else if(init_config=="stripe")
  {
    if(diff(LY/2, y) < level/2) nematicOrder = 1;
  }
  else
    throw error_msg("error: initial configuration '", init_config, "' unknown.");

  theta   = angle + noise*M_PI*(random_real() - .5);
  ftheta   = angle + noise*M_PI*(random_real() - .5);
	
  QQxx[k] = nematicOrder*(cos(2*theta));
  QQyx[k] = nematicOrder*(sin(2*theta));

  fQQxx[k] = finitorder*(cos(2*ftheta));
  fQQyx[k] = finitorder*(sin(2*ftheta));

  phi[k]  = 1.001;// nematicOrder*conc + noise*(random_real() - .5);
  totalphi += phi[k];
  fusion[k] = 0;
  // equilibrium dist
  ux[k] = uy[k] = ux_phi[k] = uy_phi[k] = 0;
  // compute totals for later checks
  ptot += phi[k];
  outS[k] = 0;
  gradS[k] = 0;
  maxstress[k] = 0;
  L1v[k] = K1;
  L2v[k] = K2;
  L3v[k] = K3;
  fusionSum[k] = 0;
}

void DryLyotropicF::ConfigureAtNode1(unsigned k,vector<float> Qxxinits, vector<float> Qyxinits)
{
  const unsigned x = GetXPosition(k);
  const unsigned y = GetYPosition(k);

  QQxx[k] = Qxxinits[k];
  QQyx[k] = Qyxinits[k];
  phi[k]  = 1.01;
  totalphi += phi[k];

  ux[k] = uy[k] = ux_phi[k] = uy_phi[k] = 0;
  // compute totals for later checks
  ptot += phi[k];
  fusion[k] = 0;
  outS[k] = 0;
  L1v[k] = K1;
  L2v[k] = K2;
  L3v[k] = K3;
  fusionSum[k] = 0;
}

void DryLyotropicF::Configure()
{
  string fileQxx = "/groups/astro/ardash/Documents/newmass/initconds/Qxx.txt";
  string fileQyx = "/groups/astro/ardash/Documents/newmass/initconds/Qxy.txt";

  string line;
  vector<float> Qxxinit = ParseInitConds(fileQxx);
  vector<float> Qyxinit = ParseInitConds(fileQyx);

  if (inits == 1){
    cout << "\n Using file to set initial conditions\n";
    for(unsigned k=0; k<DomainSize; ++k){
      ConfigureAtNode1(k, Qxxinit, Qyxinit);
    }
  }
  else{
    cout << "\n Using default initial conditions\n";
    for(unsigned k=0; k<DomainSize; ++k){
      ConfigureAtNode(k);
    }
  }

  //and do preinitialization
  cout << "Preinitialization started. ... ";
  preinit_flag = true;
  for (int i = 0; i< n_preinit; i++){
    Step();
  }
  preinit_flag = false;
  cout << "Preinitialization done. ... ";

}

vector<float>DryLyotropicF::ParseInitConds(string filename)
{
   ifstream infile;
   infile.open(filename);
   string line;
   vector<float> parsedCsv;
   while(std::getline(infile,line))
    {
      stringstream lineStream(line);
      string cell;
      size_t offset = 0;
      while(std::getline(lineStream,cell,','))
      {
        float celldbl = stof(cell,0);
        parsedCsv.push_back(celldbl);
      }
    }
   return parsedCsv; 
}

void DryLyotropicF::CalcMaxStress()
{
  std::vector<double> stressvals;
    for(unsigned kk=0; kk<DomainSize; ++kk){
      stressvals.push_back(outS[kk]);
    }
  auto min_value = *std::min_element(stressvals.begin(), stressvals.end());


  maxstress[1] = double(min_value);

}

void DryLyotropicF::CalcTotalFusion()
{
  double totalfusion = 0;
  for(unsigned kk=0; kk<DomainSize; ++kk){
    totalfusion = totalfusion + fusion[kk];
  }
  
  //auto sum_val = std::accumulate(totalfusion.begin(), totalfusion.end(),0);

  fusionSum[1] = fusionSum[1] +  totalfusion;
}

void DryLyotropicF::UpdateNematicQuantitiesAtNode(unsigned k)
{
  const auto& d = get_neighbours(k);  

  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k]; 
   
  const double del2Qxx  = laplacian(QQxx,  d, sD);
  const double dxQxx    = derivX   (QQxx,  d, sB);
  const double dyQxx    = derivY   (QQxx,  d, sB);
  const double del2Qyx  = laplacian(QQyx,  d, sD);
  const double dxQyx    = derivX   (QQyx,  d, sB);
  const double dyQyx    = derivY   (QQyx,  d, sB);

  const double p        = phi[k];
  const double dxPhi    = derivX   (phi, d, sB);
  const double dyPhi    = derivY   (phi, d, sB);
  const double del2p    = laplacian(phi, d, sD);

  const double dxQxx2 = dxQxx*dxQxx;
  const double dyQxx2 = dyQxx*dyQxx;
  const double dxQyx2 = dxQyx*dxQyx;
  const double dyQyx2 = dyQyx*dyQyx;

  const double dxdxQxx = derivXX (QQxx, d, sD);
  const double dydxQxx = derivXY (QQxx, d, sD);
  const double dxdxQyx = derivXX (QQyx, d, sD);
  const double dydxQyx = derivXY (QQyx, d, sD);
  const double dydyQxx = derivYY (QQxx, d, sD);
  const double dydyQyx = derivYY (QQyx, d, sD);
  const double dxdyQxx = derivXY (QQxx, d, sD);
  const double dxdyQyx = derivXY (QQyx, d, sD);

  //const double dxIsoS    = derivX(outS,  d, sB);
  //const double dyIsoS    = derivY(outS,  d, sB);
  //const double gradSmag = sqrt(dxIsoS*dxIsoS + dyIsoS*dyIsoS);

  //add fusion-dependent elasticity

  double mag = 5;//7.5;

  double K1now = L1v[k];
  double K2now = L2v[k];
  double K3now = L3v[k];

  //step-wise elasticity increase
  //double K1nownew = K1now;
  //double K2nownew = K2now;
  //double K3nownew = K3now;
 
  //if (time_step >= ton){
  //  K1nownew = K1*mag;
  //  K2nownew = K2*mag;
  //  K3nownew = K3*mag;
  //}

  //K1now = K1nownew;
  //K2now = K2nownew;
  //K3now = K3nownew;



  //fusion-driven, spatially-dependent elasticity increase
  //if (time_step >= ton){
  //  K1now= K1now + epsilon*fusion[k];
  //  K2now = K2now + epsilon*fusion[k];
  //  K3now= K3now + epsilon*fusion[k];

  //  if (K1now >= mag*K1){
  //    K1now = mag*K1;
  //    K2now = mag*K2;
  //    K3now = mag*K3;
  //  }
  //}


  //non-spatial elasticity increase, NO FUSION
  if (time_step > ton){
    K1now= K1now + epsilon*K1now;
    K2now = K2now + epsilon*K2now;
    K3now= K3now + epsilon*K3now;

    //K1now = K1 + epsilon*(time_step-ton);
    //K2now = K2 + epsilon*(time_step-ton);
    //K3now = K3 + epsilon*(time_step-ton);
  }

  if (K1now > mag*K1){
    K1now = mag*K1;
    K2now = mag*K2;
    K3now = mag*K3;
  }

  // calculate L elastic constants (S is q_eq)
  double const S = 1;
  double L1 = 0.25*(2*K2now + K3now-K1now)/(S*S);
  double L2 = (K1now-K2now)/(S*S);
  double L3 = 0.5*(K3now-K1now)/(S*S*S);


  

  L1v[k] = K1now;
  L2v[k] = K2now;
  L3v[k] = K3now;


  const double term = 1.0 - Qxx*Qxx - Qyx*Qyx;
  const double mixed_squared = dxQxx2 + dyQxx2 + dxQyx2 + dyQyx2;
  const double pm_mixed_squared = dxQxx2 - dyQxx2 + dxQyx2 - dyQyx2;
  const double free_energy = (L1+0.5*L2)*mixed_squared + L2*(dxQxx*dyQyx-dyQxx*dxQyx)
               + L3*(Qxx*pm_mixed_squared + Qyx*(dxQxx*dyQxx+dxQyx*dyQyx));

  const double Hxx = CC*term*Qxx + (L1+0.5*L2)*del2Qxx
           + L3*(Qxx*(dxdxQxx-dydyQxx) + 2*Qyx*dxdyQxx + dxQxx*dyQyx + dyQxx*dxQyx
              + 0.5*(dxQxx2-dyQxx2-dxQyx2+dyQyx2));
  const double Hyx = CC*term*Qyx + (L1+0.5*L2)*del2Qyx
           + L3*(Qxx*(dxdxQyx-dydyQyx) + 2*Qyx*dxdyQyx
              +dxQxx*(-dyQxx+dxQyx) + dyQyx*(-dyQxx+dxQyx));



  const double mu = 0*(surftension_on? Aphi*p*(p-1.)*(p+1.) + CC*term - Kphi*del2p : Aphi*(p-1.) + 0.0*CC*term - 0*Kphi*del2p);// 

  const double sigmaB = (free_energy +       //(!no_backflow)*
                    - 2*L1*(mixed_squared)
                    + L2*(-mixed_squared-2*dxQxx*dyQyx+2*dyQxx*dxQyx)
                    + L3*(2*Qxx*-pm_mixed_squared-4*Qyx*dxQxx*dyQxx-4*Qyx*dxQyx*dyQyx))
                    + .5*CC*term*term;
                    //+ 0*.5*Aphi*(p-1.)*(p-1.) - 0*mu*(p-1) + .5*CC*term*term;

  const double sigmaF = 1*(2*xi*( (Qxx*Qxx-1.)*Hxx + Qxx*Qyx*Hyx )
            +(L1+0.5*L2)*(-pm_mixed_squared) - Qxx*L3*mixed_squared)
            - zeta*Qxx;
            //preinit_flag? 0 :  (p > 1.0 ? 1 : 0)*zeta*Qxx ) + 0*.5*Kphi*(dyPhi*dyPhi-dxPhi*dxPhi);

  const double sigmaS = 1*(2*xi*(Qyx*Qxx*Hxx + (Qyx*Qyx-1)*Hyx)
            +(-(2*L1+L2)*(dxQxx*dyQxx+dxQyx*dyQyx) - Qyx*L3*mixed_squared))
            - zeta*Qyx;
            //- (preinit_flag? 0 :  (p > 1.0 ? 1 : 0)*zeta*Qyx ) - 0*Kphi*dxPhi*dyPhi;

  const double sigmaA = 1*(2*(Qxx*Hyx - Qyx*Hxx)
            +L3*(Qyx*-pm_mixed_squared + 2*Qxx*(dxQxx*dyQxx+dxQyx*dyQyx)));              

  // transfer to arrays  
  HHxx[k]    =  Hxx;
  HHyx[k]    =  Hyx;
  dxQQxx[k]  =  dxQxx;
  dxQQyx[k]  =  dxQyx;
  dyQQxx[k]  =  dyQxx;
  dyQQyx[k]  =  dyQyx;
  del2QQxx[k] = del2Qxx;
  del2QQyx[k] = del2Qyx;
  sigmaXX[k] =  sigmaF + sigmaB;
  sigmaYY[k] = -sigmaF + sigmaB;
  sigmaXY[k] =  sigmaS + sigmaA;
  sigmaYX[k] =  sigmaS - sigmaA;
  MU[k]      =  mu;
  outS[k]    = 0.5*(sigmaXX[k] + sigmaYY[k]) - ds;

  gradS[k] = sigmaXX[k] - sigmaYY[k];

}
void DryLyotropicF::UpdateNematicQuantities()
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateNematicQuantitiesAtNode(k);
}

void DryLyotropicF::UpdateFluidQuantitiesAtNode(unsigned k)
{
  // array placeholders for current node
  const auto& d = get_neighbours(k);

  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k];
  const double fQxx = fQQxx[k];
  const double fQyx = fQQyx[k]; 

  const double dxSxx = derivX(sigmaXX, d, sB);
  const double dySxy = derivY(sigmaXY, d, sB);
  const double dxSyx = derivX(sigmaYX, d, sB);
  const double dySyy = derivY(sigmaYY, d, sB);
  const double Fx = dxSxx + dySxy;
  const double Fy = dxSyx + dySyy;

  double fS = sqrt(fQyx*fQyx + fQxx*fQxx);
  
  FFx[k] = Fx;
  FFy[k] = Fy;

  const double p   = phi[k];
  // compute velocities
  const double vx = (FFx[k])/(friction + 1*friction*fS);
  const double vy = (FFy[k])/(friction + 1*friction*fS);

   //polarity: polarity picks the head or tail direciton of the nematic director based on ailgnment to local velocity
  double S = sqrt(Qyx*Qyx + Qxx*Qxx);
  double nx = sqrt((1+Qxx/S)/2);
  double sgn = 1;
  if (Qyx<0){
    sgn = -1;
  }
  double ny = sgn*sqrt((1-Qxx/S)/2);

  double tp = vx*nx + vy*ny;
  double btm = sqrt(vx*vx+vy*vy)*sqrt(nx*nx+ny*ny);
  double th = acos(tp/btm);

  double dirchange = 0;
  if (th < M_PI/2 || th >= 3*M_PI/2){
    dirchange = 0; //Q = (qx,qy) is already pointing in the direction of U
  }
  else{
    dirchange = M_PI; //Q = (qx,qy) is pointing in the opposite direction of U, need to change
  }

  double Qangle = atan2(ny,nx) + dirchange;

  Px[k] = cos(Qangle);
  Py[k] = sin(Qangle);

  const double px = Px[k];
  const double py = Py[k];

  // transfer to arrays
  ux[k]      =  vx;// + (p > 1.0 ? 1 : 0.5)*alpha*px;
  uy[k]      =  vy;// + (p > 1.0 ? 1 : 0.5)*alpha*py;
  ux_phi[k]  =  0;//(vx + (p > 1.0 ? 1 : 0.5)*alpha*px)*p + (p > 1.0 ? 1 : 0.5)*Vphi*px*p;
  uy_phi[k]  =  0;//(vy + (p > 1.0 ? 1 : 0.5)*alpha*px)*p + (p > 1.0 ? 1 : 0.5)*Vphi*py*p;
}
void DryLyotropicF::UpdateFluidQuantities()
{
  double sum = 0;

  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k){
    sum = sum + phi[k];
    UpdateFluidQuantitiesAtNode(k);
  }
  countphi = sum;
}

void DryLyotropicF::UpdateNematicFieldsAtNode(unsigned k, bool first)
{  
  const auto& d = get_neighbours(k);
  
  const double vx = ux[k];
  const double vy = uy[k];

  const double Qxx = QQxx[k];
  const double Qyx = QQyx[k];

  const double fQxx = fQQxx[k];
  const double fQyx = fQQyx[k];

  const double Hxx = HHxx[k];
  const double Hyx = HHyx[k];

  const double dxQxx = dxQQxx[k];
  const double dyQxx = dyQQxx[k];
  const double dxQyx = dxQQyx[k];
  const double dyQyx = dyQQyx[k];

  const double dxux = derivX(ux, d, sB);
  const double dyux = derivY(ux, d, sB);
  const double dxuy = derivX(uy, d, sB);
  const double dyuy = derivY(uy, d, sB);

  const double expansion = dxux + dyuy;
  const double shear     = .5*(dxuy + dyux);
  const double vorticity = .5*(dxuy - dyux);
  const double traceQL   = Qxx*(dxux - dyuy) + 2*Qyx*shear;

  //align to fibernectin?
  double Falign=0;
  if (time_step>tf){     //yes
    Falign=1;
  }

  //components of the Beris-Edwards equation
  const double Dxx = 1*(GammaQ*(Hxx - gam*laplacian(del2QQxx, d, sD)) - vx*dxQxx - vy*dyQxx - 2*vorticity*Qyx
    + xi*((Qxx+1)*(2*dxux-traceQL) +2*Qyx*shear -expansion)) + Falign*Kf*fQxx;

  const double Dyx = 1*(GammaQ*(Hyx - gam*laplacian(del2QQyx, d, sD)) - vx*dxQyx - vy*dyQyx + 2*vorticity*Qxx
    + xi*( Qyx*(expansion-traceQL) + 2*shear)) + Falign*Kf*fQyx;
  
  if(first)
  {
    QNxx[k] = QQxx[k] + .5*Dxx;
    QNyx[k] = QQyx[k] + .5*Dyx;
    QQxx[k] = QNxx[k] + .5*Dxx;
    QQyx[k] = QNyx[k] + .5*Dyx;
  }
  else
  {
    QQxx[k] = QNxx[k] + .5*Dxx;
    QQyx[k] = QNyx[k] + .5*Dyx;
  }


 //updating fibronectin Q
  //if (time_step< tf){
  //  fQQxx[k] = 0;
  //  fQQyx[k] = 0;
  //}
 // else {
 //  fQQxx[k] = QQxx[k-ftau];
 //  fQQyx[k] = QQyx[k-ftau];
 // }


}


void DryLyotropicF::ReadJson(){
  double curr_M = Mval[1];
  for(unsigned kk=0; kk<DomainSize; ++kk){
    const double fQxx = fQQxx[kk];
    const double fQyx = fQQyx[kk]; 

    fQQxx[kk] = fQxx;
    fQQyx[kk] = fQyx;
  }



    if (time_step > tf){

      if ((int)time_step % fntau== 0){
        cout << (int)time_step << "," << time_step-ftau << "\n";
        std::string frameno = std::to_string((int)time_step-(int)ftau);
        cout << "/groups/astro/ardash/Documents/newmass/out/frame"+frameno+".json" <<"\n";
        std::ifstream file("/groups/astro/ardash/Documents/newmass/out/frame"+frameno+".json");
        if (!file) {
          std::cerr << "Could not open the file!" << std::endl;
        }

        json j;
        file >> j;

        std::map<std::string, std::vector<double>> data_arrays;

        data_arrays["QQxx"] = j["data"]["QQxx"]["value"].get<std::vector<double>>();
        data_arrays["QQyx"] = j["data"]["QQyx"]["value"].get<std::vector<double>>();

        for(unsigned k=0; k<DomainSize; ++k){
          const double Qxx = data_arrays["QQxx"][k];
          const double Qyx = data_arrays["QQyx"][k];

          double S = sqrt(Qyx*Qyx + Qxx*Qxx);
          double nx = sqrt((1+Qxx/S)/2);
          double sgn = 1;
          if (Qyx<0){
            sgn = -1;
          }
          double ny = sgn*sqrt((1-Qxx/S)/2);

          double fS = curr_M*S;
          if (fS >= 1.0){
            fS = 1.0;
            Mval[1] = 1.0;
          }
          else{
            Mval[1] = curr_M + M;
          }

          fQQxx[k] = fS*(2.0*nx*nx-1.0);
          fQQyx[k] = 2.0*fS*nx*ny;
        }
      }
  }

}

void DryLyotropicF::UpdateNematicFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateNematicFieldsAtNode(k, first);
}

void DryLyotropicF::UpdateFluidFieldsAtNode(unsigned k, bool first)
{  
  const auto& d = get_neighbours(k);
  //const double sXX = sigmaXX(k);
  //const double sYY = sigmaYY(k);
  
  const double del2mu = laplacian(MU, d, sD);
  const double pFlux = derivX(ux_phi, d, sB) + derivY(uy_phi, d, sB);
  const double Dp = 0*(GammaPhi*del2mu - pFlux - ( conserve_phi ? (countphi-totalphi)/DomainSize : 0 ));

  //const double isostress = 0.5*(sXX+sYY);

  //get average stress value


  if (time_step != tmix){
    if (first){
      phn[k]     = phi[k]  + .5*Dp;
      phi_tmp[k] = phi[k]  +    Dp;
    }
    else {
      phi_tmp[k] = phn[k]  + .5*Dp;
    }
  }
  else{
      double valphi = conc + noise*(random_real() - .5);
      phi_tmp[k] = valphi;
      phi[k] = valphi;
      phn[k] = valphi;
      //cout << "here" << "\n";
  }
}


void DryLyotropicF::UpdateFusionFieldsAtNode(unsigned k, bool first)
{  
  const auto& d = get_neighbours(k);
  const double fus = fusion[k];
  const double isos = outS[k];

  //stress-dependent heaviside
  //const double Smax = maxstress[1];
  double Smax = smaxval;
  double hvs=0;

  if (isos < 0){
    if (isos/Smax >= thrStress){
      hvs = omega;
    } 
  }
  
  const double Df = 0*(aF*(1-fus/fC) + hvs - delta*fus); //to add Heaviside term after
  //cout << Df << '\n';

  double fusionon=0;
  if (time_step > ton ){
    fusionon = 1;
  }
  if (first){
      fusn[k]     = fusion[k];// + .5*Df*fusionon;
      fusion_tmp[k] = fusion[k];//  +    Df*fusionon;
    }
    else {
      fusion_tmp[k] = fusn[k];//  + .5*Df*fusionon;
    }

}




void DryLyotropicF::UpdateFluidFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateFluidFieldsAtNode(k, first);

  swap(phi.get_data(), phi_tmp.get_data());
}

void DryLyotropicF::UpdateFusionFields(bool first)
{
  #pragma omp parallel for num_threads(nthreads) if(nthreads)
  for(unsigned k=0; k<DomainSize; ++k)
    UpdateFusionFieldsAtNode(k, first);

  swap(fusion.get_data(), fusion_tmp.get_data());
}



void DryLyotropicF::BoundaryConditionsFields()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // channel
    case 1:
    case 2:
      QQxx.ApplyNeumannChannel();
      QQyx.ApplyNeumannChannel();
      phi .ApplyNeumannChannel();
      fusion.ApplyNeumannChannel();
      break;
    // box
    case 3:
    case 4:
      QQxx.ApplyNeumann();
      QQyx.ApplyNeumann();
      phi .ApplyNeumann();
      fusion.ApplyNeumann();
      break;
    // pbc with bdry layer
    default:
      QQxx.ApplyPBC();
      QQyx.ApplyPBC();
      phi.ApplyPBC();
      fusion.ApplyPBC();
  }
}
void DryLyotropicF::BoundaryConditionsFields2()
{
  switch(BC)
  {
    // pbc without bdry layer (nothing to do)
    case 0:
      break;
    // channel
    case 1:
    case 2:
      uy     .CopyDerivativeChannel();
      ux     .CopyDerivativeChannel();
      uy_phi .ApplyDirichletChannel(0);
      ux_phi .ApplyDirichletChannel(0);

      MU     .ApplyNeumannChannel();
      sigmaXX.ApplyNeumannChannel();
      sigmaYY.ApplyNeumannChannel();
      sigmaYX.ApplyNeumannChannel();
      sigmaXY.ApplyNeumannChannel();
      break;
    // box
    case 3:
    case 4:
      uy     .CopyDerivative();
      ux     .CopyDerivative();
      uy_phi .ApplyDirichlet(0);
      ux_phi .ApplyDirichlet(0);

      MU     .ApplyNeumann();
      sigmaXX.ApplyNeumann();
      sigmaYY.ApplyNeumann();
      sigmaYX.ApplyNeumann();
      sigmaXY.ApplyNeumann();
      break;
    // pbc with bdry layer
    default:
      ux     .ApplyPBC();
      uy     .ApplyPBC();
      ux_phi .ApplyPBC();
      uy_phi .ApplyPBC();
      MU     .ApplyPBC();
      sigmaXX.ApplyPBC();
      sigmaYY.ApplyPBC();
      sigmaYX.ApplyPBC();
      sigmaXY.ApplyPBC();
  }
}

void DryLyotropicF::Step()
{
  // boundary conditions for primary fields
  BoundaryConditionsFields();
  // predictor step
  UpdateNematicQuantities();
  UpdateFluidQuantities();
  
  BoundaryConditionsFields2();
  CalcMaxStress();
  CalcTotalFusion();
  ReadJson();

  // LB Step
  this->UpdateNematicFields(true); 
  this->UpdateFluidFields(true);  
  this->UpdateFusionFields(true); 

  // corrector steps
  for(unsigned n=1; n<=npc; ++n)
  {    
    BoundaryConditionsFields();
    UpdateNematicQuantities();
    UpdateFluidQuantities();
    BoundaryConditionsFields2();
    this->UpdateNematicFields(true); 
    this->UpdateFluidFields(true);
    this->UpdateFusionFields(true); 
  }
}

void DryLyotropicF::RuntimeChecks()
{
  // check that phi is conserved
  {
    double pcheck = 0;
    for(unsigned k=0; k<DomainSize; ++k)
        pcheck += phi[k];
    cout << "pcheck: " << pcheck << "/" << ptot << '\n';
    if(abs(ptot-pcheck)>1)
      throw error_msg("phi is not conserved (", ptot, "/", pcheck, ")");
  }
}

option_list DryLyotropicF::GetOptions()
{
  // model specific options
  opt::options_description model_options("Model options");
  model_options.add_options()
    ("GammaQ", opt::value<double>(&GammaQ),
     "Q-tensor mobility")
    ("GammaPhi", opt::value<double>(&GammaPhi),
     "binary mobility")
    ("xi", opt::value<double>(&xi),
     "tumbling/aligning parameter")
    ("friction", opt::value<double>(&friction),
     "friction from confinement")
    ("Aphi", opt::value<double>(&Aphi),
     "binary fluid bulk constant")
    ("Vphi", opt::value<double>(&Vphi),
     "Self-advection of the concentration")
    ("alpha", opt::value<double>(&alpha),
     "Selfpropulsion force")
    ("CC", opt::value<double>(&CC),
     "coupling constant")
    ("K1", opt::value<double>(&K1),
     "Frank elastic constant 1")
    ("K2", opt::value<double>(&K2),
     "Frank elastic constant 2")
    ("K3", opt::value<double>(&K3),
     "Frank elastic constant 2")
    ("Kphi", opt::value<double>(&Kphi),
     "binary gradient constant")
    ("zeta", opt::value<double>(&zeta),
     "activity parameter")
    ("npc", opt::value<unsigned>(&npc),
     "number of correction steps for the predictor-corrector method")
    ("backflow_on", opt::value<bool>(&backflow_on),
     "Backflow flag")
     ("surftension_on", opt::value<bool>(&surftension_on),
     "surface tension flag")
     ("tmix", opt::value<double>(&tmix),
     "time at which phi is mixed")
    ("n_preinit", opt::value<int>(&n_preinit),
     "number of preinitialization steps")
    ("inits", opt::value<int>(&inits),
     "use file as initial conds (1=yes, 0=no)")
    ("aF", opt::value<double>(&aF),
     "fusion growth rate")
    ("fC", opt::value<double>(&fC),
     "fusion carrying capacity")
    ("gam", opt::value<double>(&gam),
     "gam parameter for 4th deriv in Q (for numerics)")
    ("thrStress", opt::value<double>(&thrStress),
     "compressive stress threshold")
    ("omega", opt::value<double>(&omega),
     "fusion rate (Heaviside)")
    ("ton", opt::value<double>(&ton),
     "time at which fusion is on")
    ("epsilon", opt::value<double>(&epsilon),
     "fusion-elasticity coef")
    ("ds", opt::value<double>(&ds),
     "stress-ofset")
    ("delta", opt::value<double>(&delta),
     "sink term for fusion")
    ("smaxval", opt::value<double>(&smaxval),
     "max stress value (cheating)")
    ("finitorder", opt::value<double>(&finitorder),
     "initial order of fibronectin")
    ("ftau", opt::value<unsigned>(&ftau),
     "time delay in fibronectin dynamics")
    ("fM", opt::value<double>(&fM),
     "rate of increase of the fibronectin order parameter (Sf)")
    ("tf", opt::value<unsigned>(&tf),
     "time when fibronectin is turned on (tf > ftau!!!)")
    ("fntau", opt::value<unsigned>(&tf),
     "number of steps between fibronecting director update")
    ("M", opt::value<double>(&M),
     "increase in the fibronectin order")
    ("Kf", opt::value<double>(&Kf),
     "fibronectin alignment parameter");;

  // init config options
  opt::options_description config_options("Initial configuration options");
  config_options.add_options()
    ("config", opt::value<string>(&init_config),
     "initial configuration")
    ("level", opt::value<double>(&level),
     "starting thickness of the nematic region")
    ("conc", opt::value<double>(&conc),
     "starting phi concentration of nematic region")
    ("radius", opt::value<double>(&radius),
     "radius of the initial circle")
    ("angle", opt::value<double>(&angle_deg),
     "initial angle to x direction (in degrees)")
    ("noise", opt::value<double>(&noise),
     "size of initial variations");

  return { model_options, config_options };
}
