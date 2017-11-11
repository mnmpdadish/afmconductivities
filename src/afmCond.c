//
//  afmCond.c
//

#include "utilities.h"

int main(int argc, const char * argv[]) {
    printf("afmCond starting\n\n");
    
    ///////////////////////////// hamiltonian parameters (read from file) ////////////////////////////
    //hamiltonian parameters:

    double ETA=0.1;
    double t=1.0;
    double tp=-0.3;
    double tpp=0.2;
    double M=1.5;
    
    double beta=100.0;
    int nOmega = 11;
    int nK = 401;
    double amplitudeCutoff = 0.005;

    int nMu = 20;
    double muMin=-4.0; double muMax=4.0;
    
    
    // read parameters form file: /////////////////////////////
    FILE * file = fopen("model.dat", "rt");
    if(file == NULL) {printf("file %s not found", "model.dat"); exit(1);}
    printf("reading parameters from model.dat\n\n") ;
            
    readDouble(file, "tpp",  &tpp);
    readDouble(file, "tp",   &tp);
    readDouble(file, "t",    &t);
    readDouble(file, "ETA",  &ETA);
    readDouble(file, "M",    &M);
    readDouble(file, "beta", &beta);
    readDouble(file, "muMin",&muMin);
    readDouble(file, "muMax",&muMax);
    readDouble(file, "amplitudeCutoff",  &amplitudeCutoff);
    
    readInt(file, "nOmega", &nOmega);
    readInt(file, "nK",     &nK);
    readInt(file, "nMu",    &nMu);
    
    fclose(file);
    ///////////////////////////////////////////////////////////
    
    
    
    // precalculate the omega vector and the derivative of the Fermi Dirac vector:
    double energyCutoff = 2.*2.*acosh(0.25*sqrt(beta/amplitudeCutoff)) /beta; // we put an additionnal factor of 2. for the derivative cutoff
    
    double omega[nOmega];
    //double fermiDirac_dw[nOmega];
    double dfermiDirac_dw[nOmega];
    
    int n=0; for(n=0; n<nOmega; n++)
    {
        omega[n]= -energyCutoff + 2.*energyCutoff*n/(nOmega-1);
        double expBw = exp(beta*omega[n]);
        dfermiDirac_dw[n]= -beta*expBw/((expBw+1.)*(expBw+1.));
    }
    
    
    FILE *fileOut = fopen("conductivities.dat","w");
    
    printf("#mu          density      sigma1_xx       sigma2_xx       sigma1_xy       sigma2_xy");

    fprintf(fileOut, "#mu          density      ");
    fprintf(fileOut, "sigma1_xx       sigma2_xx       sigma1_xy       sigma2_xy       ");
    fprintf(fileOut, "alpha1_xx       alpha2_xx       alpha1_xy       alpha2_xy       ");
    fprintf(fileOut, "beta1_xx        beta2_xx        beta1_xy        beta2_xy\n");
    
    
    double sink[nK]; double sin2k[nK]; 
    double cosk[nK]; double cos2k[nK]; 
    
    int i=0; for(i=0; i<nK; i++)
    {
       double k = M_PI*(-1.0 + i*2.0/nK);

       //precalculate some trigo stuff:
       sink[i] = sin(k); sin2k[i] = sin(2.*k); 
       cosk[i] = cos(k); cos2k[i] = cos(2.*k);
    }
    
    
    
    int m=0; for(m=0; m<nMu; m++)
    {
       double mu = muMin + m*(muMax-muMin)/(nMu-1);
       // initialise the sums:
       double sigma1_xx = 0., sigma2_xx = 0., sigma1_xy = 0., sigma2_xy = 0.;
       double alpha1_xx = 0., alpha2_xx = 0., alpha1_xy = 0., alpha2_xy = 0.;
       double beta1_xx = 0.,  beta2_xx = 0.,  beta1_xy = 0.,  beta2_xy = 0.;
       double density = 0.;
       
       for(i=0; i<nK; i++)
       {
          int j=0; for(j=0; j<nK; j++)
          {
              ///////////////////////////// calculate eigenenergies and its derivatives ////////////////////////////
              //double kx = M_PI*(-1.0 + i*2.0/nK);
              //double ky = M_PI*(-1.0 + j*2.0/nK);

              //precalculate some trigo stuff:
              //double sinkx = sin(kx), sinky = sin(ky), sin2ky = sin(2.*ky), sin2kx = sin(2.*kx); 
              //double coskx = cos(kx), cosky = cos(ky), cos2ky = cos(2.*ky), cos2kx = cos(2.*kx);
              
              // dispersion relation (and its derivatives):
              // note: kx = k[i] and ky = k[j]
              double epsilon_k           = -2.*t*(cosk[i] + cosk[j])- 4.*tp*cosk[i]*cosk[j] - 2.*tpp*(cos2k[i] + cos2k[j]) - mu;
              double depsilon_k_dkx      =  2.*t*sink[i]             + 4.*tp*sink[i]*cosk[j] + 4.*tpp*sin2k[i];
              double depsilon_k_dky      =  2.*t*sink[j]             + 4.*tp*cosk[i]*sink[j] + 4.*tpp*sin2k[j];
              //double ddepsilon_k_dkx_dkx =  2.*t*cosk[i]           + 4.*tp*cosk[i]*cosk[j] + 8.*tpp*cos2k[i];
              double ddepsilon_k_dky_dky =  2.*t*cosk[j]             + 4.*tp*cosk[i]*cosk[j] + 8.*tpp*cos2k[j];
              double ddepsilon_k_dkx_dky =                            - 4.*tp*sink[i]*sink[j];

              double epsilon_kQ           =  2.*t*(cosk[i] + cosk[j])- 4.*tp*cosk[i]*cosk[j] - 2.*tpp*(cos2k[i] + cos2k[j]) - mu;
              double depsilon_kQ_dkx      = -2.*t*sink[i]             + 4.*tp*sink[i]*cosk[j] + 4.*tpp*sin2k[i];
              double depsilon_kQ_dky      = -2.*t*sink[j]             + 4.*tp*cosk[i]*sink[j] + 4.*tpp*sin2k[j];
              //double ddepsilon_k_dkx_dkx = -2.*t*cosk[i]           + 4.*tp*cosk[i]*cosk[j] + 8.*tpp*cos2k[i];
              double ddepsilon_kQ_dky_dky = -2.*t*cosk[j]             + 4.*tp*cosk[i]*cosk[j] + 8.*tpp*cos2k[j];
              double ddepsilon_kQ_dkx_dky =                            - 4.*tp*sink[i]*sink[j];

              //precalculate sum, diff and radical:
              double Sk          = 0.5*(  epsilon_k         +   epsilon_kQ);
              double dSk_dkx     = 0.5*( depsilon_k_dkx     +  depsilon_kQ_dkx);
              double dSk_dky     = 0.5*( depsilon_k_dky     +  depsilon_kQ_dky);
              //double ddSk_dkx_dkx= 0.5*(ddepsilon_k_dkx_dkx + ddepsilon_kQ_dkx_dkx);
              double ddSk_dky_dky= 0.5*(ddepsilon_k_dky_dky + ddepsilon_kQ_dky_dky);
              double ddSk_dkx_dky= 0.5*(ddepsilon_k_dkx_dky + ddepsilon_kQ_dkx_dky);

              double Dk          = 0.5*(  epsilon_k         -   epsilon_kQ);
              double dDk_dkx     = 0.5*( depsilon_k_dkx     -  depsilon_kQ_dkx);
              double dDk_dky     = 0.5*( depsilon_k_dky     -  depsilon_kQ_dky);
              //double ddDk_dkx_dkx= 0.5*(ddepsilon_k_dkx_dkx - ddepsilon_kQ_dkx_dkx);
              double ddDk_dky_dky= 0.5*(ddepsilon_k_dky_dky - ddepsilon_kQ_dky_dky);
              double ddDk_dkx_dky= 0.5*(ddepsilon_k_dkx_dky - ddepsilon_kQ_dkx_dky);

              double Rk, dRk_dkx, dRk_dky, ddRk_dky_dky, ddRk_dkx_dky;
              //double ddRk_dkx_dkx;
              
              if(M<=0.00001){
                  Rk = Dk; 
                  dRk_dkx = dDk_dkx;
                  dRk_dky = dDk_dky;
                  //ddRk_dkx_dkx = ddDk_dkx_dkx; 
                  ddRk_dky_dky = ddDk_dky_dky;
                  ddRk_dkx_dky = ddDk_dkx_dky;
              }
              else{
                  Rk          = sqrt( Dk*Dk + M*M );
                  dRk_dkx     = Dk*dDk_dkx/Rk;
                  dRk_dky     = Dk*dDk_dky/Rk;
                  //ddRk_dkx_dkx= ((dDk_dkx*dDk_dkx + Dk*ddDk_dkx_dkx) - (Dk*Dk*dDk_dkx*dDk_dkx) / (Rk*Rk))/Rk;
                  ddRk_dky_dky= ((dDk_dky*dDk_dky + Dk*ddDk_dky_dky) - (Dk*Dk*dDk_dky*dDk_dky) / (Rk*Rk))/Rk;
                  ddRk_dkx_dky= ((dDk_dkx*dDk_dky + Dk*ddDk_dkx_dky) - (Dk*Dk*dDk_dkx*dDk_dky) / (Rk*Rk))/Rk;
              }
              
              
              //finally calculate the eigen values and their derivatives (vertices):
              double E1_k          =   Sk         +   Rk ;
              double dE1_k_dkx     =  dSk_dkx     +  dRk_dkx;
              double dE1_k_dky     =  dSk_dky     +  dRk_dky;
              //double ddE1_k_dkx_dkx= ddSk_dkx_dkx + ddRk_dkx_dkx;
              double ddE1_k_dky_dky= ddSk_dky_dky + ddRk_dky_dky;
              double ddE1_k_dkx_dky= ddSk_dkx_dky + ddRk_dkx_dky;
              
              double E2_k          =   Sk         -   Rk ;
              double dE2_k_dkx     =  dSk_dkx     -  dRk_dkx;
              double dE2_k_dky     =  dSk_dky     -  dRk_dky;
              //double ddE2_k_dkx_dkx= ddSk_dkx_dkx - ddRk_dkx_dkx;
              double ddE2_k_dky_dky= ddSk_dky_dky - ddRk_dky_dky;
              double ddE2_k_dkx_dky= ddSk_dkx_dky - ddRk_dkx_dky;
              
              double kernel1_xx = dE1_k_dkx*dE1_k_dkx;
              double kernel2_xx = dE2_k_dkx*dE2_k_dkx;
              //double kernel1_xy = -(1./3.)*(dE1_k_dkx*dE1_k_dkx*ddE1_k_dky_dky + dE1_k_dky*dE1_k_dky*ddE1_k_dkx_dkx - 2.0*dE1_k_dkx*dE1_k_dky*ddE1_k_dkx_dky);
              //double kernel2_xy = -(1./3.)*(dE2_k_dkx*dE2_k_dkx*ddE2_k_dky_dky + dE2_k_dky*dE2_k_dky*ddE2_k_dkx_dkx - 2.0*dE2_k_dkx*dE2_k_dky*ddE2_k_dkx_dky);
              double kernel1_xy = -(2./3.)*(dE1_k_dkx*(dE1_k_dkx*ddE1_k_dky_dky - dE1_k_dky*ddE1_k_dkx_dky));
              double kernel2_xy = -(2./3.)*(dE2_k_dkx*(dE2_k_dkx*ddE2_k_dky_dky - dE2_k_dky*ddE2_k_dkx_dky));

              
              for(n=0; n<nOmega; n++)
              {
                  double complex z = omega[n] + ETA * I;
                  
                  double A1_k = -(1./M_PI)*cimag(1.0/ (z-E1_k) );
                  double A2_k = -(1./M_PI)*cimag(1.0/ (z-E2_k) );

                  double frequencyKernel1_xx = -dfermiDirac_dw[n]*kernel1_xx*A1_k*A1_k;
                  double frequencyKernel2_xx = -dfermiDirac_dw[n]*kernel2_xx*A2_k*A2_k;
                  double frequencyKernel1_xy = -dfermiDirac_dw[n]*kernel1_xy*A1_k*A1_k*A1_k;
                  double frequencyKernel2_xy = -dfermiDirac_dw[n]*kernel2_xy*A2_k*A2_k*A2_k;
                  
                  sigma1_xx += frequencyKernel1_xx;
                  sigma2_xx += frequencyKernel2_xx;
                  sigma1_xy += frequencyKernel1_xy;
                  sigma2_xy += frequencyKernel2_xy;
                              
                  alpha1_xx += beta*omega[n] * frequencyKernel1_xx;
                  alpha2_xx += beta*omega[n] * frequencyKernel2_xx;
                  alpha1_xy += beta*omega[n] * frequencyKernel1_xy;
                  alpha2_xy += beta*omega[n] * frequencyKernel2_xy;
                  
                  double omega2 = beta*beta* omega[n] * omega[n];
                  beta1_xx += omega2 * frequencyKernel1_xx;
                  beta2_xx += omega2 * frequencyKernel2_xx;
                  beta1_xy += omega2 * frequencyKernel1_xy;
                  beta2_xy += omega2 * frequencyKernel2_xy;
              }
              
              density += 1.0/(1.0+exp(beta*E1_k)) + 1.0/(1.0+exp(beta*E2_k));
          }
       }

       double f0 = 1.0/nK/nK;
       double f = (2.*energyCutoff) /(nOmega)*f0;

       printf("\n% 4.8f % 4.8f  % 4.8e % 4.8e % 4.8e % 4.8e", mu, 1.0-f0*density,  f*sigma1_xx, f*sigma2_xx, f*sigma1_xy, f*sigma2_xy);
       fprintf(fileOut,"% 4.8f % 4.8f  ", mu, 1.0-f0*density);
       fprintf(fileOut,"% 4.8e % 4.8e % 4.8e % 4.8e ", f*sigma1_xx, f*sigma2_xx, f*sigma1_xy, f*sigma2_xy);
       fprintf(fileOut,"% 4.8e % 4.8e % 4.8e % 4.8e ", f*alpha1_xx, f*alpha2_xx, f*alpha1_xy, f*alpha2_xy);
       fprintf(fileOut,"% 4.8e % 4.8e % 4.8e % 4.8e ", f*beta1_xx,  f*beta2_xx,  f*beta1_xy,  f*beta2_xy);
       fprintf(fileOut,"\n");
    }
    fclose(fileOut);
 
#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KBOLD "\033[1m"
   
    printf("\n\x1B[32m\033[1mafmCond over.\x1B[0m\n\n");
    return 0;
}
