//
//  main.cpp
//  oneBody
//

#define _USE_MATH_DEFINES
#include "utilities.h"

using namespace std;


int main(int argc, const char * argv[]) {
    printf("afmCond starting\n\n");
    
    
    
    ///////////////////////////// hamiltonian parameters (read from file) ////////////////////////////
    //hamiltonian parameters:
    double MU=0.0;
    double ETA=0.1;
    double t=1.0;
    double tp=0.0;
    double tpp=0.0;
    double M=0.0;
    
    double beta=100.0;
    int nOmega = 21;
    int nK = 401;
    double amplitudeCutoff = 0.005;
    
    //read file para.dat:
    if (not exists("model.dat")) {printf("ERROR: couldn't find file 'model.dat'\n\n"); exit(1);}
    printf("reading parameters from model.dat\n\n") ;
    ifstream file;
    file.open("model.dat");

    readNumber(file,"MU",MU); 
    readNumber(file,"ETA",ETA);
    readNumber(file,"t",t);
    readNumber(file,"tp",tp);
    readNumber(file,"tpp",tpp);
    readNumber(file,"M",M);
    readNumber(file,"beta",beta);
    
    readNumber(file,"nOmega",nOmega);
    readNumber(file,"nK",nK);
    readNumber(file,"amplitudeCutoff",amplitudeCutoff);
    
    file.close();
    
    
    
    // precalculate the omega vector and the derivative of the Fermi Dirac vector:
    double energyCutoff = 2.*2.*acosh(sqrt(beta/amplitudeCutoff)) /beta; // we put an additionnal factor of 2. for the derivative cutoff
    
    double omega[nOmega];
    //double fermiDirac_dw[nOmega];
    double dfermiDirac_dw[nOmega];
    for(int n=0; n<nOmega; n++)
    {
        omega[n]= -energyCutoff + 2.*energyCutoff/(nOmega-1);
        double expBw = exp(beta*omega[n]);
        dfermiDirac_dw[n]= -beta*expBw/((expBw+1.)*(expBw+1.));
    }
    
    // initialise the sums:
    double sigma1_xx = 0., sigma2_xx = 0., sigma1_xy = 0., sigma2_xy = 0.0;
    double alpha1_xx = 0., alpha2_xx = 0., alpha1_xy = 0., alpha2_xy = 0.0;
    double beta1_xx  = 0.,  beta2_xx = 0.,  beta1_xy = 0.,  beta2_xy = 0.0;
    double density = 0.;
        
    for(int i=0; i<nK; i++)
    for(int j=0; j<nK; j++)
    {
        ///////////////////////////// calculate eigenenergies and its derivatives ////////////////////////////
        double kx = M_PI*(-1.0 + i*1.0/(nK-1));
        double ky = M_PI*(-1.0 + j*1.0/(nK-1));

        //precalculate some trigo stuff:
        double sinkx = sin(kx), sinky = sin(ky), sin2ky = sin(2.*ky), sin2kx = sin(2.*kx); 
        double coskx = cos(kx), cosky = cos(ky), cos2ky = cos(2.*ky), cos2kx = cos(2.*kx);
        //double cos2kyQy = cos(2.*(ky+M_PI)), cos2kxQx = cos(2.*(kx+M_PI)); // no need to define since cos(a+2*pi) = cos(a)

        // dispersion relation (and its derivatives):
        double epsilon_k           = -2.*t*(coskx + cosky)- 4.*tp*coskx*cosky - 2.*tpp*(cos2kx + cos2ky) - MU;
        double depsilon_k_dkx      =  2.*t*sinkx          + 4.*tp*sinkx*cosky - 4.*tpp*sin2kx;
        double depsilon_k_dky      =  2.*t*sinky          + 4.*tp*coskx*sinky - 4.*tpp*sin2ky;
        //double ddepsilon_k_dkx_dkx =  2.*t*coskx          + 4.*tp*coskx*cosky - 8.*tpp*cos2kx;
        double ddepsilon_k_dky_dky =  2.*t*cosky          + 4.*tp*coskx*cosky - 8.*tpp*cos2kx;
        //double ddepsilon_k_dkx_dky =                      - 4.*tp*sinkx*sinky;

        // dispersion relation (and its derivatives):
        double epsilon_kQ           =  2.*t*(coskx + cosky)- 4.*tp*coskx*cosky - 2.*tpp*(cos2kx + cos2ky) - MU;
        double depsilon_kQ_dkx      = -2.*t*sinkx          + 4.*tp*sinkx*cosky - 4.*tpp*sin2kx;
        double depsilon_kQ_dky      = -2.*t*sinky          + 4.*tp*coskx*sinky - 4.*tpp*sin2ky;
        //double ddepsilon_kQ_dkx_dkx = -2.*t*coskx          + 4.*tp*coskx*cosky - 8.*tpp*cos2kx;
        double ddepsilon_kQ_dky_dky = -2.*t*cosky          + 4.*tp*coskx*cosky - 8.*tpp*cos2kx;
        //double ddepsilon_kQ_dkx_dky =                      - 4.*tp*sinkx*sinky;

        //precalculate sum, diff and radical:
        double Sk          = 0.5*(  epsilon_k         +   epsilon_kQ);
        double dSk_dkx     = 0.5*( depsilon_k_dkx     +  depsilon_kQ_dkx);
        //double dSk_dky     = 0.5*( depsilon_k_dky     +  depsilon_kQ_dky);
        //double ddSk_dkx_dkx= 0.5*(ddepsilon_k_dkx_dkx + ddepsilon_kQ_dkx_dkx);
        double ddSk_dky_dky= 0.5*(ddepsilon_k_dky_dky + ddepsilon_kQ_dky_dky);
        //double ddSk_dkx_dky= 0.5*(ddepsilon_k_dkx_dky + ddepsilon_kQ_dkx_dky);

        double Dk          = 0.5*(  epsilon_k         -   epsilon_kQ);
        double dDk_dkx     = 0.5*( depsilon_k_dkx     -  depsilon_kQ_dkx);
        double dDk_dky     = 0.5*( depsilon_k_dky     -  depsilon_kQ_dky);
        //double ddDk_dkx_dkx= 0.5*(ddepsilon_k_dkx_dkx - ddepsilon_kQ_dkx_dkx);
        double ddDk_dky_dky= 0.5*(ddepsilon_k_dky_dky - ddepsilon_kQ_dky_dky);
        //double ddDk_dkx_dky= 0.5*(ddepsilon_k_dkx_dky - ddepsilon_kQ_dkx_dky);

        double Rk          = sqrt( Dk*Dk + M*M );
        double dRk_dkx     = Dk*dDk_dkx/Rk;
        //double dRk_dky     = Dk*dDk_dky/Rk;
        //double ddRk_dkx_dkx= ((dDk_dkx*dDk_dkx + Dk*ddDk_dkx_dkx) - (Dk*Dk*dDk_dkx*dDk_dkx) / (Rk*Rk))/Rk;
        double ddRk_dky_dky= ((dDk_dky*dDk_dky + Dk*ddDk_dky_dky) - (Dk*Dk*dDk_dky*dDk_dky) / (Rk*Rk))/Rk;
        //double ddRk_dkx_dky= ((dDk_dkx*dDk_dky + Dk*ddDk_dkx_dky) - (Dk*Dk*dDk_dkx*dDk_dky) / (Rk*Rk))/Rk;

        //finally calculate the eigen values and their derivatives (vertices):
        double E1_k          =   Sk         +   Rk ;
        double dE1_k_dkx     =  dSk_dkx     +  dRk_dkx;
        double ddE1_k_dky_dky= ddSk_dky_dky + ddRk_dky_dky;
        
        double E2_k          =   Sk         -  Rk ;
        double dE2_k_dkx     =  dSk_dkx     -  dRk_dkx;
        double ddE2_k_dky_dky= ddSk_dky_dky - ddRk_dky_dky;
        
        double kernel1_xx = dE1_k_dkx*dE1_k_dkx;
        double kernel1_xy = dE1_k_dkx*dE1_k_dkx*ddE1_k_dky_dky;
        
        double kernel2_xx = dE2_k_dkx*dE2_k_dkx;
        double kernel2_xy = dE2_k_dkx*dE2_k_dkx*ddE2_k_dky_dky;

        
        for(int n=0; n<nOmega; n++)
        {
            complex<double> z(omega[n],ETA);
            
            double A1_k = -0.5*imag(1.0/ (z-E1_k) );
            double A2_k = -0.5*imag(1.0/ (z-E2_k) );

            double frequencyKernel1_xx = -dfermiDirac_dw[n]*kernel1_xx*A1_k*A1_k;
            double frequencyKernel2_xx = -dfermiDirac_dw[n]*kernel2_xx*A2_k*A2_k;
            
            double frequencyKernel1_xy = -dfermiDirac_dw[n]*kernel1_xy*A1_k*A1_k*A1_k;
            double frequencyKernel2_xy = -dfermiDirac_dw[n]*kernel2_xy*A2_k*A2_k*A2_k;
            
            
            sigma1_xx += frequencyKernel1_xx;
            sigma2_xx += frequencyKernel2_xx;
            
            sigma1_xy += frequencyKernel1_xy;
            sigma2_xy += frequencyKernel2_xy;
            
                        
            alpha1_xx += omega[n] * frequencyKernel1_xx;
            alpha2_xx += omega[n] * frequencyKernel2_xx;
            
            alpha1_xy += omega[n] * frequencyKernel1_xy;
            alpha2_xy += omega[n] * frequencyKernel2_xy;
            
            
            double omega2 = omega[n] * omega[n];
            beta1_xx += omega2 * frequencyKernel1_xx;
            beta2_xx += omega2 * frequencyKernel2_xx;
            
            beta1_xy += omega2 * frequencyKernel1_xy;
            beta2_xy += omega2 * frequencyKernel2_xy;
            
        }
        
        density += 1.0/(1.0+exp(beta*E1_k)) + 1.0/(1.0+exp(beta*E2_k));
    }


    double f0 = 1.0/nK/nK;
    double f = f0/nOmega;

    FILE *fileOut = fopen("conductivities.dat","w");
    printf("% 4.8f  \n", 1.0-f0*density);
    fprintf(fileOut,"% 4.8f  ", 1.0-f0*density);
    fprintf(fileOut,"% 4.8e % 4.8e % 4.8e % 4.8e ", f*sigma1_xx, f*sigma2_xx, f*sigma1_xy, f*sigma2_xy);
    fprintf(fileOut,"% 4.8e % 4.8e % 4.8e % 4.8e ", f*alpha1_xx, f*alpha2_xx, f*alpha1_xy, f*alpha2_xy);
    fprintf(fileOut,"% 4.8e % 4.8e % 4.8e % 4.8e ", f*beta1_xx,  f*beta2_xx,  f*beta1_xy,  f*beta2_xy);
    fprintf(fileOut,"\n");
    fclose(fileOut);
    
    printf("afmCond over.\n");
    return 0;
}
