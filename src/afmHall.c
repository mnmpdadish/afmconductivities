//
//  afmCond.c
//

#include "utilities.h"


// 2x2 matrix indices:
//
// 0 1
// 2 3

// A=B
void matrix2x2Copy(double A[], double const B[]){
    int n;
    for(n=0;n<4;n++) A[n] = B[n];
}

void matrix2x2Scale(double A[], double const factor){
    int n;
    for(n=0;n<4;n++) A[n] *= factor;
}


// A = B*C
void two_matrix2x2Multiplication(double A[], double const B[], double const C[]){
    A[0] = B[0]*C[0] + B[1]*C[2];
    A[1] = B[0]*C[1] + B[1]*C[3];
    A[2] = B[2]*C[0] + B[3]*C[2];
    A[3] = B[2]*C[1] + B[3]*C[3];
}

// A = B*C*D
void three_matrix2x2Multiplication(double A[], double const B[], double const C[], double const D[]){
    double dummy1[4];
    two_matrix2x2Multiplication(dummy1, C, D);
    two_matrix2x2Multiplication(A, B, dummy1);
}

// A = B*C*D*E
void four_matrix2x2Multiplication(double A[], double const B[], double const C[], double const D[], double const E[]){
    double dummy1[4];
    double dummy2[4];
    two_matrix2x2Multiplication(dummy1, B, C);
    two_matrix2x2Multiplication(dummy2, D, E);
    two_matrix2x2Multiplication(A, dummy1, dummy2);
}



// A = B^-1
void matrix2x2Inverse(double * A, double const * B){
    double det_inverse = 1.0/(B[0]*B[3] - B[1]*B[2]);
    A[0] =  B[3]*det_inverse;
    A[1] = -B[1]*det_inverse;
    A[2] = -B[2]*det_inverse;
    A[3] =  B[0]*det_inverse;
}




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
    //double energyCutoff = 2.*2.*acosh(0.25*sqrt(beta/amplitudeCutoff)) /beta; // we put an additionnal factor of 2. for the derivative cutoff
    
    //double omega[nOmega];
    //double fermiDirac_dw[nOmega];
    //double dfermiDirac_dw[nOmega];
    
    /*
    int n=0; for(n=0; n<nOmega; n++)
    {
        omega[n]= -energyCutoff + 2.*energyCutoff*n/(nOmega-1);
        double expBw = exp(beta*omega[n]);
        dfermiDirac_dw[n]= -beta*expBw/((expBw+1.)*(expBw+1.));
    }
    */
    
    
    FILE *fileOut = fopen("conductivities.dat","w");
    
    printf("#mu          density      sigma_xx        sigma_xy        sigma_xx_bubble ");

    fprintf(fileOut, "#mu          density      ");
    fprintf(fileOut, "sigma1_xx       sigma_xy        ");
    fprintf(fileOut, "sigma_xx_bubble ");
    fprintf(fileOut, "sigma_xy_tri1   sigma_xy_tri2   sigma_xy_rect1  sigma_xy_rect2  \n");
    
    
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
       double sigma_xx = 0., sigma_xy = 0., sigma_xx_bubble = 0.; 
       double sigma_xy_tri1 = 0.,   sigma_xy_tri2 = 0., sigma_xy_rect1 = 0., sigma_xy_rect2 = 0.;
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
              
              double zeta_k           = -4.*tp*cosk[i]*cosk[j] - 2.*tpp*(cos2k[i] + cos2k[j]) - mu;
              double dzeta_k_dkx      =  4.*tp*sink[i]*cosk[j] + 4.*tpp*sin2k[i];
              double dzeta_k_dky      =  4.*tp*cosk[i]*sink[j] + 4.*tpp*sin2k[j];
              double ddzeta_k_dky_dky =  4.*tp*cosk[i]*cosk[j] + 8.*tpp*cos2k[j];
              double ddzeta_k_dkx_dky = -4.*tp*sink[i]*sink[j];
              
              double xi_k           = -2.*t*(cosk[i] + cosk[j]);
              double dxi_k_dkx      =  2.*t*sink[i];
              double dxi_k_dky      =  2.*t*sink[j];
              double ddxi_k_dky_dky =  2.*t*cosk[j];
              double ddxi_k_dkx_dky =  0.0;
              
              //defining multiple 2x2 matrices:
              double lambda_x[4], lambda_y[4], lambda_yy[4], lambda_xy[4], A[4], t_plus_self[4], to_invert[4];
              
              t_plus_self[0] = zeta_k + M;   t_plus_self[1] = xi_k;
              t_plus_self[2] = xi_k;         t_plus_self[3] = zeta_k - M;
              
              lambda_x[0] = dzeta_k_dkx;  lambda_x[1] = dxi_k_dkx;
              lambda_x[2] = dxi_k_dkx;    lambda_x[3] = dzeta_k_dkx;
              
              lambda_y[0] = dzeta_k_dky;  lambda_y[1] = dxi_k_dky;
              lambda_y[2] = dxi_k_dky;    lambda_y[3] = dzeta_k_dky;
              
              lambda_yy[0] = ddzeta_k_dky_dky;  lambda_yy[1] = ddxi_k_dky_dky;
              lambda_yy[2] = ddxi_k_dky_dky;    lambda_yy[3] = ddzeta_k_dky_dky;

              lambda_xy[0] = ddzeta_k_dkx_dky;  lambda_xy[1] = ddxi_k_dkx_dky;
              lambda_xy[2] = ddxi_k_dkx_dky;    lambda_xy[3] = ddzeta_k_dkx_dky;
              
              // dummy1 = omega-t
              //matrix2x2Copy(dummy1,t);
              //matrix2x2Scale(dummy1,1.0);
              //dummy1[0]+=mu; // +omega=0;
              //dummy1[3]+=mu; // +omega=0;
              
              // dummy2 = (omega-t)^2 +eta^2
              two_matrix2x2Multiplication(to_invert, t_plus_self, t_plus_self);
              to_invert[0]+=ETA*ETA;
              to_invert[3]+=ETA*ETA;
              
              // A= 1./((omega-t)^2 +eta^2)
              matrix2x2Inverse(A, to_invert);
              matrix2x2Scale(A,ETA/M_PI);
              //bref,...
              
              //Calculating Reza's 2x2 matrices and stuff:
              double Ax[4], Ay[4], Axy[4], Ayy[4], tAx[4], tAy[4];
              two_matrix2x2Multiplication(Ax,  A, lambda_x);
              two_matrix2x2Multiplication(Ay,  A, lambda_y);
              two_matrix2x2Multiplication(Ayy, A, lambda_yy);
              two_matrix2x2Multiplication(Axy, A, lambda_xy);
              two_matrix2x2Multiplication(tAx, t_plus_self, Ax);
              two_matrix2x2Multiplication(tAy, t_plus_self, Ay);
              
              double bubble[4], tri1[4], tri2[4], rect1[4], rect2[4], rect3[4], rect4[4];
              two_matrix2x2Multiplication(bubble, Ax, Ax);
              three_matrix2x2Multiplication(tri1, Ax, Ay, Axy);
              three_matrix2x2Multiplication(tri2, Ax, Ax, Ayy);
              four_matrix2x2Multiplication(rect1, tAx, Ay, Ax, Ay);
              four_matrix2x2Multiplication(rect2, tAy, Ax, Ax, Ay);
              four_matrix2x2Multiplication(rect3, tAx, Ay, Ay, Ax);
              four_matrix2x2Multiplication(rect4, tAy, Ax, Ay, Ax);
              
              sigma_xx_bubble+= bubble[0]+ bubble[3]; //trace
              sigma_xy_tri1  += tri1[0]  + tri1[3]; //trace
              sigma_xy_tri2  += tri2[0]  + tri2[3]; //trace
              sigma_xy_rect1 += rect1[0] + rect1[3] + rect4[0] + rect4[3] - (rect2[0] + rect2[3] + rect3[0] + rect3[3]);
              sigma_xy_rect2 += M*rect1[0] - M*rect1[3] - M*rect2[0] + M*rect2[3];
              //sigma_xy_rect1 += M*rect1[0] - M*rect1[3]; //not exactly trace, still trying to understand (but it works)...
              //sigma_xy_rect2 += M*rect2[0] - M*rect2[3]; //same
              
              
              

              //other approach:
              //calculation toward eigenvalues basis:
              double epsilon_k           = xi_k           + zeta_k;
              double depsilon_k_dkx      = dxi_k_dkx      + dzeta_k_dkx;
              double depsilon_k_dky      = dxi_k_dky      + dzeta_k_dky;
              double ddepsilon_k_dky_dky = ddxi_k_dky_dky + ddzeta_k_dky_dky;
              double ddepsilon_k_dkx_dky = ddxi_k_dkx_dky + ddzeta_k_dkx_dky;

              double epsilon_kQ           = -xi_k           + zeta_k;
              double depsilon_kQ_dkx      = -dxi_k_dkx      + dzeta_k_dkx;
              double depsilon_kQ_dky      = -dxi_k_dky      + dzeta_k_dky;
              double ddepsilon_kQ_dky_dky = -ddxi_k_dky_dky + ddzeta_k_dky_dky;
              double ddepsilon_kQ_dkx_dky = -ddxi_k_dkx_dky + ddzeta_k_dkx_dky;


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
              double kernel1_xy = -(2./3.)*(dE1_k_dkx*(dE1_k_dkx*ddE1_k_dky_dky - dE1_k_dky*ddE1_k_dkx_dky));
              double kernel2_xy = -(2./3.)*(dE2_k_dkx*(dE2_k_dkx*ddE2_k_dky_dky - dE2_k_dky*ddE2_k_dkx_dky));


              ///////////////////////
              double complex z = 0.0 + ETA * I; // omega=0.0;
                  
              double A1_k = -(1.0/M_PI)*cimag(1.0/ (z-E1_k) );
              double A2_k = -(1.0/M_PI)*cimag(1.0/ (z-E2_k) );
              
              sigma_xx += kernel1_xx*A1_k*A1_k      + kernel2_xx*A2_k*A2_k;
              sigma_xy += kernel1_xy*A1_k*A1_k*A1_k + kernel2_xy*A2_k*A2_k*A2_k;
              
              density += 1.0/(1.0+exp(beta*E1_k)) + 1.0/(1.0+exp(beta*E2_k));
          }
       }
       double f0 = 1.0/nK/nK;
       double c1 = 1.0;
       double c2 = 1.0;
       double c3 = M_PI/ETA;

       printf("\n% 4.8f % 4.8f  % 4.8e % 4.8e % 4.8e", mu, 1.0-f0*density,  f0*sigma_xx*c1, f0*sigma_xy*c2, f0*sigma_xx_bubble*c1); 
       fprintf(fileOut,"% 4.8f % 4.8f  ", mu, 1.0-f0*density);
       fprintf(fileOut,"% 4.8e % 4.8e ", f0*sigma_xx*c1, f0*sigma_xy*c2);
       fprintf(fileOut,"% 4.8e ", f0*sigma_xx_bubble*c1);
       fprintf(fileOut,"% 4.8e % 4.8e % 4.8e % 4.8e ", f0*sigma_xy_tri1*c2, f0*sigma_xy_tri2*c2, f0*sigma_xy_rect1*c3, f0*sigma_xy_rect2*c3); //6,7,8,9,10
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
