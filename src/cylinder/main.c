/*                                    cre. on Feb. 12, 2004 by S.D.@SP8


*/
#include <stdio.h>  /* for printf,fgets NULL */
#include <stdlib.h>  /* for exit, malloc */
#include <ctype.h> /* for tolower */
#include <string.h> /* for strcpy */
#include <math.h>   /* to use sqrt */
#include <time.h>
#include "const.h" /* for R4OVER and so on */
/*===========================================================================*/
void main(){
char infile[80], outfile[80], logfile[80];
FILE *instr, *outstr, *logstr;
time_t time_pointer;

double At[10];
double *Mtrx1, *Mtrx2, *Mtrx3;
double *Ginv;
int ndata;

int status, n,i,j,k;
int n1, n2;


      strcpy(logfile,"main.log");
      logstr=fopen(logfile,"w");

n=3;
  Ginv = calloc(n*n,sizeof(double));
  Mtrx1 = calloc(n*n,sizeof(double));
  Mtrx2 = calloc(n*n,sizeof(double));
  Mtrx3 = calloc(n*n,sizeof(double));

/*
 for(i=0;i<n; i++){  
 for(j=0;j<n; j++){
    *(Mtrx1+n*i+j) = 10*i+j ;
    printf("[%d,%d]=%f\n",i,j,*(Mtrx1+n*i+j));
}}
*/

/*
n1=n;
n2=n;
   for(i=0;i<n1;i++){
     for(j=0;j<n2;j++){ printf("[%d][%d]=%f\n",i,j,*(Mtrx1+n1*i+j));
                     }
   }
*/

i=0;j=0;
*(Mtrx1+n*i+j) = 1.;
i=0;j=1;
*(Mtrx1+n*i+j) = 2.;
i=0;j=2;
*(Mtrx1+n*i+j) = 3.;
i=1;j=0;
*(Mtrx1+n*i+j) = -2.;
i=1;j=1;
*(Mtrx1+n*i+j) = .5;
i=1;j=2;
*(Mtrx1+n*i+j) = -3.;
i=2;j=0;
*(Mtrx1+n*i+j) = 3.;
i=2;j=1;
*(Mtrx1+n*i+j) = -.5;
i=2;j=2;
*(Mtrx1+n*i+j) = 1.;

i=0;j=0;
*(Mtrx2+n*i+j) = -.5;
i=0;j=1;
*(Mtrx2+n*i+j) = 3.;
i=0;j=2;
*(Mtrx2+n*i+j) = 5.;
i=1;j=0;
*(Mtrx2+n*i+j) = 6.;
i=1;j=1;
*(Mtrx2+n*i+j) = -2.;
i=1;j=2;
*(Mtrx2+n*i+j) = 1.;
i=2;j=0;
*(Mtrx2+n*i+j) = -2.;
i=2;j=1;
*(Mtrx2+n*i+j) = 4.;
i=2;j=2;
*(Mtrx2+n*i+j) = -1.;



Mprod(logstr, n, Mtrx1, Mtrx2, Mtrx3);



    dump(Mtrx3, n, n);

         free(Mtrx1);
         free(Mtrx2);
         free(Mtrx3);
}
/*===========================================================================*/
int Mprod(FILE *logstr,  int N,  double *A, double *B, double *AB){
/*  Product AB of rectangular matrices A and B
*/
int status, i, j, l;


 for(i=0;i<N; i++){  
 for(j=0;j<N; j++){
    *(AB + N*i +j) = 0.;
 for(l=0;l<N; l++){
    *(AB + N*i +j) = *(AB + N*i +j) + *(A + N*i +l) * *(B + N*l + j);
}}}


      return(0);
}
/*===========================================================================*/
int dump(double *A, int n1, int n2){
int i,j;

/* print */
/*   for(i=0;i<n1*n2;i++){ printf("aa[%d]=%d\n",i,aa[i]); } */

   for(i=0;i<n1;i++){
     for(j=0;j<n2;j++){ printf("[%d][%d]=%f\n",i,j,*(A+n1*i+j));
/*also works for(j=0;j<n2;j++){ printf("aa[%d][%d]=%d\n",i,j,aa[3*i+j]); */
                     }
   }
   
   return(0);
   
#if(0)
the followings don't work
   for(i=0;i<n1;i++){ 
   for(j=0;j<n2;j++){ printf("aa[%d][%d]=%d\n",i,j,aa[i][j]);
		     }
   }
#endif
}
/*===========================================================================*/
int NNptdist_ISR(FILE *logstr, double pt, double rs, double *invdifcr, char *kind){
/*
#  ISR, B. Alper et al Nucl. Phys. B100(1975)291-301
#    High Pt (>1GeV) fit based on  Table 4
#
#
#  E d^3 sigma/dp^3 = A(1-pt/pb)^m /(pt^2+M^2)^n
#
#
#
# From the top^[$B!"^[(B each term corresponds to
#  pi+, pi-, K+, K-, proton, anti-proton 
#
#  pb = pbeam = sqrt(s)/2
#
y1pp(x,pb) = 6.9*(1-x/pb)**(11.0) /(x*x+ 0.86*0.86)**(3.85)
y1pm(x,pb) = 7.4*(1-x/pb)**(11.9) /(x*x+ 0.89*0.89)**(3.89)
y1kp(x,pb) = 9.9*(1-x/pb)**(9.0) /(x*x+ 1.30*1.30)**(4.36)
y1km(x,pb) =10.4*(1-x/pb)**(12.2) /(x*x+ 1.33*1.33)**(4.38)
y1pr(x,pb) = 52*(1-x/pb)**(7.3) /(x*x+ 1.35*1.35)**(5.19)
y1ap(x,pb) = 9.0*(1-x/pb)**(14.0) /(x*x + 1.08*1.08)**(4.55)
*/

double pb, x;
double y1pp, y1pm, y1kp, y1km, y1pr, y1ap;

       pb = rs/2.;
       x = pt;

       *invdifcr = 0.;

   if(x > pb){ 
     /* fprintf(logstr,"E...NNptdist_ISR pt=%e > rs/2. =%e\n",pt,pb); */
      return(1);
    }

  if(strcmp(kind,"charged") == 0){

     *invdifcr = 6.9*pow((1-x/pb),11.0) /pow((x*x+ 0.86*0.86),3.85)
               + 7.4*pow((1-x/pb),11.9) /pow((x*x+ 0.89*0.89),3.89)
               + 9.9*pow((1-x/pb),9.0) /pow((x*x+ 1.30*1.30),4.36)
               + 10.4*pow((1-x/pb),12.2) /pow((x*x+ 1.33*1.33),4.38)
               + 52*pow((1-x/pb),7.3) /pow((x*x+ 1.35*1.35),5.19)
               + 9.0*pow((1-x/pb),14.0) /pow((x*x + 1.08*1.08),4.55);

     return(0);
   }
  else if(strcmp(kind,"pip") == 0){
     *invdifcr = 6.9*pow((1-x/pb),11.0) /pow((x*x+ 0.86*0.86),3.85);
     return(0);
   }
  else if(strcmp(kind,"pim") == 0){
     *invdifcr = 7.4*pow((1-x/pb),11.9) /pow((x*x+ 0.89*0.89),3.89);
     return(0);
   }
  else if(strcmp(kind,"kp") == 0){
     *invdifcr = 9.9*pow((1-x/pb),9.0) /pow((x*x+ 1.30*1.30),4.36);
     return(0);
   }
  else if(strcmp(kind,"km") == 0){
     *invdifcr = 10.4*pow((1-x/pb),12.2) /pow((x*x+ 1.33*1.33),4.38);
     return(0);
   }
  else if(strcmp(kind,"pr") == 0){
     *invdifcr = 52*pow((1-x/pb),7.3) /pow((x*x+ 1.35*1.35),5.19);
     return(0);
   }
  else if(strcmp(kind,"ap") == 0){
     *invdifcr = 9.0*pow((1-x/pb),14.0) /pow((x*x + 1.08*1.08),4.55);
     return(0);
   }

}
/*===========================================================================*/
int rapdens0(FILE *logstr,  double rshat, double pt, double *rho0){
/*   normalized rapidity density at y = 0 
    based on notes on 030906 meeting in Hiroshima

*/
double proms = 0.938271998;
double pims = 0.1395669;
double amt2, thymax, ymax;

     amt2 = pt*pt + pims*pims;
     thymax = sqrt(1. - 4.*amt2/rshat/rshat);
     ymax = log(rshat*(1. + thymax)/2./sqrt(amt2));

     if(ymax - thymax > 1.e-31){
          *rho0 = 1./2./(ymax - thymax);}
     else{ *rho0 = 0.;}

/*    fprintf(logstr,"%f %f %e\n",pt,rshat,(*rho0));*/

     return(0);
}
/*===========================================================================*/
int avchmul(FILE *logstr,  double sh, double *nch){
/*              Average Charged Multiplicity of the newly produced particles in NN inelastic collisions
            copied from URASiMA anti1.01 hhinel.f

        input:   shat   of  the secondary system in GeV^2
        output:   nch 

c
c       sh: CM energy ** 2 of secondary particles
*/

double c0=0.88, c1=0.44, c2=0.118, ain=-1.0, ajn=2.0;
double cb0, cb1, avm;
      

      cb0=c0-2.*c1*ain+c2*(6.*pow(ain,2.)-2.*ajn);
      cb1=c1-4.*c2*ain;
      *nch =cb0+cb1*log(sh)+c2*pow(log(sh),2.);  /* ! avarage multiplicity*/
 
     return(0);
}
/*===========================================================================*/
int eat_table(FILE *logstr, FILE *instr, double *arry1, int *nd){

int ndata, i, lskip;
char buffer[80], buffer2[80], *a;
double a1, a2, a3, a4;



/*---skip header lines ---*/
    lskip=1;
    i = 0;
    while( i < lskip ){
      fgets(buffer,sizeof(buffer),instr);
      i++;
    }

/*************** Data Line Loop *****************/
/*---loop of reading/printing---*/
    ndata = 0;
    while( ( fgets(buffer,sizeof(buffer),instr) )!=NULL ){

     i = strlen(buffer);
     if( i <= 1 ) continue;  /* skip blank lines */

        strcpy(buffer2, buffer);  /* to preserve buffer */

        a = strtok(buffer2," ");
        sscanf(a,"%lf",&a1);  
        a = strtok(NULL," ");
        sscanf(a,"%lf",&a2);  
        a = strtok(NULL," ");
        sscanf(a,"%lf",&a3);
        a = strtok(NULL," ");
        sscanf(a,"%lf",&a4);

        *(arry1 + ndata) = a2;

        ndata++;


      } /* data line loop */

     *nd = ndata;
     return(0);  
}








