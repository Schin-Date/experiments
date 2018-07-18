/*                                       created in Feb. 04-12, 2003 by S.D.@SP8
    V03 030312: nuclear profile, basic functions -> ws.o
    V02 030212: routines integrated in NN_coll_prob.

        This is a program to calculate various probabilities,
      averages, dispersions in AB collisions by Glauber theory.



*/
#include <stdio.h>  /* for printf,fgets NULL */
#include <stdlib.h>  /* for exit, malloc */
#include <ctype.h> /* for tolower */
#include <string.h> /* for strcpy */
#include <math.h>   /* to use sqrt */
#include <time.h>
#include "const.h" /* for R4OVER and so on */
/*===========================================================================*/
/* global variables */
    double sigmNN = 4.1; /* fm^2 : elementary NN cross section */
    double RAmax = 10.5; /* fm : maximum distance from the center of nucleus A */
    double RBmax = 10.5; /* fm : maximum distance from the center of nucleus B */
/*   For values of Rmax, see nuclear profile function in ws-ta-au.ps        */ 
/*--------------------------------------------------------------------------*/
void main(){
char logfile[80], infile[80], outfile[80];
FILE *logstr, *instr, *outstr;
time_t time_pointer;
int     seed = 566387;
double  getrandom();
void    setrandom(int);

int i, j, n, m, status;
double a, b, c, x, y, z;

char method[80], shape[20], approx[20];
int Ninteg;
double bimp, tA, sigmNNbI, PNBbI;
double sigmABbI_A, sigmABbI_B, sigmABbI_opt;
double Wnm[210][210];
double avNcbI;

      strcpy(logfile,"main.log");
      logstr=fopen(logfile,"w");
      strcpy(outfile,"main.out");
      outstr=fopen(outfile,"w");

      time(&time_pointer);
      fprintf(logstr,"present time = %.24s\n",ctime(&time_pointer));

      setrandom(seed);
      fprintf(logstr,"seed = %d\n",seed);

       bimp = 13.4;
       bimp = 0.;
       strcpy(shape,"disk");
       strcpy(shape,"WS");  /* Woods-Saxon */
       strcpy(shape,"WS_tbl"); /* Woods-Saxon table interpolation */
       strcpy(method,"mont");
       strcpy(method,"mesh");
       Ninteg = 10000;
       Ninteg = 1000;
       n=197; m=197;
       fprintf(logstr,"%s %d %s bI=%.1f fm\n", method, Ninteg, shape, bimp);
       printf("%s %d %s bI=%.1f fm\n", method, Ninteg, shape, bimp);


/*--- read tA(b) table ---*/
      status =  eat_tbl_profile(logstr, bimp, &tA, "read");
            if(status != 0){  
             fprintf(logstr,"E...main< eat_tbl_profile"
                 "status = %d\n",status);
             exit(1); 
            }


/*--- compute the main subject ---*/
       status = NN_coll_prob(logstr, "Au", "Au", bimp, n, m,
                &sigmNNbI, &PNBbI, Wnm, method, Ninteg, shape);
            if(status != 0){  
             fprintf(logstr,"E...main< NN_coll_prob"
                 "status = %d\n",status);
             exit(1); 
            }
 
       fprintf(logstr,"N...sigmNNbI = %e PNBbI = %e\n",sigmNNbI, PNBbI);
       printf("N...sigmNNbI = %e PNBbI = %e\n",sigmNNbI, PNBbI);

       status = AB_coll_prob(logstr, "Au", "Au", sigmNNbI, PNBbI,
                 &sigmABbI_A, &sigmABbI_B, &sigmABbI_opt );
            if(status != 0){  
             fprintf(logstr,"E...main< AB_coll_prob:"
                 "status = %d\n",status);
             exit(1); 
            }

       avNcbI = 197.*197. * sigmNNbI / sigmABbI_A;

       fprintf(logstr,"N...sigmABbI A = %e  B = %e opt = %e <Nc>(A) = %e\n", 
                     sigmABbI_A, sigmABbI_B, sigmABbI_opt, avNcbI);
       fprintf(logstr,"output table:\n"
                   " <Nc>(A)  A  B   bI  [1st line]\n");
       fprintf(outstr,"%e 197. 197. %.3f\n", avNcbI, bimp);

       fprintf(logstr,"n m P(n,m) sum [only P> 1./4.10^6 are printed]\n");
     z = 0.;
        for(i=0;i<n+1;i++){
        for(j=0;j<m+1;j++){

        x = *(*Wnm+(m+1)*i+j);
        y = x / 197./197./sigmNNbI;
    if(i>0 & j>0){        z = z + y;   }
        if( y > 1./4.e6 ){
          fprintf(outstr,"%d %d %e %e\n",i,j,y,z);
	}
    }}

        time(&time_pointer);
        fprintf(logstr,"present time = %.24s\n",ctime(&time_pointer));
}
/*===========================================================================*/
int AB_coll_prob(FILE *logstr, char *A, char *B, double sigmNNbI, double PNBbI, 
                 double *sigmABbI_A, double *sigmABbI_B, double *sigmABbI_opt){
 
int IAmass, IBmass;

     if(strcmp(A,"Au") == 0 & strcmp(B,"Au") == 0){
        IAmass = 197;  IBmass = 197;
     }     
         else{ 
         printf("E...AB_coll_prob: A:'%s' B:'%s'\n",A,B);
         fprintf(logstr,"E...AB_coll_prob: A:'%s' B:'%s'\n",A,B);
         return(1);
         }                 

/* Approx. A */
       *sigmABbI_A = 1. - pow((1. - sigmNNbI),(double)(IAmass*IBmass));

/* Approx. B */
       *sigmABbI_B = 1. - pow(PNBbI,(double)IAmass);

/* Approx. opt */
       *sigmABbI_opt = 1. - exp(- sigmNNbI*(double)(IAmass*IBmass));

       return(0);
}
/*===========================================================================*/
int NN_coll_prob(FILE *logstr, char *A, char *B, double bI, int n, int m, 
       double *sigmNNbI, double *PNBbI, double *Wnm,
                 char *method, int Ninteg, char *shape){

     if(strcmp(method,"mesh") == 0){
          NN_coll_prob_mesh(logstr, A, B, bI, n, m, 
                            sigmNNbI, PNBbI, Wnm, Ninteg, shape);
         }
                                /* Ninteg= 10^4 heavy on kanbai */

     else if(strcmp(method,"mont") == 0){
          NN_coll_prob_mont(logstr, A, B, bI, n, m, 
                            sigmNNbI, PNBbI, Wnm, Ninteg, shape);
         }

     else{
         printf("E...NN_coll_prob: method:'%s'\n",method);
         fprintf(logstr,"E...NN_coll_prob: method:'%s'\n",method);
         return(1);
         }                 

     return(0);
}
/*---------------------------------------------------------------------------*/
int NN_coll_prob_mont(FILE *logstr, char *A, char *B, double bI, int n, int m, 
      double *sigmNNbI, double *PNBbI, double *Wnm, int Ninteg, char *shape){

double  getrandom();
int IAmass, IBmass;
int i,j,in,im,status;
double sx, sy, sxmin, sxmax, symax,  sA, sB, wAm, wBn;
double x,y,z, sumsNN, sumPN, *wAm_arr, *sum;
double tA, tB;

        *sigmNNbI = 0.;
        *PNBbI = 1.;

        for(i=0;i<n+1;i++){
        for(j=0;j<m+1;j++){
               *(Wnm + (m+1)*i + j)    = 0.;
    }}

         if( bI > RAmax + RBmax ){ return(0); }

         if(strcmp(A,"Au") == 0 & strcmp(B,"Au") == 0){
            IAmass = 197;  IBmass = 197;
         }     
             else{ 
             printf("E...NN_coll_prob_mont: A:'%s' B:'%s'\n",A,B);
             fprintf(logstr,"E...NN_coll_prob_mont: A:'%s' B:'%s'\n",A,B);
             return(1);
             }                 

    /* Assume the vector bI directs in the x-direction */
         if( bI==0. || bI <= fabs(RAmax - RBmax) ){
             if(RAmax - RBmax < 0.){
                sxmin = -RAmax; sxmax = RAmax;
                symax = RAmax;
	      }
              else{
                sxmin = -bI -RBmax; sxmax = -bI + RBmax;
                symax = RBmax;
	      }
         }          
         else{
               sxmin = -RAmax;  sxmax = -bI + RBmax;
               sx = fabs(RBmax*RBmax - RAmax*RAmax - bI*bI)/2./bI;
               symax =  sqrt(RAmax*RAmax - sx*sx);
	 }

         sumsNN =0.;
         sumPN = 0.;
         sum = calloc((n+1)*(m+1),sizeof(double));
         wAm_arr = calloc(m+1,sizeof(double));

      for(i = 0; i <Ninteg; i++){
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ integration loop */
   /* generation of uniform random */
        sx = sxmin + getrandom() * (sxmax - sxmin);
        sy = getrandom() * symax;

         sA = sqrt(sx*sx + sy*sy);
         sB = sqrt((bI+sx)*(bI+sx) + sy*sy);
         status = nuclear_profile(logstr, A, sA, &tA, shape);
              if(status != 0){  
               fprintf(logstr,"E...NN_coll_prob_mont A< nuclear_profile: "
                   "status = %d\n",status);
               fprintf(logstr," ...sA=%e sx=%e sxmin=%e sxmax=%e\n"
                              "          sy=%e symax=%e\n",
                    sA,sx,sxmin,sxmax,sy,symax);
               exit(1); 
              }
         status = nuclear_profile(logstr,B , sB, &tB, shape);
              if(status != 0){  
               fprintf(logstr,"E...NN_coll_prob_mont B< nuclear_profile:"
                   " status = %d\n",status);
               fprintf(logstr," ...sB=%e sx=%e sxmin=%e sxmax=%e\n"
                              "    bI=%e sy=%e symax=%e\n",
                    sA,sx,sxmin,sxmax,bI,sy,symax);
               exit(1); 
              }

     if( tA*tB != 0 ){
             
        wAm = 0.; wBn = 0.; 

        x = (double)(IAmass-1) * log(sigmNN * tA);
        if( x > -32. ){ wAm = exp(x); }

        x = (double)(IBmass-1) * log(sigmNN * tB);
        if( x > -32. ){ wBn = exp(x); }


        for(im=IAmass - 1; im>m; im--){
              nikou(logstr, IAmass, im, &x, &y);
             y = x + (double)(im-1)*log(sigmNN * tA)  
                   + (double)(IAmass-im)*log(1.-sigmNN * tA);
            if( y > -32. ){ wAm = wAm + exp(y); }
	}
        for(in=IBmass - 1; in>n; in--){
              nikou(logstr, IBmass, in, &x, &y);
             y = x + (double)(in-1)*log(sigmNN * tB) 
                  + (double)(IBmass-in)*log(1.-sigmNN * tB);
            if( y > -32. ){ wBn = wBn + exp(y); }
	}

        for(im=m; im>-1; im--){
              nikou(logstr, IAmass, im, &x, &y);
             y = x + (double)(im-1)*log(sigmNN * tA) 
                  + (double)(IAmass-im)*log(1.-sigmNN * tA);
             if( y > -32. ){ 
               wAm = wAm + exp(y); 
               *(wAm_arr + im) = wAm;
             }
	 }


        for(in=n; in>-1; in--){
              nikou(logstr, IBmass, in, &x, &y);
             y = x + (double)(in-1)*log(sigmNN * tB) 
                  + (double)(IBmass-in)*log(1.-sigmNN * tB);
             if( y > -32. ){ wBn = wBn + exp(y); }

        for(im=m; im>-1; im--){
             wAm = *(wAm_arr + im);

        *(sum+ (m+1)*in + im) 
           = *(sum+ (m+1)*in + im) +  tA * tB * wBn * wAm;

         }
         }

         sumsNN =sumsNN + tA*tB;
         sumPN = sumPN + tA*(1. - pow(1.- sigmNN * tB,(double)IBmass));
         } /* tA* tB != 0 */


      } /* end of integration loop */


        for(in=0; in<n+1; in++){
        for(im=0; im<m+1; im++){

     *(Wnm + (m+1)*in + im)    = 
        sigmNN * 2.* (sxmax-sxmin)* symax * (*(sum+ (m+1)*in + im))/ (double)Ninteg;
         }
         }

     *sigmNNbI = sigmNN * 2.* (sxmax-sxmin)* symax * sumsNN / (double)Ninteg;
     *PNBbI = 1. - 2.* (sxmax-sxmin)* symax * sumPN / (double)Ninteg;


     free(sum);
     free(wAm_arr);
     return(0);

}      
/*---------------------------------------------------------------------------*/
int NN_coll_prob_mesh(FILE *logstr, char *A, char *B, double bI, int n, int m, 
      double *sigmNNbI, double *PNBbI, double *Wnm, int Ninteg, char *shape){

int IAmass, IBmass;
int i,j,in,im,status;
double sx, sy, sxmin, sxmax, symax,  sA, sB, sborder, wAm, wBn;
double Rupp, Rdwn, dsx, dsy;
double x,y,z, sumsNN, sumPN, *wAm_arr, *sum;
double sumysNN, sumyPN, *sumy;
double tA, tB;

        *sigmNNbI = 0.;
        *PNBbI = 1.;

        for(i=0;i<n+1;i++){
        for(j=0;j<m+1;j++){
               *(Wnm + (m+1)*i + j)    = 0.;
    }}

         if( bI > RAmax + RBmax ){ return(0); }
/*     .........................................  */
         if(strcmp(A,"Au") == 0 & strcmp(B,"Au") == 0){
            IAmass = 197;  IBmass = 197;
         }     
             else{ 
             printf("E...NN_coll_prob_mont: A:'%s' B:'%s'\n",A,B);
             fprintf(logstr,"E...NN_coll_prob_mont: A:'%s' B:'%s'\n",A,B);
             return(1);
             }                 

    /* Assume the vector bI directs in the x-direction */
         if( bI==0. || bI <= fabs(RAmax - RBmax) ){
             if(RAmax - RBmax < 0.){
                sxmin = -RAmax; sxmax = RAmax;
                sborder = RAmax;
	      }
              else{
                sxmin = -bI -RBmax; sxmax = -bI + RBmax;
                sborder = - RAmax;
	      }
         }          
         else{
               sxmin = -RAmax;  sxmax = -bI + RBmax;
               sborder = (RBmax*RBmax - RAmax*RAmax - bI*bI)/2./bI;
	 }
         Rupp = RBmax; Rdwn = RAmax; 

         sumsNN =0.;
         sumPN = 0.;
         sum = calloc((n+1)*(m+1),sizeof(double));
         sumy = calloc((n+1)*(m+1),sizeof(double));
         wAm_arr = calloc(m+1,sizeof(double));

      dsx = (sxmax - sxmin)/(double)Ninteg;
      for(i = 0; i <Ninteg; i++){
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ sx-integration loop */
         sx = sxmin + dsx/2. + (double)i*dsx;
         if( sx > sborder ){
               symax =  sqrt(RBmax*RBmax - (bI+sx)*(bI+sx));
         }       
         else{
               symax =  sqrt(RAmax*RAmax - sx*sx);
        }
         dsy = symax/(double)Ninteg;

         sumysNN =0.;
         sumyPN = 0.;
         for(in=0; in<n+1; in++){
         for(im=0; im<m+1; im++){
          *(sumy + (m+1)*in + im)    = 0.;
          }
          }
     for(j=0;j<Ninteg;j++){
/*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ sy-integration loop */
         sy =  dsy/2. + (double)j*dsy;
         sA = sqrt(sx*sx + sy*sy);
         sB = sqrt((bI+sx)*(bI+sx) + sy*sy);
         status = nuclear_profile(logstr, A, sA, &tA, shape);
              if(status != 0){  
               fprintf(logstr,"E...NN_coll_prob_mesh A< nuclear_profile:"
                   " status = %d\n",status);
               fprintf(logstr," ...sA=%e sx=%e sxmin=%e sxmax=%e\n"
                              "          sy=%e symax=%e\n",
                    sA,sx,sxmin,sxmax,sy,symax);
               exit(1); 
              }
         status = nuclear_profile(logstr,B , sB, &tB, shape);
              if(status != 0){  
               fprintf(logstr,"E...NN_coll_prob_mesh B< nuclear_profile:"
                   " status = %d\n",status);
               fprintf(logstr," ...sB=%e sx=%e sxmin=%e sxmax=%e\n"
                              "    bI=%e sy=%e symax=%e\n",
                    sA,sx,sxmin,sxmax,bI,sy,symax);
               exit(1); 
              }
     if( tA*tB != 0 ){

        wAm = 0.; wBn = 0.; 

        x = (double)(IAmass-1) * log(sigmNN * tA);
        if( x > -32. ){ wAm = exp(x); }

        x = (double)(IBmass-1) * log(sigmNN * tB);
        if( x > -32. ){ wBn = exp(x); }

        for(im=IAmass - 1; im>m; im--){
              nikou(logstr, IAmass, im, &x, &y);
             y = x + (double)(im-1)*log(sigmNN * tA)  
                   + (double)(IAmass-im)*log(1.-sigmNN * tA);
            if( y > -32. ){ wAm = wAm + exp(y); }
	}
        for(in=IBmass - 1; in>n; in--){
              nikou(logstr, IBmass, in, &x, &y);
             y = x + (double)(in-1)*log(sigmNN * tB) 
                  + (double)(IBmass-in)*log(1.-sigmNN * tB);
            if( y > -32. ){ wBn = wBn + exp(y); }
	}

        for(im=m; im>-1; im--){
              nikou(logstr, IAmass, im, &x, &y);
             y = x + (double)(im-1)*log(sigmNN * tA) 
                  + (double)(IAmass-im)*log(1.-sigmNN * tA);
             if( y > -32. ){ 
               wAm = wAm + exp(y); 
               *(wAm_arr + im) = wAm;
             }
	 }

        for(in=n; in>-1; in--){
              nikou(logstr, IBmass, in, &x, &y);
             y = x + (double)(in-1)*log(sigmNN * tB) 
                  + (double)(IBmass-in)*log(1.-sigmNN * tB);
             if( y > -32. ){ wBn = wBn + exp(y); }

        for(im=m; im>-1; im--){
             wAm = *(wAm_arr + im);

        *(sumy+ (m+1)*in + im) 
           = *(sumy+ (m+1)*in + im) +  dsy * tA * tB * wBn * wAm;

         }
         }

         sumysNN = sumysNN + dsy*tA*tB;
         sumyPN = sumyPN + dsy*tA*(1. - pow(1.- sigmNN * tB,(double)IBmass));

          } /* tA*tB != 0 */


       } /* y-integration */

         sumsNN = sumsNN + dsx * 2.* sumysNN;
         sumPN = sumPN + dsx * 2.* sumyPN;

        for(in=0; in<n+1; in++){
        for(im=0; im<m+1; im++){
        *(sum+ (m+1)*in + im) 
           = *(sum+ (m+1)*in + im) +  dsx * 2.* (*(sumy+ (m+1)*in + im));
      }
      }


     } /* x-integration */


     *sigmNNbI = sigmNN * sumsNN;
     *PNBbI =  1. - sumPN;

        for(in=0; in<n+1; in++){
        for(im=0; im<m+1; im++){
     *(Wnm + (m+1)*in + im)    = 
        sigmNN * (*(sum+ (m+1)*in + im));
         }
         }


     free(sum);
     free(sumy);
     free(wAm_arr);
     return(0);

}
/*===========================================================================*/
