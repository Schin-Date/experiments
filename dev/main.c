                                  /*   created on Jun. 05, 1997 by S.D. @SP8 
29Jan01 revising
* moon_pos combined with main.
29Jan01 reviewing
06Jan99 revised

*/
#include <stdio.h>  /* for printf,fgets NULL */
#include <stdlib.h>  /* for exit, malloc */
#include <ctype.h> /* for tolower */
#include <string.h> /* for strcpy */
#include <math.h>   /* to use sqrt */
#include <time.h>
#include "/users/sr/schin/clib/head/const.h" /* for R4OVER and so on */
/*===========================================================================*/
void main(){
double year,month,day,hrs,mins,sec,ut,UT,days_97,tstep,abmin,abmin0;
char date[]="yyyy/mm/dd hh:mm:ss",date0[]="yyyy/mm/dd hh:mm:ss";
char dend[]="yyyy/mm/dd hh:mm:ss", dstrt[]="yyyy/mm/dd hh:mm:ss";
double xdays;

FILE *outstr;
char outfile[30];
char buffer[80];

double julian, stime_rad, g0dist, gra_rad, gdec_rad;
double tra_rad, tdec_rad, tdist; 

double Sra_rad, Sdec_rad;
double e_sun[3], e_moon[3];

double q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,a=0,b=0,c=0,d=0;

double rglov=6.37814;/*Mm*/
double r0sun=149597.870; /*Mm*/
double rrsun=1.016; /* rsun/r0sun    VALUE for June 20 - July 20*/
double r0mrgl=60.2682; /* r0moon/rglov */
double U0sun=8.87210e8; /* m^2/s^2 */
double U0moon=1.27554e4; /* m^2/s^2 */

double rads; 
rads = PAI/180.;

/* set defaults for date and location */
/*      strcpy(date, "1998/02/18 00:00:00");*/ /*!!!starting time in DT + UT !!!*/
/*      strcpy(dstrt, "1998/06/17 00:00:00");*/ /*!!!starting time in DT + UT !!!*/
/*      strcpy(dend, "1998/07/03 00:00:00");*/ /*!!!end time in DT + UT !!!*/
/*      strcpy(date0, "1998/06/22 00:00:00");*/ /*!!!origin time in DT + UT !!!*/

      strcpy(dstrt, "2000/05/03 00:00:00"); /*!!!starting time in DT + UT !!!*/
      strcpy(dend, "2000/05/12 00:00:00"); /*!!!end time in DT + UT !!!*/
      strcpy(date0, "2000/05/06 00:00:00"); /*!!!origin time in DT + UT !!!*/
      tstep=0.5; /* time step in hours */
      UT = +9.;
      abmin0 = 2880.;  /* definition of time=0 in min */
      strcpy(outfile,"98_08-cel");
      strcpy(outfile,"test.out");

      printf("\n\ncalculation starts from: ");
      printf("%s",dstrt);
      sprintf(buffer,"%s",dstrt);
      inp_edt(buffer);
      if( (strlen(buffer)-strlen(dstrt)) == 1 &&
           buffer[strlen(dstrt)]=='q' ) exit(0);
      strcpy(dstrt,buffer);
      strcpy(date,dstrt);

      printf("    calculation ends at: ");
      printf("%s",dend);
      sprintf(buffer,"%s",dend);
      inp_edt(buffer);
      strcpy(dend,buffer);

      printf("        origin of date: ");
      printf("%s",date0);
      sprintf(buffer,"%s",date0);
      inp_edt(buffer);
      strcpy(date0,buffer);

      printf("time step in hours = ");
      printf(        "%f",tstep);
      sprintf(buffer,"%f",tstep);
      inp_edt(buffer);
      sscanf(buffer,"%lf",&tstep);

      days_from(dstrt,date0,&xdays); abmin0=xdays*24.*60.;

      printf("output file name: ");
      printf("%s",outfile);
      sprintf(buffer,"%s",outfile);
      inp_edt(buffer);
      strcpy(outfile,buffer);

/* resolve the date into numbers */
      if( sscanf(dstrt,"%lf/%lf/%lf %lf:%lf:%lf",
            &year,&month,&day,&hrs,&mins,&sec) != 6 ) exit(1);

/*** reset for the time loop ***/
      strcpy(date,dstrt);
      days_from(dend,date,&xdays); 
      abmin = 0.;

/*** confirm the setting ***/
      printf("\n");
      printf("outfile=\"%s\"\n",outfile);
      printf("calculation for %d days from %s to %s\n",(int)(xdays+1.),dstrt, dend);
      printf("origin of time = %s (abmin0 = %f min)\n",date0,abmin0);
      printf("time step = %f hours\n",tstep);

      printf("go");
      sprintf(buffer,"go");
      inp_edt(buffer);
      if(strstr(buffer,"q")!= NULL) exit(0);

      outstr=fopen(outfile,"w");

/***** the time loop *****/
  while(xdays >= 0){

    /* days from the begining of 1997 */
      days_from(date,"1997/01/01 00:00:00",&x);
      days_97 = (double)(int)(x+1.);      
    /* resolve the date into numbers */
      if( sscanf(date,"%lf/%lf/%lf %lf:%lf:%lf",
            &year,&month,&day,&hrs,&mins,&sec) != 6 ) exit(1);
    /* evaluate ut <-- !!! give DT for this variable !!!*/
      ut= hrs - UT + mins/60. + sec/3600.;          
      if( ut < 0 ){ days_97 = days_97 - 1.; ut = ut + 24.; }

    /*--- calculate moon position ---*/       
    /* moon position at Greenwich on equator */
      if( moon_pos(0., 0., days_97, ut
              ,&julian, &stime_rad, &g0dist, &gra_rad, &gdec_rad
              , &tra_rad, &tdec_rad, &tdist) != 0 ){
        printf("E...moon_pos!=0 for days=%f, ut=%f, glong=0, glat=0\n\n"
                ,days_97, ut);
        exit(1);   }  

    /*--- calculate the sun position ---*/       
      if( sun_pos(days_97, ut, &Sra_rad, &Sdec_rad) !=0 ){
        printf("E...sun_pos!=0 for days=%f, ut=%f\n\n"
                ,days_97, ut);
        exit(1);   }  

    /*--- convert to the local vectors ---*/
        G2SP8(Sra_rad, Sdec_rad, ut, days_97, e_sun );
        G2SP8(gra_rad, gdec_rad, ut, days_97, e_moon );

        x = e_moon[0]; y = e_moon[1]; z = e_moon[2]; 
        legendr(2,z,&q); legendr(3,z,&r);
        u = pow( r0mrgl/g0dist, 3.) * q
          + pow( r0mrgl/g0dist, 4.) * r / r0mrgl;
        u = u * U0moon / (r0mrgl*r0mrgl);

        x = e_sun[0]; y = e_sun[1]; z = e_sun[2]; 
        legendr(2,z,&q); legendr(3,z,&r);
        v = pow( 1/rrsun, 3.) * q
          + pow( 1/rrsun, 4.) * r * (rglov/r0sun);
        v = v * U0sun * pow(rglov/r0sun,2);

     /* the tide potential */
        u =  u+v;      /* U_T in m**2/s**2 */

#if(0)
     /*  the elastisity constant */
        u = u * 5.e-6;   /* mu in s^2/m */
#endif

        printf("%g %13e\n", 
                (abmin+abmin0)/(24.*60),u);

        fprintf(outstr,"%g %13e\n", 
                (abmin+abmin0)/(24.*60),u);


     /* renew loop variables */

       strcpy(dstrt,date);
       x = tstep / 24.;
       days_add(date,dstrt,x);
 
       days_from(dend,date,&xdays); 

       abmin = abmin + tstep*60;



  } /*while(xdays => 0)*/

}
/*===========================================================================*/
/*     Transform a direction unit vector given by RA and DEC in the Gframe
   into SP8 local coodinate at a given DT on a day counted from start of 97. 
*/
int G2SP8(double ra_rad, double dec_rad, double dt, double days97, double e[]){

double longi=134.5, latit=34.9;
/*     ^^^^^^^^^^^^^^^^^^^^^^^   */
double julian,time;
double gst, cphi, sphi, cmu, smu;
double x,y,z,v;

double rads;
rads = PAI/180.;

/* evaluate gst */
      julian = -1096.5 + days97 + dt/24.;
      time = julian / 36525.;
      gst = 100.46 + 36000.77 * time + 15*dt;
      gst = fmod(gst,360);
      if(gst<0) gst=gst+360;  /* in degree */

/* sin and cos of angles */
      sphi = sin( rads*( gst + longi ) );
      cphi = cos( rads*( gst + longi ) );
      smu = sin( rads*( 90. - latit) );
      cmu = cos( rads*( 90. - latit) );

/* the unit vector given by RA and DEC in the Gframe*/
      x = cos(dec_rad)*cos(ra_rad);
      y = cos(dec_rad)*sin(ra_rad);
      z = sin(dec_rad);

/* the unit vector e in the local frame */
      e[0] = cmu*cphi*x + cmu*sphi*y - smu*z;
      e[1] = -sphi*x + cphi*y ;
      e[2] = smu*cphi*x +  smu*sphi*y + cmu*z;
     
      v = e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
      v = sqrt( v );

       e[0] = e[0] / v;
       e[1] = e[1] / v;
       e[2] = e[2] / v;

#ifdef DEBUG
       if(v > 1 )
          printf("W...norm( e ) = %e > 1\n", v);
#endif
/*
        {printf("E...norm e =%e\n"
         "e = %e, %e, %e\n"
         "x,y,z = %e, %e, %e\n"
         "cphi=%e, sphi=%e, cmu=%e smu=%e\n"
         ,v,e[0],e[1],e[2],x,y,z,cphi,sphi,cmu,smu  );}
*/
     
     return(0);
}
/*===========================================================================*/
#if(0)
int Grwich(char *date, double UT,


double glong, double glat, double days, double ut,
             double *juli, double *stime_rad, 
             double *g0dist, double *gra_rad, double *gdec_rad,
             double *tra_rad, double *tdec_rad, double *tdist){
double julian, stime, time;
double xc,yc,zc,rc,x,y,z;
double long_moon,lat_moon,para_moon,rm;
double gra,gdec,tra,tdec,tp;
double rads;
rads = PAI/180;
#endif
/*===========================================================================*/
int moon_pos(double glong, double glat, double days, double ut,
             double *juli, double *stime_rad, 
             double *g0dist, double *gra_rad, double *gdec_rad,
             double *tra_rad, double *tdec_rad, double *tdist){
double julian, stime, time;
double xc,yc,zc,rc,x,y,z;
double long_moon,lat_moon,para_moon,rm;
double gra,gdec,tra,tdec,tp;
double rads;
rads = PAI/180;


/*printf("glong=%f    ",glong);*/


julian = -1096.5 + days + ut/24;
time = julian / 36525;
/*printf("time=%f",time);*/

stime = 100.46 + 36000.77 * time + glong + 15*ut;
/*printf("    stime=%f\n",stime);*/

stime = fmod(stime,360);
if(stime<0) stime=stime+360;

/* geocentric longitude */
long_moon = rads*218.32 + rads*481267.883*time;
long_moon = long_moon + 6.29*rads * sin(rads*fmod(134.9+477198.85*time,360) );
long_moon = long_moon - 1.27*rads * sin(rads*fmod(259.2-413335.38*time,360) );
long_moon = long_moon + 0.66*rads * sin(rads*fmod(235.7+890534.23*time,360) );
long_moon = long_moon + 0.21*rads * sin(rads*fmod(269.9+954397.7*time,360) );
long_moon = long_moon - 0.19*rads * sin(rads*fmod(357.5+35999.05*time,360) );
long_moon = long_moon - 0.11*rads * sin(rads*fmod(186.6+966404.05*time,360) );
long_moon = fmod(long_moon,2*PAI);
if(long_moon<0) long_moon = long_moon + 2.*PAI;

/* geocentric latitude */
lat_moon = 5.13*rads * sin( rads*fmod(93.3+483202.03*time,360) );
lat_moon = lat_moon + 0.28*rads * sin( rads*fmod(228.2+960400.87*time,360) );
lat_moon = lat_moon - 0.28*rads * sin( rads*fmod(318.3+6003.18*time,360) );
lat_moon = lat_moon - 0.17*rads * sin( rads*fmod(217.6-407332.2*time,360) );

/* geocentric RA and DEC */
x = cos(lat_moon)*cos(long_moon);
y = .9175 * cos(lat_moon)*sin(long_moon) - .3978 * sin(lat_moon);
z = .3978 * cos(lat_moon)*sin(long_moon) + .9175 * sin(lat_moon);
gra=atan(y/x); if(x<0) gra = gra+PAI; if(y<0||x>0) gra=gra+2*PAI;
gra=fmod(gra,2*PAI); if(gra<0) gra=gra+2.*PAI;
gdec=asin(z);

/* geocentric horizontal parallax */
para_moon = 0.9508*rads;
para_moon = para_moon + 0.0518*rads*cos(rads*fmod(134.9+477198.85*time,360) );
para_moon = para_moon + 0.0095*rads*cos(rads*fmod(259.2-413335.38*time,360) );
para_moon = para_moon + 0.0078*rads*cos(rads*fmod(235.7+890534.23*time,360) );
para_moon = para_moon + 0.0028*rads*cos(rads*fmod(269.9+954397.7*time,360) );
rm=1/sin(para_moon);

/* topocentric coodinates */
xc = rm*x -cos(rads*glat)*cos(rads*stime);
yc = rm*y -cos(rads*glat)*sin(rads*stime);
zc = rm*z -sin(rads*glat);
rc=sqrt(xc*xc+yc*yc+zc*zc);
tra=atan(yc/xc);  if(xc<0) tra = tra+PAI; if(yc<0||xc>0) tra=tra+2*PAI;
tra=fmod(tra,2*PAI); if(tra<0) tra=tra+2*PAI;
tdec=asin(zc/rc);
tp=asin(1/rc);

#if(0)
/* print result */
printf("              Julian: %f\n",julian+2451545);
printf("       siderial time: %f\n",stime/15);
printf("     geocentric long: %f\n",long_moon/rads);
printf("      geocentric lat: %f\n",lat_moon/rads);
printf(" geocentric parallax: %f\n",para_moon/rads);
printf("            distance: %f\n",rm);
printf("       geocentric RA: %f\n",gra/rads/15);
printf("      geocentric DEC: %f\n",gdec/rads);
printf("      topocentric RA: %f\n",tra/rads/15);
printf("     topocentric DEC: %f\n",tdec/rads);
printf("topocentric parallax: %f\n",tp/rads);
printf("topocentric distance: %f\n",rc);
#endif

         *juli = julian;
         *stime_rad = stime*rads;
         *g0dist = rm;
         *gra_rad = gra;
         *gdec_rad = gdec;
         *tra_rad = tra;
         *tdec_rad = tdec;
         *tdist = rc;

         return(0);
}
/*===========================================================================*/
/*   Give DT in hours to the variable dt 
     Give days since start of 1997       */
int sun_pos(double days, double dt, double *ra_rad, double *dec_rad){
/*double x,y,z,x0,y0,z0,phi0,dphi;*/
double x,y,z,x0,y0,z0,dphi;
double tspan;

double day0_97, gst0, ra0, dec0;
double cjku_dmin[3]={23., 26., 21.448}, cjku, period;
/*cjku_dmin[0]=25.; cjku_dmin[1]=0.; cjku_dmin[2]=0.;*/
/*period = 365.25636*24.;*/
period = (365.*24. +5.) + 49./60.;
dmin2deg(&y,cjku_dmin[0],cjku_dmin[1],cjku_dmin[2]);
cjku = ( y/180. )*PAI;

      day0_97 = 150; /* MAY 30, 1997 0:0:0DT*/
      hmin2deg( &gst0, 16.,30.,10.9);
      hmin2deg( &ra0, 4.,27.,36.46);
      dmin2deg( &dec0, 21.,44.,23.7);

#if(0)
      day0_97 = 2; /* JAN 2, 1997 0:0:0DT*/
      hmin2deg( &gst0,  6., 46., 40.7);
      hmin2deg( &ra0, 18., 50., 33.44);
      dmin2deg( &dec0, -22., -55., -42.1);
#endif

      ra0 = ( ra0/180. )*PAI;
      dec0 = ( dec0/180. )*PAI;

/* n0 in G frame */
      x0 = cos( dec0 )*cos( ra0 );
      y0 = cos( dec0 )*sin( ra0 );
      z0 = sin( dec0 );
      /*if(days == 150) printf("D...Gframe x=%f, y=%f, z=%f\n",x0,y0,z0);*/

/* G -> S0 frame */
      x = x0;
      y =  cos(cjku)*y0 + sin(cjku)*z0;
      z = -sin(cjku)*y0 + cos(cjku)*z0;
      x0=x;y0=y;z0=z;
/*
      phi0 = atan(y0/x0); 
      if(x0<0) phi0=phi0+PAI; if(y0<0||x0>0) phi0=phi0+2*PAI;
      phi0=fmod(phi0,2*PAI); if(phi0<0) phi0=phi0+2*PAI;
*/
/* time span in hours from time0 */
      tspan = days*24. + dt - day0_97*24.;
/* delta phi at the time measured from time0 */
      dphi = (tspan/period)*2.*PAI;
      /*if(days == 150) printf("D...dphi=%f\n",dphi);*/

/* direction vector in the S0 frame */
      x = x0 * cos(dphi) - y0 * sin(dphi);
      y = x0 * sin(dphi) + y0 * cos(dphi);
      z = z0;
      x0=x;y0=y;z0=z;

/* -> G frame */
      x = x0;
      y =  cos(cjku)*y0 - sin(cjku)*z0;
      z =  sin(cjku)*y0 + cos(cjku)*z0;
      /*if(days == 150) printf("D...Gframe x=%f, y=%f, z=%f\n",x,y,z);*/

/* evaluate RA and DEC of the sun at the time */
      *ra_rad = atan(y/x); 
      if(x<0) *ra_rad=*ra_rad+PAI; 
      if(y<0||x>0) *ra_rad=*ra_rad+2.*PAI;
      *ra_rad=fmod(*ra_rad,2.*PAI); if(*ra_rad<0) *ra_rad=*ra_rad+2.*PAI;

      *dec_rad = asin(z);

      /*if(days == 150) printf("D...ra0=%f, dec0=%f, ra=%f, dec=%f\n"
               ,ra0,dec0,*ra_rad,*dec_rad);*/

      return(0);
}
/*===========================================================================*/
int legendr(int n, double z, double *P){

   *P=0.;
   if( n == 0 ){ *P = 1.;  return(0); }
   if( fabs(z) > 1. ){ printf("E...legendr: z = %g\n",z); return(1); }
   if( n == 1 ){ *P = z;  return(0); }
   else if( n == 2 ){ *P = (3.*z*z -1)/2.;  return(0); }
   else if( n == 3 ){ *P = (5.*z*z*z -3.*z)/2.;  return(0); }
   else{ printf("E...legendr: n = %d\n",n); return(1); }

}
