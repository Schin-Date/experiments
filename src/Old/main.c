                                  /*   created on Jun. 05, 1997 by S.D. @SP8 
21Nov94

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
double glong=134.5, glat=34.9;
double year,month,day,hrs,mins,sec,ut,UT,days_97,tstep,abmin,abmin0;
char date[]="yyyy/mm/dd hh:mm:ss";
char dend[]="yyyy/mm/dd hh:mm:ss", date0[]="yyyy/mm/dd hh:mm:ss";
char place[]="E000.0/N00.0";    /* Greenwich on equator*/
int ndays;
int days_in_month[12]={31,28,31,30,31,30,31,31,30,31,30,31};
char monthnam[12][3]=
    {"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"};

FILE *outstr;
char outfile[30];
char *data_sun, *data_moon;
char title[30],buffer[80];
double dat_ra_moon_hmin[3], dat_dec_moon_dmin[3];
double dat_ra_moon, dat_dec_moon;
double dat_gst_hmin[3], dat_ra_sun_hmin[3], dat_dec_sun_dmin[3];
double dat_gst, dat_ra_sun, dat_dec_sun;

double julian, stime_rad, g0dist, gra_rad, gdec_rad;
double tra_rad, tdec_rad, tdist; 

double Sra_rad, Sdec_rad;
double e_sun[3], e_moon[3];

double q=0,r=0,s=0,t=0,u=0,v=0,w=0,x=0,y=0,z=0,a=0,b=0,c=0,d=0;
int i,j,k;
double rglov=6.37814;/*Mm*/
double r0sun=149597.870; /*Mm*/
double rrsun=1.016; /* rsun/r0sun    VALUE for June 20 - July 20*/
double r0mrgl=60.2682; /* r0moon/rglov */
double U0sun=8.87210e8; /* m^2/s^2 */
double U0moon=1.27554e4; /* m^2/s^2 */


/*int ifmon;*/

double rads; 
rads = PAI/180.;

/* set defaults for date and location */
      strcpy(date, "1998/02/18 00:00:00"); /*!!!starting time in DT + UT !!!*/
      strcpy(dend, "1998/02/28 00:00:00"); /*!!!end time in DT + UT !!!*/
      strcpy(date0, "1998/02/23 00:00:00"); /*!!!origin time in DT + UT !!!*/
      tstep=1.; /* time step in hours */
      UT = 0.;
      UT = +9.;
      ndays = 8;
      abmin0 = 2880.;  /* definition of time=0 in min */
      strcpy(outfile,"98_01-cel");

      printf("\n\ncalculation starts from: ");
      printf("%s",date);
      sprintf(buffer,"%s",date);
      inp_edt(buffer);
      if( (strlen(buffer)-strlen(date)) == 1 &&
           buffer[strlen(date)]=='q' ) exit(0);
      strcpy(date,buffer);

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
/*printf("buffer=\"%s\"\n",buffer);*/
      sscanf(buffer,"%lf",&tstep);

      days_from(date,date0,&abmin0); abmin0=abmin0*24.*60.;
      days_from(dend,date,&x); ndays = (int) (x+0.5);

      printf("output file name: ");
      printf("%s",outfile);
      sprintf(buffer,"%s",outfile);
      inp_edt(buffer);
      strcpy(outfile,buffer);

      printf("\n");
      printf("outfile=\"%s\"\n",outfile);
      printf("calculation for %d days from %s to %s\n",ndays,date, dend);
      printf("origin of time = %s (abmin0 = %f min)\n",date0,abmin0);
      printf("time step = %f hours\n",tstep);


      printf("go");
      sprintf(buffer,"go");
      inp_edt(buffer);
      if(strstr(buffer,"q")!= NULL) exit(0);


      outstr=fopen(outfile,"w");

/* resolve date and location into numbers */
      if( sscanf(date,"%lf/%lf/%lf %lf:%lf:%lf",
            &year,&month,&day,&hrs,&mins,&sec) != 6 ) exit(1);
      if( sscanf(place,"E%lf/N%lf",&glong,&glat) != 2) exit(1);
      /*      if( year != 1997 ) exit(1);*/

      data_moon = "dat_nasa_moon";
      data_sun = "dat_nasa_sun";
      sprintf(title,"%-5d",(int)year);
      printf("\n %s\n",title);
      printf("     local longitude: %f",glong);
      printf("   local latitude: %f   UT = %d\n",glat,UT);
      printf("     date =  %s\n",date);


/* evaluate ut <-- !!! give DT for this variable !!!*/
      ut= hrs - UT + mins/60. + sec/3600.;          
      abmin = 0.;
/* !!!this ut can be nagative and the day assignment
      in days_97 to call moon_pos etc. can be wrong !!!!*/


/*      ifmon = 0; */
/*********************/
/* prepare title to read the moon data and count days since start of 1997*/
  while( ndays >0 ){  /* loop of the monthes */    i=0;

      if( month<1 || month >12 ){ printf("E...month=%f\n\n",month); exit(1);  }

      strncpy(buffer,monthnam[(int)month-1],3);
      strcat(buffer,"  ");


      x = 0;
      for( j = 0; j<(int)month-1; j++ ) x = x + days_in_month[j];

      days_97 = day + x;

      days_97 = days_97 + (year - 1997)*365. ;
       

/*****
      days_from(date,"1997/01/01 00:00:00",&x);  j = (int) (x+1.);
      days_97 = (double)j;      
****/
/**************** valid only for one month shift *****/
/***
         if(ifmon == 1){
           days_97=days_97+days_in_month[(int)month-2];
           ifmon = 0; }
***/
      
/*
*/
#ifdef DEBUG
printf("D...month = %d %.3s %f\n",(int)month,monthnam[(int)month-1],month);
printf("D...month=%d, ndays=%d days_97=%f\n",(int)month,ndays,days_97);
/*printf("D...entering to the loop of days \n"); */
#endif

   while((int)day+i<=days_in_month[(int)month-1] && ndays>0){     

    /* make the day in the month */
     sprintf(&buffer[4],"%2d",(int)day+i);

            while(ut<24.){

    /*--- calculate moon position ---*/       
      if( moon_pos(glong, glat, days_97+i, ut
              ,&julian, &stime_rad, &g0dist, &gra_rad, &gdec_rad
              , &tra_rad, &tdec_rad, &tdist) != 0 ){
        printf("E...moon_pos!=0 for days=%f, ut=%f, glong=%f, glat=%f\n\n"
                ,days_97+i, ut, glong, glat);
        exit(1);   }  

    /*--- calculate the sun position ---*/       
      if( sun_pos(days_97+i, ut, &Sra_rad, &Sdec_rad) !=0 ){
        printf("E...sun_pos!=0 for days=%f, ut=%f\n\n"
                ,days_97+i, ut);
        exit(1);   }  

    /*--- convert to the local vectors ---*/
        G2SP8(Sra_rad, Sdec_rad, ut, days_97+i, e_sun );
        G2SP8(gra_rad, gdec_rad, ut, days_97+i, e_moon );

        x = e_moon[0]; y = e_moon[1]; z = e_moon[2]; 
        a = atan(y/x); 
        if(x<0) a=a+PAI; if(y<0&&x>0) a=a+2.*PAI;
        a=fmod(a,2.*PAI); if(a<0) a=a+2.*PAI;
        b = sqrt(1.-z*z);

        legendr(2,z,&q); legendr(3,z,&r);
        u = pow( r0mrgl/g0dist, 3.) * q
          + pow( r0mrgl/g0dist, 4.) * r / r0mrgl;
        u = u * U0moon / (r0mrgl*r0mrgl);
        /*printf("D... u=%g, g0dist=%g",u,g0dist);*/

        x = e_sun[0]; y = e_sun[1]; z = e_sun[2]; 
        c = atan(y/x); 
        if(x<0) c=c+PAI; if(y<0&&x>0) c=c+2.*PAI;
        c=fmod(c,2.*PAI); if(c<0) c=c+2.*PAI;
        d = sqrt(1.-z*z);


        legendr(2,z,&q); legendr(3,z,&r);
        v = pow( 1/rrsun, 3.) * q
          + pow( 1/rrsun, 4.) * r * (rglov/r0sun);
        v = v * U0sun * pow(rglov/r0sun,2);
        /*printf(", v=%g\n",v);*/


/* arbitral adjustment to draw figures */
        b = 0.0005 + 0.00018*(b-0.35)/0.65;
        d = 0.0003 + 0.00018*(d-0.2)/0.8;
        u =  u+v;      /* U_T in m**2/s**2 */

#ifdef DEBUG
        printf("%g %13e, %13e, %13e, %13e %13e\n", 
                abmin+abmin0,a,b,c,d,u);
#endif

        printf("%g %13e\n", 
                (abmin+abmin0)/(24.*60),u);

        fprintf(outstr,"%g %13e\n", 
                (abmin+abmin0)/(24.*60),u);


#ifdef RDATA
    /*---read the sun data---*/
      strcpy(&title[5],"SUN");
      if( read_sun(data_sun,title,buffer, 
/*          dat_gst_hmin, dat_ra_sun_hmin, dat_dec_sun_dmin) != 0 ) ;*/
          dat_gst_hmin, dat_ra_sun_hmin, dat_dec_sun_dmin) != 0 ) continue;
         hmin2deg(&y,dat_gst_hmin[0],dat_gst_hmin[1],dat_gst_hmin[2]);
         dat_gst = ( y/180. )*PAI;
         hmin2deg(&y,dat_ra_sun_hmin[0],dat_ra_sun_hmin[1],dat_ra_sun_hmin[2]);
         dat_ra_sun = ( y/180. )*PAI;
         dmin2deg(&y,dat_dec_sun_dmin[0],
                     dat_dec_sun_dmin[1],dat_dec_sun_dmin[2]);
         dat_dec_sun = ( y/180. )*PAI;

    /*--- read moon data ---*/
      strcpy(&title[5],"MOON");
      if( read_moon(data_moon,title,buffer,dat_ra_moon_hmin, dat_dec_moon_dmin)
          != 0 ) continue;
    /* convert the data to radians */
      hmin2deg(&y,dat_ra_moon_hmin[0],dat_ra_moon_hmin[1],dat_ra_moon_hmin[2]);
      dat_ra_moon = ( y/180. )*PAI;
      dmin2deg(&y,dat_dec_moon_dmin[0],
                  dat_dec_moon_dmin[1],dat_dec_moon_dmin[2]);
      dat_dec_moon = ( y/180. )*PAI;
#endif

  /*----- printout the result -----*/
#ifdef DIFFER
    /*--- the differences between the data and calculation for the sun ---*/
      x = Sra_rad - dat_ra_sun;
      y = Sdec_rad - dat_dec_sun;
      if( fabs(x) > q ) q = fabs(x);
      if( fabs(y) > r ) r = fabs(y);
      printf("%s %13e %13e\n"  ,buffer,x,y);

    /*--- the differences between the data and calculation for moon---*/
      x= gra_rad - dat_ra_moon;
      y= gdec_rad - dat_dec_moon;
      a= tra_rad - dat_ra_moon;
      b= tdec_rad - dat_dec_moon;
      if( fabs(x) > q ) q = fabs(x);
      if( fabs(y) > r ) r = fabs(y);
      if( fabs(a) > s ) s = fabs(a);
      if( fabs(b) > t ) t = fabs(b);
      printf("%s %13e %13e %13e %13e\n"  ,buffer,x,y,a,b);
#endif
      /*printf("%s %d: %13e, %13e, %13e, %13e\n",
           buffer,(int)days_97+i, dat_ra_sun, dat_dec_sun,Sra_rad,Sdec_rad);*/
    /*  printf("%s %f %13e %13e %13e\n"  
             ,buffer, julian+2451545, dat_gst, stime_rad, stime_rad-dat_gst);*/
    /*  printf("%s %13e\n"  ,buffer,g0dist*rglov);*//* moon distance */


       ut = ut + tstep;  abmin = abmin + tstep*60;
                          }/* time step loop */

       ut = 0;
       i = i + 1;
       ndays = ndays - 1;       

   printf("D...ndays=%d, i=%d, month=%d day=%d\n",
             ndays,i,(int)month,(int)day);



   }/* the day loop */

   month = month + 1; day = 1; 

/*****************!!!!!!!!!!!!!!!!!!!!!**********/
   if( ndays == 0 ) exit(0);
   printf("W...the month is going to be changed.\n");
   printf("    PREPARE THE PROGRAM FOR THIS CASE.\n");
/*   printf("    =====> STOP\n");*/
   printf("D...ndays=%d, i=%d, month=%d day=%d\n",
             ndays,i,(int)month,(int)day);
/*   exit(0);*/
/*****************!!!!!!!!!!!!!!!!!!!!!**********/


  }/* loop of the monthes*/


#ifdef DIFFER
      printf("max abs: %13e %13e\n" , q,r);
      printf("max:  %13e %13e %13e %13e\n" , q,r,s,t);
#endif


#if(0)
    /*--- output of moon position calculator moon_pos ---*/       
printf("              Julian: %f\n",julian+2451545);
printf("       siderial time: %f\n",stime_rad/15);
printf("            distance: %f\n",g0dist);
printf("       geocentric RA: %f\n",gra_rad/rads/15);
printf("      geocentric DEC: %f\n",gdec_rad/rads);
printf("      topocentric RA: %f\n",tra_rad/rads/15);
printf("     topocentric DEC: %f\n",tdec_rad/rads);
printf("topocentric distance: %f\n",tdist);
#endif
}
/*===========================================================================*/
/*     Transform a direction unit vector given by RA and DEC in the Gframe
   into SP8 local coodinate at a given DT on a day counted from start of 97. 
*/
int G2SP8(double ra_rad, double dec_rad, double dt, double days97, double e[]){

double longi=134.5, latit=34.9;

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
#if(0)
/*===========================================================================*/
int chk_sun(){
char *data_sun, *data_moon;
char title[10],date[10],buffer[80];
double gst_hmin[3], ra_sun_hmin[3], dec_sun_dmin[3];
double ra_moon_hmin[3], dec_moon_dmin[3];
double cjku_dmin[3]={23, 26, 21.448}, cjku;

double gst, ra_sun, dec_sun, ra_moon, dec_moon;

double deg;
int i,j,k, day;
double x,y,z,dlta;

    dmin2deg(&deg,cjku_dmin[0],cjku_dmin[1],cjku_dmin[2]);
    cjku = ( deg/180 )*PAI;
    

    data_sun = "dat_nasa_sun";
    data_moon = "dat_nasa_moon";

    strcpy(title,"1997 SUN");


    for(i = 0; i<10; i++){
      day = 20 + i;
      sprintf(date,"MAY %2d",day);
      if( read_sun(data_sun,title,date, gst_hmin, ra_sun_hmin, dec_sun_dmin)
          != 0 ) continue;

         hmin2deg(&deg,gst_hmin[0],gst_hmin[1],gst_hmin[2]);
         gst = ( deg/180 )*PAI;
         hmin2deg(&deg,ra_sun_hmin[0],ra_sun_hmin[1],ra_sun_hmin[2]);
         ra_sun = ( deg/180 )*PAI;
         dmin2deg(&deg,dec_sun_dmin[0],dec_sun_dmin[1],dec_sun_dmin[2]);
         dec_sun = ( deg/180 )*PAI;

       dlta = cos(cjku)*sin(dec_sun)
             - sin(cjku)*cos(dec_sun)*sin(ra_sun);
         printf("%s delta = %e\n",date,dlta);

   }/* for( i */

}
/*===========================================================================*/
int handy_position(){

char *data_sun, *data_moon;
char title[10],date[10],buffer[80];
double gst_hmin[3], ra_sun_hmin[3], dec_sun_dmin[3];
double ra_moon_hmin[3], dec_moon_dmin[3];

double gst, ra_sun, dec_sun, ra_moon, dec_moon;

double deg;


    data_sun = "dat_nasa_sun";
    data_moon = "dat_nasa_moon";


 while(1){
    printf("title: ");   gets(title);
    if( title[0]=='\0') exit(0);

 while(1){
    printf("date: ");   gets(date);
    if( date[0]=='\0') break;

   if( strstr(title,"SUN")!=0 ){
     if( read_sun(data_sun,title,date, gst_hmin, ra_sun_hmin, dec_sun_dmin)
         != 0 ) continue;

        printf("%s   ",date);
        printf("GST_hmin = %2.0f:%2.0f:%3.1f   ", 
                                 gst_hmin[0],gst_hmin[1],gst_hmin[2]);
        printf("RA_hmin = %2.0f:%2.0f:%4.2f   ", 
                         ra_sun_hmin[0],ra_sun_hmin[1],ra_sun_hmin[2]);
        printf("Dec_dmin = %2.0f:%2.0f:%3.1f\n", 
                      dec_sun_dmin[0],dec_sun_dmin[1],dec_sun_dmin[2]);
     hmin2deg(&deg,gst_hmin[0],gst_hmin[1],gst_hmin[2]);
     gst = ( deg/180 )*PAI;
     hmin2deg(&deg,ra_sun_hmin[0],ra_sun_hmin[1],ra_sun_hmin[2]);
     ra_sun = ( deg/180 )*PAI;
     dmin2deg(&deg,dec_sun_dmin[0],dec_sun_dmin[1],dec_sun_dmin[2]);
     /*printf("D...deg = %e\n",deg);*/
     dec_sun = ( deg/180 )*PAI;
        printf("...in radians: GST = %e, RA=%e, Dec = %e\n\n"
                    ,gst,ra_sun,dec_sun);
   }/*SUN*/

   else if( strstr(title,"MOON")!=0 ){
     if( read_moon(data_moon,title,date, ra_moon_hmin, dec_moon_dmin)
         != 0 ) continue;

        printf("%s MOON  ",date);
        printf("RA_hmin = %2.0f:%2.0f:%4.2f   ", 
                         ra_moon_hmin[0],ra_moon_hmin[1],ra_moon_hmin[2]);
        printf("Dec_dmin = %2.0f:%2.0f:%3.1f\n", 
                      dec_moon_dmin[0],dec_moon_dmin[1],dec_moon_dmin[2]);
     hmin2deg(&deg,ra_moon_hmin[0],ra_moon_hmin[1],ra_moon_hmin[2]);
     ra_moon = ( deg/180 )*PAI;
     dmin2deg(&deg,dec_moon_dmin[0],dec_moon_dmin[1],dec_moon_dmin[2]);
     /*printf("D...deg = %e\n",deg);*/
     dec_moon = ( deg/180 )*PAI;
        printf("...in radians: RA=%e, Dec = %e\n\n"
                    ,ra_moon,dec_moon);
   }/*MOON*/

   else break;

  }/* while(1)*/

 }/*the first while(1)*/

}
if((int)days_97+i==150) 
printf("ra: %f;%f;%f = %e\n",
dat_ra_sun_hmin[0],dat_ra_sun_hmin[1],dat_ra_sun_hmin[2],dat_ra_sun);
#endif
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
