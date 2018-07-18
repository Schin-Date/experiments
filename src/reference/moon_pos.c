                                  /*   created on Jun. 03, 1997 by S.D. @SP8 */
#include <stdio.h>  /* for printf,fgets NULL */
#include <stdlib.h>  /* for exit, malloc */
#include <ctype.h> /* for tolower */
#include <string.h> /* for strcpy */
#include <math.h>   /* to use sqrt */
#include <time.h>
#include "/users/sr/schin/clib/head/const.h" /* for R4OVER and so on */
/*===========================================================================*/
int handy_moon(){
/* void main(){*/
double glong=134.5, glat=34.9;
double year,month,day,days,hrs,mins,sec,ut,UT;
char *date,*place;
double x,y,z;

double julian, stime_rad, g0dist, gra_rad, gdec_rad;
double tra_rad, tdec_rad, tdist; 
double rads;
rads = PAI/180;


date = "1997/05/20 00:00:00";
place = "E000.0/N00.0";      /* Greenwich on equator*/
UT = +9;
UT = 0;


#ifdef DEBUG2
sscanf(date,"%lf/%lf/%lf %lf:%lf:%lf",
      &year,&month,&day,&hrs,&mins,&sec);
printf("year=%f,month=%f,day=%f,hrs=%f,mins=%f,sec=%f\n",
        year,month,day,hrs,mins,sec);
#endif

if( sscanf(date,"%lf/%lf/%lf %lf:%lf:%lf",
      &year,&month,&day,&hrs,&mins,&sec) != 6 ) return(1);
if( sscanf(place,"E%lf/N%lf",&glong,&glat) != 2) return(1);


days = day;
if( year != 1997 ) exit(1);
if(month > 1 ) days=days+31;  if(month > 2 ) days=days+28;
if(month > 3 ) days=days+31;  if(month > 4 ) days=days+30;
if(month > 5 ) days=days+31;  if(month > 6 ) days=days+30;
if(month > 7 ) days=days+31;  if(month > 8 ) days=days+31;
if(month > 9 ) days=days+30;  if(month > 10 ) days=days+31;
if(month > 11 ) days=days+30; if(month > 12 ) return(1);
/*days=31+28+31+30+20;hrs = 0; mins=-1;  */ /*1997 MAY 20 UT=-1min*/


ut= hrs - UT + mins/60;          


printf("     local longitude: %f\n",glong);
printf("      local latitude: %f\n",glat);
printf("                days: %f\n",days);
printf("             hrs(UT): %f\n",hrs);
printf("                mins: %f\n",mins);

    moon_pos(glong, glat, days, ut
              ,&julian, &stime_rad, &g0dist, &gra_rad, &gdec_rad
              , &tra_rad, &tdec_rad, &tdist); 

printf("              Julian: %f\n",julian+2451545);
printf("       siderial time: %f\n",stime_rad/15);
printf("            distance: %f\n",g0dist);
printf("       geocentric RA: %f\n",gra_rad/rads/15);
printf("      geocentric DEC: %f\n",gdec_rad/rads);
printf("      topocentric RA: %f\n",tra_rad/rads/15);
printf("     topocentric DEC: %f\n",tdec_rad/rads);
printf("topocentric distance: %f\n",tdist);

       return(0);

}
/*===========================================================================*/
/*     Give the DT time in the place of ut.  */
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
