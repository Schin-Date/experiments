#include	<stdio.h>
#include	<stdlib.h>
#include	<time.h>
#include	<string.h>
/*-----------------------------------------------------------------------------
                                                 June 24, 1998 S.D.@SP8
                 Approximate Geocentric Moon Position

    Purpose:  To evaluate approximate geocentric coordinates of the Moon

     Usage :   In the main routine, prepare

                   time_t  tint;
                   struct tm *tblock;
                   char atime[128], cntl[128];
                   int rtn;

            and specify one of 

                   atime = 'yyyy/mm/dd hh:mm:ss'   when cntl="char",
                   tint = an integer               when cntl="int"
                   or a set  (*tblock).tm_year,... when cntl="tm"

            to call

                   rtn = texprss(atime, &tblock, &tint, cntl);
             
            to obtain values of the time in other expressions.
                In the cases of cntl="char" or "tm", you need to sandwitch
            your main program between a peir,

                   tblock = calloc(1,sizeof(struct tm));
                   free( tblock );

            See an example attached at the end of this file.
-------------------------------------------------------------------------------*/
int agmoon()

-------------------------------------------------------------------------------*/

int texprss(char *atime, struct tm **tblock, time_t *tint, char *cntl)
{
	char buffer[128], buffer2[128], *tok;


        if( strstr(cntl,"in")!=NULL){

                *tblock = localtime(tint);
                sprintf(atime,"19%02d/%02d/%02d %02d:%02d:%02d",
                 (*tblock)->tm_year, (*tblock)->tm_mon + 1, (*tblock)->tm_mday,
                 (*tblock)->tm_hour, (*tblock)->tm_min, (*tblock)->tm_sec);
                return(0);
	}

        if( strstr(cntl,"tm")!=NULL || strstr(cntl,"bl")!=NULL ){
                          /*printf("D...tkb_yr=%d,  mkt=%d\n",
                            (**tblock).tm_mon,mktime(*tblock));*/
		*tint=mktime(*tblock);
                sprintf(atime,"19%02d/%02d/%02d %02d:%02d:%02d",
                 (*tblock)->tm_year, (*tblock)->tm_mon + 1, (*tblock)->tm_mday,
                 (*tblock)->tm_hour, (*tblock)->tm_min, (*tblock)->tm_sec);
                return(0);
	}

        if( strstr(cntl,"ch")!=NULL){

    		strcpy(buffer2,atime);
 	   	tok=strtok(buffer2,"/");
		strcpy(buffer, tok);
		(*tblock)->tm_year=atoi(buffer)-1900;
		tok=strtok(NULL, "/");
		strcpy(buffer, tok);
		(*tblock)->tm_mon=atoi(buffer)-1;
		tok=strtok(NULL, " ");
		strcpy(buffer, tok);
		(*tblock)->tm_mday=atoi(buffer);
		tok=strtok(NULL, ":");
		strcpy(buffer, tok);
		(*tblock)->tm_hour=atoi(buffer);
		tok=strtok(NULL, ":");
		strcpy(buffer, tok);
		(*tblock)->tm_min=atoi(buffer);
		tok=strtok(NULL, ",");
		strcpy(buffer, tok);
		(*tblock)->tm_sec=atoi(buffer);
		*tint=mktime(*tblock);
                return(0);
	      }
        printf("E...texpress: cntl = '%s' not understood.\n");
        return(1);
}

#ifdef TEST
/*=============================================================================
                       An example of the main routine
*/
main()
{
	time_t	tint;
   	struct tm *tblock;
        char atime[128], cntl[128], buff[128];
        int rtn;

        printf("original time expression{char|tm|int} = ");
        gets(cntl);

        if( strstr(cntl,"in")!=NULL){
            printf("\ntint = ");
            gets(buff); sscanf(buff,"%d",&tint);

            rtn = texprss(atime, &tblock, &tint, cntl);
            
            printf("-> atime = '%s'\n",atime);
            printf("-> tblock : 19%02d/%02d/%02d %02d:%02d:%02d\n",
                 tblock->tm_year, tblock->tm_mon + 1, tblock->tm_mday,
                 tblock->tm_hour, tblock->tm_min, tblock->tm_sec);
            /*printf("D...mkt=%d\n",mktime(tblock));*/
	  }
         
        if( strstr(cntl,"tm")!=NULL || strstr(cntl,"bl")!=NULL ){

           tblock = calloc(1,sizeof(struct tm));
            printf("yyyy = ");
            gets(buff); sscanf(buff,"%d",&(tblock->tm_year));
            (*tblock).tm_year = (*tblock).tm_year - 1900;
            printf("mm = ");
            gets(buff); sscanf(buff,"%d",&(tblock->tm_mon));
            tblock->tm_mon = tblock->tm_mon -1;
            printf("dd = ");
            gets(buff); sscanf(buff,"%d",&(tblock->tm_mday));
            printf("hh = ");
            gets(buff); sscanf(buff,"%d",&(tblock->tm_hour));
            printf("mm = ");
            gets(buff); sscanf(buff,"%d",&(tblock->tm_min));
            printf("ss = ");
            gets(buff); sscanf(buff,"%d",&(tblock->tm_sec));
            /*printf("D...mkt=%d\n",mktime(tblock));*/

            rtn = texprss(atime, &tblock, &tint, cntl);
         
            printf("\ntblock : 19%02d/%02d/%02d %02d:%02d:%02d\n",
                 tblock->tm_year, tblock->tm_mon + 1, tblock->tm_mday,
                 tblock->tm_hour, tblock->tm_min, tblock->tm_sec);
            printf("-> atime = '%s'\n",atime);
            printf("-> tint = %d\n",tint);

            free( tblock );
	  }


        if( strstr(cntl,"ch")!=NULL){
           tblock = calloc(1,sizeof(struct tm));

            printf("\natime(yyyy/mm/dd hh:mm:ss) = ");
            gets(atime);

            rtn = texprss(atime, &tblock, &tint, cntl);

        printf("-> tblock : 19%02d/%02d/%02d %02d:%02d:%02d\n->(int)%d\n",
                 tblock->tm_year, tblock->tm_mon + 1, tblock->tm_mday,
                 tblock->tm_hour, tblock->tm_min, tblock->tm_sec,
                 tint);

            free( tblock );

	  }
}
#endif
