/*
 * Assignment 
 * 
 * Adapt these routines to handle distributed arrays 
 * 
 * Created by G.P. Brandino, I. Girotto, R. Gebauer
 * Last revision: March 2016
 */

#include <stdio.h>
#include "utilities.h"

int FileExists(const char *filename)
{    
   FILE *fp = fopen (filename, "r");
   if (fp!=NULL) fclose (fp);
   return (fp!=NULL);
}


void plot_data_1d( char* name, int n1, int n2, int n3, int dir, double* data)
{
    int i1, i2, i3;
    FILE *fp;
    int num = 1;    
    char buf[256];
    int index;

    snprintf(buf, sizeof(buf), "%s_%d.dat", name, num); 
    while (FileExists(buf))
          {
          num++;
          snprintf(buf, sizeof(buf), "%s_%d.dat", name, num);
          }
    fp = fopen (buf, "w");

 /*
  * HINT: Assuming you sliced your system along i3, iif idir.eq.1 or idir.eq.2 you 
  *       need to understand which process holds the line you are inrested in. If idir.eq.3
  *       you need to take the correct slice of the line from each process
  */

    if ( dir == 1)
        {
        i2=n2/2-1;
        i3=n3/2-1;
        for (i1 = 0; i1 < n1; ++i1)
            {
            index = index_f(i1,i2,i3,n1,n2,n3);  
            fprintf(fp, " %14.6f ", data[index] );
            }
        fprintf(fp, "\n");
        }
    else if ( dir == 2)
        {
        i1=n1/2-1;
        i3=n3/2-1;
        for (i2 = 0; i2 < n2; ++i2)
            {
            index = index_f(i1,i2,i3,n1,n2,n3);    
            fprintf(fp, " %14.6f ", data[index] );
            }
        fprintf(fp, "\n");
        }
    else if ( dir == 3)
        {
        i1=n1/2-1;
        i2=n2/2-1;
        for (i3 = 0; i3 < n3; ++i3)
            {
            index = index_f(i1,i2,i3,n1,n2,n3);    
            fprintf(fp, " %14.6f ", data[index] );
            }
        fprintf(fp, "\n");
        }
    else
        fprintf(stderr, " Wrong value for argument 5 in plot_data_2d \n");

    fclose(fp);
}

void plot_data_2d( char* name, int n1, int n2, int n3, int dir, double* data)
{
    int i1, i2, i3;
    FILE *fp;
    int num = 1;    
    char buf[256];
    int index;

    snprintf(buf, sizeof(buf), "%s_%d.dat", name, num); 
    while (FileExists(buf))
          {
          num++;
          snprintf(buf, sizeof(buf), "%s_%d.dat", name, num);
          }
    fp = fopen (buf, "w");
  /*
   * HINT: Assuming you sliced your system along i3, iif idir==1 or idir==2 you 
   *       need to take the correct slice of the plane from each process. If idir==3, you need  
   *       to understand which process holds the plane you are inrested in. 
   */
    if ( dir == 1)
        {
        i1=n1/2-1;
        for (i2 = 0; i2 < n2; ++i2)
            {
            for (i3 = 0; i3 < n3; ++i3)
                {
                index = index_f(i1,i2,i3,n1,n2,n3);  
                fprintf(fp, " %14.6f ", data[index] );
                }
            fprintf(fp, "\n");
            }
        }
    else if ( dir == 2)
        {
        i2=n2/2-1;
        for (i1 = 0; i1 < n1; ++i1)
            {
            for (i3 = 0; i3 < n3; ++i3)
                {
                index = index_f(i1,i2,i3,n1,n2,n3);
                fprintf(fp, " %14.6f ", data[index] );
                }
            fprintf(fp, "\n");
            }
        }
    else if ( dir == 3)
        {
        i3=n3/2-1;
        for (i1 = 0; i1 < n1; ++i1)
            {
            for (i2 = 0; i2 < n2; ++i2)
                {
                index = index_f(i1,i2,i3,n1,n2,n3);
                fprintf(fp, " %14.6f ", data[index] );
                }
            fprintf(fp, "\n");
            }
        }
    else
        fprintf(stderr, " Wrong value for argument 5 in plot_data_2d \n");

    fclose(fp);
}
