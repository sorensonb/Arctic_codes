/*----------------------------------------------------------------------------/
 *
 *
 * SYNTAX:
 *  $ gcc -o misr_exec netcdf_helper_newest.c -I/usr/include/ -lnetcdf -lm
 *
 * MODIFICATIONS:
 *  Shawn L. Jaker <shawn.jaker@und.edu>        - 2019/02/26
 *      Written 
 *
 * --------------------------------------------------------------------------*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <dirent.h>
//#include "/usr/include/netcdf.h"
#include "./inc/hdf.h"
#include "./inc/netcdf.h"
#include <math.h>

// Grid numbers for a 2.5x5.0 grid
#define NLON  72
#define NLAT  72
// Grid numbers for a 1x1 grid
//#define NLON  360
//#define NLAT  180

#define  FIELD_SIZE     1024         /* maximum length of all the field names */
#define  MAX_DIMS         20         /* maximum number of dimensions in sds  */
#define  MAX_SDS         200         /* maximum number of sds  */
#define  MAX_VDATA       200         /* maximum number of vdata  */
#define  IN_TAIL0      "hdf"          /* possible input filename extension */
#define  IN_TAIL1      "HDF"          /* possible input filename extension */
#define  MAX_CHAR       800          /* maximum characters in a line */
#define  MAX_WORD       15           /* maximum characters in a word */

#define  PI             3.1415926535897932384626433

#define DEBUG 0

//#define MAX_DIMS 3

////////////////////////////////////////////////////////////////////////////
///* read HDF data size  assume 2-D*/
///* input hdf filename, variable string name  */
///* output data size  */
//////////////////////////////////////////////////////////////////////////////
void  readhdfdim(char *filen, char *varname, int *lsize, int *readstatus)
{


   //calibration paramters
   float64  gain, gain_err, offset, offset_err;

   //paramters for reading sds
   int32   sd_id, sds_id, file_id, sds_idx;
   int32   rank, attributes, num_type;
   int32   cal_data_type, num_element ;
   int32   status, i,j;
   int32   dim_sizes[MAX_VAR_DIMS];
   int32   start[MAX_DIMS], edges[MAX_DIMS];
   char    sds_name[64];

//printf("%s\n", varname);

/*
 * Open connections to original (input)  HDF file
*/
     file_id = Hopen (filen, DFACC_READ, 0);
     if ( file_id == -1 ) {
        file_id=0;
        printf("readhdfdm: Error opening input file : %s\n", filen);
        printf("Returning to main\n");
        *readstatus = -1;
        return;
        //exit(-1);
     }

    sd_id = SDstart(filen, DFACC_READ);
    if (sd_id == -1) {
        printf ("SDstart failed.\n");
        HEprint (stdout, 0);
        exit (-1);
    }

   num_element = 1;
   sds_idx = SDnametoindex (sd_id, varname);
   if(sds_idx == -1) {
       printf ("SDnametoindex failed.\n");
       exit (-1);
   }

   sds_id = SDselect (sd_id, sds_idx);
   status = SDgetinfo(sds_id, sds_name, &rank, dim_sizes, &num_type, &attributes);

   for (j = 0; j < rank; j++) {
        num_element *= dim_sizes[j];
   }

   //assume two dimensions
   *lsize = num_element;

   status = SDendaccess(sds_id);

   status = SDend(sd_id);

}


////////////////////////////////////////////////////////////////////////////
///* read HDF data
///* input hdf filename, variable string name  */
///* output data for the given input variable string name, and size  */
//////////////////////////////////////////////////////////////////////////////
void  readhdf(char *filen, char *varname, int *llsize,
              float64  *gain, float64 *offset, float *tmp, int *readstatus)
{


   //calibration paramters
   float64   gain_err, offset_err;

   //data in different format
   char8   *c8;
   int8    *i8;
   uint8   *ui8;
   int16   *i16;
   uint16  *ui16;
   int32   *i32;
   uint32  *ui32;
   float   *f32;
   float64 *f64;

   //paramters for reading sds
   int32   sd_id, sds_id, file_id, sds_idx;
   int32   rank, attributes, num_type;
   int32   cal_data_type, num_element ;
   int32   status, i,j;
   int32   dim_sizes[MAX_VAR_DIMS];
   int32   start[MAX_DIMS], edges[MAX_DIMS];
   char    sds_name[64];


     file_id = Hopen (filen, DFACC_READ, 0);
     if ( file_id == -1 ) {
        file_id = 0;
        printf("readhdf: Error opening input file : %s\n", filen);
        printf("Returning to main\n");
        *readstatus = -1;
        return;
        //exit(-1);
     }

    sd_id = SDstart(filen, DFACC_READ);
    if (sd_id == -1) {
        printf ("SDstart failed.\n");
        HEprint (stdout, 0);
        exit (-1);
    }

   num_element = 1;
   sds_idx = SDnametoindex (sd_id, varname);
   if(sds_idx == -1) {
       printf ("SDnametoindex failed.\n");
       exit (-1);
   }
   sds_id = SDselect (sd_id, sds_idx);
   status = SDgetinfo(sds_id, sds_name, &rank, dim_sizes, &num_type, &attributes);
   status = SDgetcal(sds_id, gain, &gain_err, offset, &offset_err, &cal_data_type );


   for (j = 0; j < rank; j++) {
      num_element *= dim_sizes[j];
      edges[j]  = dim_sizes[j];
      start[j]  = 0;
   }

   //assume one dimensions
   *llsize = dim_sizes[0];
//   *llsize1 = dim_sizes[1];

   switch (num_type)
   {
      case DFNT_FLOAT32: f32 = (float *) malloc( num_element * sizeof(float)); break;
      case DFNT_FLOAT64: f64 = (float64 *) malloc( num_element * sizeof(float64)); break;
      case DFNT_INT8: ui8 = (uint8 *) malloc( num_element * sizeof(uint8)); break;
      case DFNT_UINT8: ui8 = (uint8 *) malloc( num_element * sizeof(uint8)); break;
      case DFNT_INT16: i16 = (int16 *) malloc( num_element * sizeof(int16)); break;
      case DFNT_UINT16: ui16 = (uint16 *) malloc( num_element * sizeof(uint16)); break;
      case DFNT_INT32: i32 = (int32 *) malloc( num_element * sizeof(int32)); break;
      case DFNT_UINT32: ui32 = (uint32 *) malloc( num_element * sizeof(uint32)); break;
      case DFNT_CHAR8: c8 = (char8 *) malloc( num_element * sizeof(char8)); break;
      default: printf("not valid data type\n"); break;
   }


   switch (num_type)
   {
      case DFNT_FLOAT32: status = SDreaddata (sds_id, start, NULL, edges, f32); break;
      case DFNT_FLOAT64: status = SDreaddata (sds_id, start, NULL, edges, f64); break;
      case DFNT_INT8:   status = SDreaddata (sds_id, start, NULL, edges, ui8); break;
      case DFNT_UINT8:   status = SDreaddata (sds_id, start, NULL, edges, ui8); break;
      case DFNT_INT16:   status = SDreaddata (sds_id, start, NULL, edges, i16); break;
      case DFNT_UINT16:  status = SDreaddata (sds_id, start, NULL, edges, ui16); break;
      case DFNT_INT32:   status = SDreaddata (sds_id, start, NULL, edges, i32); break;
      case DFNT_UINT32:  status = SDreaddata (sds_id, start, NULL, edges, ui32); break;
      case DFNT_CHAR8:   status = SDreaddata (sds_id, start, NULL, edges, c8); break;
      default: printf("not valid data type\n"); break;
   }

   if ( status != 0 )
      printf("\n SDreaddata failed on data set %s. \n", sds_name);

   for(i=0; i< num_element; i++) {

     switch (num_type)
       {
           case DFNT_FLOAT32: tmp[i]= f32[i]; break;
           case DFNT_FLOAT64: tmp[i]= f64[i]; break;
           case DFNT_INT8:    tmp[i]= ui8[i]; break;
           case DFNT_UINT8:   tmp[i]= ui8[i]; break;
           case DFNT_INT16:   tmp[i]= i16[i]; break;
           case DFNT_UINT16:  tmp[i]= ui16[i]; break;
           case DFNT_INT32:   tmp[i]= i32[i]; break;
           case DFNT_UINT32:  tmp[i]= ui32[i]; break;
           case DFNT_CHAR8:   tmp[i]= c8[i]; break;
           default: printf("not valid data type\n"); break;
      }

   }  //endfor

   status = SDendaccess(sds_id);

   status = SDend(sd_id);

}


int main(int argc, char *argv[]) {
   
    // Check command line arguments
    if(argc != 3)
    {
        //printf("SYNTAX: ./misr_exec data_file_name average_file\n");
        printf("SYNTAX: ./autoceres_exec year month\n");
        return 0;
    }

    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
    
    char      data1[] = "/home/shared/CERES/SSF/SSF_Terra-FM1_Ed3A/";
    char      data2[] = "/home/shared/CERES/SSF/SSF_Terra-FM2_Ed3A/";
    char      data3[] = "/home/shared/CERES/SSF/SSF_Aqua-FM3_Ed3A/";
    char      data4[] = "/home/shared/CERES/SSF/SSF_Aqua-FM4_Ed3A/";
    char      tfile[] = "CER_SSF_Aqua-FM4-MODIS_Edition3A_300301.2005032923";
    char      temp_path[] = "/home/shared/CERES/SSF/SSF_Terra-FM1_Ed3A/CER_SSF_Terra-FM1-MODIS_Edition3A_304305.2015071721";
    char      *outfile_name,*met_check;
    char      *check_date,*input_date,*make_data1,*make_data2;
    char*     fm2_dates[9];

    float     avg_swf,tlat,avg_lwf;
    float     *lat,*lon,*swf,*lwf,*cld_land,*cld_ocean;
    float64   gain, gain_err, ofst, offset_err;
    //float     **ceres_lwf;     // average AOD array
    //float     **ceres_swf;     // average AOD array
    float     ceres_swf[NLAT][NLON];     // average AOD array
    float     ceres_lwf[NLAT][NLON];     // average AOD array
    float     latgrid[NLAT];     // Latitude array
    float     longrid[NLON];     // Latitude array
    float     minlat,maxlat,minlon,maxlon,tempvar;

    int       i,j,nelement,llsize0,llsize1,llsize2,fm2_check;
    int       status1,status2; // read status
    int       pathlen,slen,lat_index,lon_index,fm2_index;
    int       swf_count,lwf_count;
    //int       **ceres_lwf_cc;      // count array
    //int       **ceres_swf_cc;      // count array
    int       ceres_swf_cc[NLAT][NLON];      // count array
    int       ceres_lwf_cc[NLAT][NLON];      // count array
    int32     llsize[2];    /* lat lon size */

    FILE *fout;

    struct dirent *de;
    DIR *dr;

    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

    status1 = 0;
    status2 = 0;

    ////float misr_aod_unc[NLAT][NLON];     // average AOD array

    //check_month = strndup(tfile+40,6);
    //printf("%s\n",check_month);
    //free(check_month);

    ////// Generate the output file name   
    ////outfile_name = (char*)malloc((strlen(avg_file_start)+10)*sizeof(char)+1);
    ////strcpy(outfile_name,avg_file_start);
    ////strcat(outfile_name,argv[1]);
    ////strcat(outfile_name,argv[2]);
    ////strcat(outfile_name,".txt");
    ////input_date = (char*)malloc((6*sizeof(char))+1);
    ////strcpy(input_date,argv[1]);
    ////strcat(input_date,argv[2]);
    
    // Fill up array to hold the dates where data from FM2 should be used
    fm2_dates[0] = "200102";
    fm2_dates[1] = "200103";
    fm2_dates[2] = "200104";
    fm2_dates[3] = "200108";
    fm2_dates[4] = "200109";
    fm2_dates[5] = "200110";
    fm2_dates[6] = "200601";
    fm2_dates[7] = "200602";
    fm2_dates[8] = "999999";
    fm2_index = 0;

    minlat=-90.;
    maxlat=-90.;
    minlon=0.;
    maxlon=360.;

    // FIGURE OUT HOW TO DEAL WITH INDEXING FOR ONLY LOOKING AT NORTH OF 30 N
    // Allocate arrays
    // Initialize arrays 
    for(i=0;i<NLAT;i++)
    {
        // Fill in the lat grid array
        latgrid[i] = -88.75+2*i;
        for(j=0;j<NLON;j++)
        {
            ceres_swf[i][j]=-999.;
            ceres_swf_cc[i][j]=-9;
            ceres_lwf[i][j]=-999.;
            ceres_lwf_cc[i][j]=-9;
        }
    }

    // Fill up the longrid array
    for(j=0;j<NLON;j++)
    {
        longrid[i] = 2.5+5*i;
    }


    if(fm2_check==-1)
    {
        pathlen = strlen(data1);
        dr = opendir(data1);
    }
    else
    {
        pathlen = strlen(data2);
        dr = opendir(data2);
    }

    if(dr==NULL)
    {
        printf("Could not open directory\n");
        return 0;   
    }

    // Open output file
    fout = fopen("ceres_tester_output.txt","w");
    if(fout == NULL)
    {
        printf("ERROR: Error opening output file. Exiting\n");
        return 0;
    }

    // Read values from the file
    readhdfdim(temp_path,"Colatitude of CERES FOV at surface", &nelement, &status1);
    // If the reading here failed, assume that there is a
    // problem with the entire file. Don't attempt to do
    // anything more with this file.
    if(status1 == -1)
    {
        printf("Problem reading Colatitude data for file %s\n",temp_path);
        status1=0;
    }
    else
    {
        lat = (float *) malloc( nelement * sizeof(float));
        readhdf(temp_path,"Colatitude of CERES FOV at surface", &llsize0,  &gain, &ofst, lat,&status1);

        readhdfdim(temp_path,"Longitude of CERES FOV at surface", &nelement,&status1);
        lon = (float *) malloc( nelement * sizeof(float));
        readhdf(temp_path,"Longitude of CERES FOV at surface", &llsize1,  &gain, &ofst, lon,&status1);

        //llsize[0] = llsize0;
        //printf("Lat size = %d Lon size = %d\n",llsize0,llsize1);

        readhdfdim(temp_path,"CERES SW TOA flux - upwards", &nelement,&status1);
        swf = (float *) malloc( nelement * sizeof(float));
        readhdf(temp_path,"CERES SW TOA flux - upwards", &llsize2,  &gain, &ofst, swf,&status1);

        readhdfdim(temp_path,"CERES LW TOA flux - upwards", &nelement,&status1);
        lwf = (float *) malloc( nelement * sizeof(float));
        readhdf(temp_path,"CERES LW TOA flux - upwards", &llsize0,  &gain, &ofst, lwf,&status1);

        readhdfdim(temp_path,"PSF-wtd MOD04 cloud fraction land", &nelement,&status1);
        cld_land = (float *) malloc( nelement * sizeof(float));
        readhdf(temp_path,"PSF-wtd MOD04 cloud fraction land", &llsize0,  &gain, &ofst, cld_land,&status1);

        readhdfdim(temp_path,"Clear area percent coverage at subpixel resolution", &nelement,&status1);
        cld_ocean = (float *) malloc( nelement * sizeof(float));
        readhdf(temp_path,"Clear area percent coverage at subpixel resolution", &llsize0,  &gain, &ofst, cld_ocean,&status1);

        //readhdfdim(temp_path,"PSF-wtd MOD04 cloud fraction ocean", &nelement,&status1);
        //cld_ocean = (float *) malloc( nelement * sizeof(float));
        //readhdf(temp_path,"PSF-wtd MOD04 cloud fraction ocean", &llsize0,  &gain, &ofst, cld_ocean,&status1);

        for (j=0; j<llsize0; j++) 
        {
            tlat=90.-lat[j]; 
            if((swf[j]<10000.) && (lwf[j]<10000.))
            {
                fprintf(fout,"%f %f %f %f %f %f\n",tlat,lon[j],swf[j],lwf[j],cld_land[j],cld_ocean[j]);
            }
        }
        free(lat);
        free(lon);
        free(swf);
        free(lwf);
        free(cld_land);
        free(cld_ocean);
    }

    fclose(fout);
  
    return 0;
   
}
