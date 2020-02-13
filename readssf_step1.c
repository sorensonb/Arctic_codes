/* programmed by Jianglong Zhang @UND 2017
   for reading SSF data
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "./inc/hdf.h"
#include "./inc/netcdf.h"

#define  FIELD_SIZE     1024         /* maximum length of all the field names */
#define  MAX_DIMS         20         /* maximum number of dimensions in sds  */
#define  MAX_SDS         200         /* maximum number of sds  */
#define  MAX_VDATA       200         /* maximum number of vdata  */
#define  IN_TAIL0      "hdf"          /* possible input filename extension */
#define  IN_TAIL1      "HDF"          /* possible input filename extension */

#define  MAX_CHAR       800          /* maximum characters in a line */
#define  MAX_WORD       15           /* maximum characters in a word */

#define  PI             3.1415926535897932384626433


////////////////////////////////////////////////////////////////////////////
///* read HDF data size  assume 2-D*/
///* input hdf filename, variable string name  */
///* output data size  */
//////////////////////////////////////////////////////////////////////////////
void  readhdfdim(char *filen, char *varname, int *lsize)
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
        printf("Error opening input file : %s\n", filen);
        exit(-1);
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
void  readhdf(char *filen, char *varname, int *llsize0,
              float64  *gain, float64 *offset, float *tmp)
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
        printf("Error opening input file : %s\n", filen);
        exit(-1);
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
   *llsize0 = dim_sizes[0];
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






////////////////////////////////////////////////////////////////////////////
//// Main program
//// Input:
//// Output:
/////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{

   FILE    *meta_fp, *sds_fp, *vd_fp, *an_fp;

   int32   sd_id, sds_id, dim_id, file_id;
   int32   index, attr_index, sds_idx, dim_index;
   int32   nindex, sds_index[MAX_SDS];

   int32   rank, attributes, n_datasets, n_file_attrs, num_type, count;
   int32   n_records, nparm, temparm, cal_data_type, num_element ;
   int32   status, i,j, ifield;

   int32   dim_sizes[MAX_VAR_DIMS];
   int32   start[MAX_DIMS],  stride[MAX_DIMS], edges[MAX_DIMS];

   int32   interlace, vdata_size, field_index, field_size, field_order;
   int32   field_type, vdata_idx, vdata_ref, vdata_id;

   char    fields[FIELD_SIZE], attr_name[64], sds_name[64];

   float64  gain, gain_err, ofst, offset_err;

   char8   *c8;
   int8    *i8;
   uint8   *ui8;
   int16   *i16;
   uint16  *ui16;
   int32   *i32;
   uint32  *ui32;
   float   *f32;
   float64 *f64;

/*--------------------------------------------------------------------*/

   float   *tmpt;

   float   *lper;
   float   *lw;
   float   *sw;
   float   *swr;
   float   *aod;
   float   *cld;
   float   *clr;
   float   *raz;
   float   *vza;
   float   *sza;
   float   *lon;
   float   *lat;
   float   *totr;


   int       nelement,llsize0,llsize1;
   int32     llsize[2];    /* lat lon size */
   int32     datasize[3];    /* data size */
   uint16  *datas;

   FILE      *fp,*tmp;
   char      infile[128],outfile[128],hours[128];

/*--------------------------------------------------------------------*/


    if (argc != 3) {
        printf("Usage: %s <dir>  <CERES SSF_file>  \n\n", argv[0]);
        exit(-1);
     }

//   read ssf data
     strcpy(infile,argv[1]);
     strcat(infile,argv[2]);
     if (Hishdf(infile) == 0) {
       printf("Error: %s is not a valid HDF file, or file not found\n\n", argv[1]);
       exit(-1);
     }

   readhdfdim(infile,"Colatitude of CERES FOV at surface", &nelement);
   tmpt = (float *) malloc( nelement * sizeof(float));
   readhdf(infile,"Colatitude of CERES FOV at surface", &llsize0,  &gain, &ofst, tmpt);
   llsize[0] = llsize0;
   lat = (float *) malloc( llsize[0] * sizeof(float*));
   for (j=0; j<llsize[0]; j++) {
       lat[j]=tmpt[j];
   }
   free(tmpt);

   readhdfdim(infile,"Longitude of CERES FOV at surface", &nelement);
   tmpt = (float *) malloc( nelement * sizeof(float));
   readhdf(infile,"Longitude of CERES FOV at surface", &llsize0,  &gain, &ofst, tmpt);
   llsize[0] = llsize0;
   lon = (float *) malloc( llsize[0] * sizeof(float*));
   for (j=0; j<llsize[0]; j++) {
       lon[j]=tmpt[j];
   }
   free(tmpt);


   readhdfdim(infile,"CERES viewing zenith at surface", &nelement);
   tmpt = (float *) malloc( nelement * sizeof(float));
   readhdf(infile,"CERES viewing zenith at surface", &llsize0,  &gain, &ofst, tmpt);
   llsize[0] = llsize0;
   vza = (float *) malloc( llsize[0] * sizeof(float*));
   for (j=0; j<llsize[0]; j++) {
       vza[j]=tmpt[j];
   }
   free(tmpt);

   readhdfdim(infile,"CERES solar zenith at surface", &nelement);
   tmpt = (float *) malloc( nelement * sizeof(float));
   readhdf(infile,"CERES solar zenith at surface", &llsize0,  &gain, &ofst, tmpt);
   llsize[0] = llsize0;
   sza = (float *) malloc( llsize[0] * sizeof(float*));
   for (j=0; j<llsize[0]; j++) {
       sza[j]=tmpt[j];
   }
   free(tmpt);



   readhdfdim(infile,"CERES relative azimuth at surface", &nelement);
   tmpt = (float *) malloc( nelement * sizeof(float));
   readhdf(infile,"CERES relative azimuth at surface", &llsize0,  &gain, &ofst, tmpt);
   llsize[0] = llsize0;
   raz = (float *) malloc( llsize[0] * sizeof(float*));
   for (j=0; j<llsize[0]; j++) {
       raz[j]=tmpt[j];
   }
   free(tmpt);


   readhdfdim(infile,"Clear area percent coverage at subpixel resolution", &nelement);
   tmpt = (float *) malloc( nelement * sizeof(float));
   readhdf(infile,"Clear area percent coverage at subpixel resolution", &llsize0,  &gain, &ofst, tmpt);
   llsize[0] = llsize0;
   clr = (float *) malloc( llsize[0] * sizeof(float*));
   for (j=0; j<llsize[0]; j++) {
       clr[j]=tmpt[j];
   }
   free(tmpt);


   readhdfdim(infile,"PSF-wtd MOD04 cloud fraction ocean", &nelement);
   tmpt = (float *) malloc( nelement * sizeof(float));
   readhdf(infile,"PSF-wtd MOD04 cloud fraction ocean", &llsize0,  &gain, &ofst, tmpt);
   llsize[0] = llsize0;
   cld = (float *) malloc( llsize[0] * sizeof(float*));
   for (j=0; j<llsize[0]; j++) {
       cld[j]=tmpt[j];
   }
   free(tmpt);

   readhdfdim(infile,"PSF-wtd MOD04 effective optical depth average ocean (0.550)", &nelement);
   tmpt = (float *) malloc( nelement * sizeof(float));
   readhdf(infile,"PSF-wtd MOD04 effective optical depth average ocean (0.550)", &llsize0,  &gain, &ofst, tmpt);
   llsize[0] = llsize0;
   aod = (float *) malloc( llsize[0] * sizeof(float*));
   for (j=0; j<llsize[0]; j++) {
       aod[j]=tmpt[j];
   }
   free(tmpt);


   readhdfdim(infile,"CERES TOT filtered radiance - upwards", &nelement);
   tmpt = (float *) malloc( nelement * sizeof(float));
   readhdf(infile,"CERES TOT filtered radiance - upwards", &llsize0,  &gain, &ofst, tmpt);
   llsize[0] = llsize0;
   totr = (float *) malloc( llsize[0] * sizeof(float*));
   for (j=0; j<llsize[0]; j++) {
       totr[j]=tmpt[j];
   }
   free(tmpt);

   readhdfdim(infile,"CERES SW filtered radiance - upwards", &nelement);
   tmpt = (float *) malloc( nelement * sizeof(float));
   readhdf(infile,"CERES SW filtered radiance - upwards", &llsize0,  &gain, &ofst, tmpt);
   llsize[0] = llsize0;
   swr = (float *) malloc( llsize[0] * sizeof(float*));
   for (j=0; j<llsize[0]; j++) {
       swr[j]=tmpt[j];
   }
   free(tmpt);

   readhdfdim(infile,"CERES SW TOA flux - upwards", &nelement);
   tmpt = (float *) malloc( nelement * sizeof(float));
   readhdf(infile,"CERES SW TOA flux - upwards", &llsize0,  &gain, &ofst, tmpt);
   llsize[0] = llsize0;
   sw = (float *) malloc( llsize[0] * sizeof(float*));
   for (j=0; j<llsize[0]; j++) {
       sw[j]=tmpt[j];
   }
   free(tmpt);

   readhdfdim(infile,"CERES LW TOA flux - upwards", &nelement);
   tmpt = (float *) malloc( nelement * sizeof(float));
   readhdf(infile,"CERES LW TOA flux - upwards", &llsize0,  &gain, &ofst, tmpt);
   llsize[0] = llsize0;
   lw = (float *) malloc( llsize[0] * sizeof(float*));
   for (j=0; j<llsize[0]; j++) {
       lw[j]=tmpt[j];
//printf("%d %f\n", j,lw[j]);
   }
   free(tmpt);

   readhdfdim(infile,"Percentage of CERES FOV with MODIS land aerosol", &nelement);
   tmpt = (float *) malloc( nelement * sizeof(float));
   readhdf(infile,"Percentage of CERES FOV with MODIS land aerosol", &llsize0,  &gain, &ofst, tmpt);
   llsize[0] = llsize0;
   lper = (float *) malloc( llsize[0] * sizeof(float*));
   for (j=0; j<llsize[0]; j++) {
       lper[j]=tmpt[j];
   }
   free(tmpt);

//printf("%d\n", llsize[0]);

    strcpy(hours, argv[2]);

   strncpy(outfile, argv[2], 49);
   strcat(outfile, ".clr");

   tmp=fopen(outfile,"a+");

   for (j=0; j<llsize[0]; j++) {

//      if(aod[j]>=0 && cld[j] <=5 && clr[j] >=95 && lper[j]<0.0000000001 && swr[j] < 50000 && totr[j] < 50000)
      if(clr[j] >=95 && lper[j]<0.0000000001 && swr[j] < 50000 && totr[j] < 50000)
         fprintf(tmp, "%f %f %f %f %f %16.10f %16.10f %f %f %f %f %f %f   %c%c\n", lat[j], lon[j], sza[j],vza[j],raz[j], swr[j],totr[j],sw[j],lw[j],aod[j],cld[j],clr[j], lper[j], hours[49],hours[50]);
   }

   fclose(tmp);

   free(lat);
   free(lper);
   free(lw);
   free(sw);
   free(swr);
   free(aod);
   free(cld);
   free(clr);
   free(raz);
   free(vza);
   free(sza);
   free(lon);
   free(totr);

return(0);
}





