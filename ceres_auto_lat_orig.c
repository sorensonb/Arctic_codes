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

#define NLON  360
#define NLAT  180

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
    char      avg_file_start[] = "/home/bsorenson/HighLatitudeStudy/CERES/average_files/TERRA/N20to90/ceres_avg_flux_";
    char      tfile[] = "CER_SSF_Aqua-FM4-MODIS_Edition3A_300301.2005032923";
    char      *temp_path,*outfile_name;
    char      *check_date,*input_date;
    char*     fm2_dates[9];

    float     *tmpt;
    float     avg_swf,tlat,avg_lwf;
    float     *lat,*swf,*lwf;
    float64   gain, gain_err, ofst, offset_err;
    float     ceres_lwf[NLAT];     // average AOD array
    float     ceres_swf[NLAT];     // average AOD array
    float     minlat;

    int       i,j,nelement,llsize0,llsize1,fm2_check;
    int       pathlen,slen,lat_index,lon_index,fm2_index;
    int       swf_count,lwf_count;
    int       ceres_lwf_cc[NLAT];      // count array
    int       ceres_swf_cc[NLAT];      // count array
    int32     llsize[2];    /* lat lon size */

    FILE *fout;

    struct dirent *de;
    DIR *dr;

    // = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

    ////float misr_aod_unc[NLAT][NLON];     // average AOD array

    //check_month = strndup(tfile+40,6);
    //printf("%s\n",check_month);
    //free(check_month);

    // Generate the output file name   
    outfile_name = (char*)malloc((strlen(avg_file_start)+10)*sizeof(char)+1);
    strcpy(outfile_name,avg_file_start);
    strcat(outfile_name,argv[1]);
    strcat(outfile_name,argv[2]);
    strcat(outfile_name,".txt");
    input_date = (char*)malloc((6*sizeof(char))+1);
    strcpy(input_date,argv[1]);
    strcat(input_date,argv[2]);
    
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

    // Determine if the current month is in FM1 or FM2 for Terra
    fm2_check=-1;
    i=0;
    while(i<8)
    {
        if(strcmp(input_date,fm2_dates[i])==0)
        {
            fm2_check=i;
            i=8;
        }
        else
        {
            i++;
        }
    }

    minlat=-20.;

    // FIGURE OUT HOW TO DEAL WITH INDEXING FOR ONLY LOOKING AT NORTH OF 30 N
    // Initialize arrays 
    for(i=0;i<NLAT;i++)
    {
        ceres_lwf[i]=-999.;
        ceres_swf[i]=-999.;
        ceres_lwf_cc[i]=-9;
        ceres_swf_cc[i]=-9;
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
    // Loop over all files in the directory
    while((de=readdir(dr)) != NULL)
    {
        //printf("Current directory %s\n",de->d_name);
        slen = strlen(de->d_name);
        // Ignore '.' and '..'
        if(slen>2)
        {
            // Extract year and month from the directory name
            ///temp_dir = (char*)malloc(slen*sizeof(char)+1);
            ///strcpy(temp_dir,de->d_name);
            check_date = strndup(de->d_name+41,6);
            //printf("%s %s %s\n",tmonth,de->d_name,check_month);
            //dir_year = strtok(temp_dir,".");
            //dir_month = strtok(NULL,".");
            //dir_day = strtok(NULL,"\n");
            //if((strcmp(tmonth,check_month)==0))

            //temp_path = (char*)malloc((pathlen+slen)*sizeof(char)+1);
            //if(strcmp(check_date,fm2_dates[fm2_index])==0)
            //{
            //    strcpy(temp_path,data2);
            //    fm2_index++;
            //}
            //else
            //{
            //}

            if(strcmp(input_date,check_date)==0)
            {
                //printf("    Good file: %s\n",de->d_name);

                // Find the complete path
                temp_path = (char*)malloc((pathlen+slen)*sizeof(char)+1);
                if(fm2_check==-1)
                {
                    strcpy(temp_path,data1);
                }
                else
                {
                    strcpy(temp_path,data2);
                }
                strcat(temp_path,de->d_name);
                //strcpy(temp_path,data1);

                printf("%s\n",temp_path);

                readhdfdim(temp_path,"Colatitude of CERES FOV at surface", &nelement);
                lat = (float *) malloc( nelement * sizeof(float));
                readhdf(temp_path,"Colatitude of CERES FOV at surface", &llsize0,  &gain, &ofst, lat);

                readhdfdim(temp_path,"CERES SW TOA flux - upwards", &nelement);
                swf = (float *) malloc( nelement * sizeof(float));
                readhdf(temp_path,"CERES SW TOA flux - upwards", &llsize0,  &gain, &ofst, swf);

                readhdfdim(temp_path,"CERES LW TOA flux - upwards", &nelement);
                lwf = (float *) malloc( nelement * sizeof(float));
                readhdf(temp_path,"CERES LW TOA flux - upwards", &llsize0,  &gain, &ofst, lwf);

                llsize[0] = llsize0;
                //lat = (float *) malloc( llsize[0] * sizeof(float*));
                //for (j=0; j<llsize[0]; j++) {
                for (j=0; j<llsize[0]; j++) 
                {
                    //lat[j]=tmpt[j];
                    tlat=90.-lat[j]; 
                    if(tlat>minlat) 
                    {
                    //if((tlat>50.) && (tlat<51.)) {
                        lat_index = (int)tlat+90;
                     //   lon_index = (int)lon[index]+180;
                        if(swf[j]<10000.) 
                        {
                            if(ceres_swf_cc[lat_index]==-9) {
                                ceres_swf[lat_index]=swf[j];
                                //avg_swf=swf[j];
                                ceres_swf_cc[lat_index]=1;
                                //swf_count=1;
                            }
                            else {
                //                printf("Old avg = %d %f %f\n",swf_count,tlat,avg_swf);
                                ceres_swf[lat_index] = ((ceres_swf[lat_index]*ceres_swf_cc[lat_index])+
                                    swf[j])/(ceres_swf_cc[lat_index]+1);
                                //avg_swf = ((avg_swf*swf_count)+swf[j])/(swf_count+1);
                                ceres_swf_cc[lat_index]++;
                                //swf_count++;
                //                printf("New avg = %d %f %f\n",swf_count,tlat,avg_swf);
                            }
                        }
                        if(lwf[j]<10000.) 
                        {
                            if(ceres_lwf_cc[lat_index]==-9) {
                                //avg_lwf=lwf[j];
                                ceres_lwf[lat_index]=lwf[j];
                                ceres_lwf_cc[lat_index]=1;
                                //lwf_count=1;
                            }
                            else {
                //                printf("Old avg = %d %f %f\n",swf_count,tlat,avg_swf);
                                ceres_lwf[lat_index] = ((ceres_lwf[lat_index]*
                                    ceres_lwf_cc[lat_index])+lwf[j])/(ceres_lwf_cc[lat_index]+1);
                                //avg_lwf = ((avg_lwf*lwf_count)+lwf[j])/(lwf_count+1);
                                ceres_lwf_cc[lat_index]++;
                                //lwf_count++;
                //                printf("New avg = %d %f %f\n",swf_count,tlat,avg_swf);
                            }
                        }
                    }
                }
                free(lat);
                free(swf);
                free(lwf);

                /////////max_index = pathlen+slen;
                ///////temp_str = (char*)malloc((pathlen+slen)*sizeof(char)+2);
                ///////strcpy(temp_str,data_path);
                ///////strcat(temp_str,de->d_name);
                ///////strcat(temp_str,"/");
                ///////printf("%s\n",temp_str);

                ///////// Search this directory for files
                ///////DIR *dr_2 = opendir(temp_str);
                ///////if(dr_2==NULL)
                ///////{
                ///////    printf("Could not open directory %s\n",temp_str);
                ///////    return 0;   
                ///////}
                ///////    // Loop over all directories 
                //////printf("Average SWF between 55 and 56 N for %s = %f %d\n",tmonth,ceres_swf[145],ceres_swf_cc[145]);
                //////printf("Average LWF between 55 and 56 N for %s = %f %d\n",tmonth,ceres_lwf[145],ceres_lwf_cc[145]);
                free(temp_path);
            }
            free(check_date);
        }
    }
    closedir(dr);

    //printf("Average SWF between 50 and 51 N for %s = %f %d\n",tmonth,avg_swf,swf_count);
    //printf("Average LWF between 50 and 51 N for %s = %f %d\n",tmonth,avg_lwf,lwf_count);

    fout = fopen(outfile_name,"w");
    if( fout != NULL)
    {
        for(i=(minlat+90);i<NLAT;i++)
        {
            lat_index = i-90;
            //printf("%d %f %d %f %d\n",lat_index,ceres_lwf[i],ceres_lwf_cc[i],ceres_swf[i],ceres_swf_cc[i]);
            fprintf(fout,"%d %f %d %f %d\n",lat_index,ceres_lwf[i],ceres_lwf_cc[i],ceres_swf[i],ceres_swf_cc[i]);
        //    fprintf(fout,"%d %d %f %d\n",lat_index,lon_index,misr_aod[i][j],misr_cc[i][j]);
        } 
        fclose(fout);
        printf("Saved file %s\n",outfile_name);
    }
    else
    {
        printf("Error opening %s. Ensure directory tree is intact\n",outfile_name);
    }
   
   
        
    /////char path[]  = "/data/MISR_v23/l5ftl01.larc.nasa.gov/misrl2l3/MISR/MIL2ASAE.003/2006.05.05/MISR_AM1_AS_AEROSOL_P014_O033935_F13_0023.nc";
    ///char data_path[]  = "/data/MISR_v23/l5ftl01.larc.nasa.gov/misrl2l3/MISR/MIL2ASAE.003/";
    /////char average_file[]  = "/home/bsorenson/HighLatitudeStudy/MISR/Updated_Data/30to90/misr_avg_aod_newData.txt";
    ///char meta_group[] = "/4.4_KM_PRODUCTS/";
    /////char data_group[] = "/4.4_KM_PRODUCTS/AUXILIARY/";
    ///char aod_name[]  = "Aerosol_Optical_Depth";
    /////char aod_raw_name[]  = "Aerosol_Optical_Depth_Raw";
    ///char aod_unc_name[]  = "Aerosol_Optical_Depth_Uncertainty_Raw";
    ///char lat_name[]  = "Latitude";
    ///char lon_name[]  = "Longitude";
    ///char flag_name[] = "Aerosol_Retrieval_Screening_Flags";
    ///char avg_file_start[] = "/home/bsorenson/HighLatitudeStudy/MISR/Updated_Data/30to90/aod_no_flags/misr_avg_aod_newData";
    ///char *temp_str,*dir_year,*dir_month,*temp_dir,*file_name,*outfile_name;
    ///char *dir_day;

    ///int i, j, col, row, max_index,slen,pathlen,double_len,s2len;
    ///int lat_index,lon_index;
    ///float tlat,tlon,t_avg;
    ///size_t max_aod,max_lat,max_lon,col_loop,row_loop,index;

    /////char testmonth[] = "01";     // would be passed from outside
    /////char testyear[] = "2004";    // would be passed from outside

    ///struct dirent *de;
    ///struct dirent *de_2;

    ///float misr_aod[NLAT][NLON];     // average AOD array
    ///////float misr_aod_unc[NLAT][NLON];     // average AOD array
    ///int misr_cc[NLAT][NLON];      // count array

    ///FILE *fout;

    /////  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    ///// Generate the output file name   
    ///outfile_name = (char*)malloc((strlen(avg_file_start)+17)*sizeof(char)+1);
    ///////outfile_name = (char*)malloc((strlen(avg_file_start)+11)*sizeof(char)+1);
    ///strcpy(outfile_name,avg_file_start);
    ///strcat(outfile_name,"_");
    ///strcat(outfile_name,argv[1]);
    ///strcat(outfile_name,argv[2]);
    ///strcat(outfile_name,"_noraw.txt");
    ///////strcat(outfile_name,".txt");

    /////col=4032;
    /////row=496;  
    /////max_index=col*row;
   
    ///// Initialize arrays 
    ///for(i=0;i<NLAT;i++)
    ///{
    ///    for(j=0;j<NLON;j++)
    ///    {
    ///        misr_aod[i][j]=-999.;
    ///        ////misr_aod_unc[i][j]=-999.;
    ///        misr_cc[i][j]=-9;
    ///    }
    ///}

    ///// Open the main data storage directory
    ///pathlen = strlen(data_path);
    ///DIR *dr = opendir(data_path);
    ///if(dr==NULL)
    ///{
    ///    printf("Could not open directory\n");
    ///    return 0;   
    ///}
    ///// Loop over all directories 
    ///while((de=readdir(dr)) != NULL)
    ///{
    ///    printf("Current directory %s\n",de->d_name);
    ///    slen = strlen(de->d_name);
    ///    // Ignore '.' and '..'
    ///    if(slen>2)
    ///    {
    ///        // Extract year and month from the directory name
    ///        temp_dir = (char*)malloc(slen*sizeof(char)+1);
    ///        strcpy(temp_dir,de->d_name);
    ///        dir_year = strtok(temp_dir,".");
    ///        dir_month = strtok(NULL,".");
    ///        dir_day = strtok(NULL,"\n");
    ///        if((strcmp(dir_year,argv[1])==0))
    ///        {
    ///            if(strcmp(dir_month,argv[2])==0)
    ///            //if((strcmp(dir_month,argv[2])==0) && (strcmp(dir_day,"14")==0))
    ///            {
    ///                //max_index = pathlen+slen;
    ///                temp_str = (char*)malloc((pathlen+slen)*sizeof(char)+2);
    ///                strcpy(temp_str,data_path);
    ///                strcat(temp_str,de->d_name);
    ///                strcat(temp_str,"/");
    ///                printf("%s\n",temp_str);

    ///                // Search this directory for files
    ///                DIR *dr_2 = opendir(temp_str);
    ///                if(dr_2==NULL)
    ///                {
    ///                    printf("Could not open directory %s\n",temp_str);
    ///                    return 0;   
    ///                }
    ///                // Loop over all directories 
    ///                while((de_2=readdir(dr_2)) != NULL)
    ///                {
    ///                    s2len = strlen(de_2->d_name);
    ///                    //printf("Made it here: %s\n",de_2->d_name);
    ///                    if(s2len>2)
    ///                    {
    ///                        double_len = strlen(temp_str);
    ///                        //max_index = double_len+s2len;
    ///                        file_name = (char*)malloc((double_len+s2len)*sizeof(char)+2);
    ///                        strcpy(file_name,temp_str);
    ///                        strcat(file_name,de_2->d_name);
    ///                        printf("%s\n",file_name);

    ///                        //size_t* unc_dims = (size_t*) get_variable_dims_by_name(file_name, data_group, aod_unc_name);        // NOTE! This is not returning int!!! 
    ///                        //float*  unc_raw  = (float*)  get_variable_data_by_name(file_name, data_group, aod_unc_name);
    ///                        ////size_t* aod_dims = (size_t*) get_variable_dims_by_name(file_name, data_group, aod_raw_name);        // NOTE! This is not returning int!!! 
    ///                        ////float*  aod_raw  = (float*)  get_variable_data_by_name(file_name, data_group, aod_raw_name);
    ///                        size_t* aod_dims = (size_t*) get_variable_dims_by_name(file_name, meta_group, aod_name);        // NOTE! This is not returning int!!! 
    ///                        float*  aod      = (float*)  get_variable_data_by_name(file_name, meta_group, aod_name);
    ///                        size_t* lat_dims = (size_t*) get_variable_dims_by_name(file_name, meta_group, lat_name);        // NOTE! This is not returning int!!! 
    ///                        float*  lat      = (float*)  get_variable_data_by_name(file_name, meta_group, lat_name);
    ///                        size_t* lon_dims = (size_t*) get_variable_dims_by_name(file_name, meta_group, lon_name);        // NOTE! This is not returning int!!! 
    ///                        float*  lon      = (float*)  get_variable_data_by_name(file_name, meta_group, lon_name);

    ///                        // Compare the sizes of each dimension
    ///                        //size_t unc_width  = unc_dims[0];
    ///                        //size_t unc_height = unc_dims[1];
    ///                        //size_t unc_levels = unc_dims[2];
    ///                        size_t aod_width  = aod_dims[0];
    ///                        size_t aod_height = aod_dims[1];
    ///                        size_t aod_levels = aod_dims[2];
    ///                        size_t lat_width  = lat_dims[0];
    ///                        size_t lat_height = lat_dims[1];
    ///                        size_t lat_levels = lat_dims[2];
    ///                        size_t lon_width  = lon_dims[0];
    ///                        size_t lon_height = lon_dims[1];
    ///                        size_t lon_levels = lon_dims[2];

    ///                        //max_unc = unc_width*unc_height;
    ///                        max_aod = aod_width*aod_height;
    ///                        max_lat = lat_width*lat_height;
    ///                        max_lon = lon_width*lon_height;
    ///            
    ///                        if((max_aod!=max_lat) | (max_aod!=max_lon) | (max_lat!=max_lon))
    ///                        {
    ///                            printf("Different sizes. aod=%d lat=%d lon=%d\n",max_aod,max_lat,max_lon);
    ///                        }

    ///                        // Read data from this filename
    ///                        // Get the data 
    ///                        ///float* aod_raw = (float*) get_variable_data_by_name(file_name, data_group,aod_raw_name);
    ///                        //////float* aod_unc_raw = (float*) get_variable_data_by_name(path, data_group,  \
    ///                        //////                 aod_raw_unc_name);
    ///                        ///float* lat = (float*) get_variable_data_by_name(file_name, meta_group, lat_name);
    ///                        ///float* lon = (float*) get_variable_data_by_name(file_name, meta_group, lon_name);
    ///                        /////printf("after lon\n");
    ///                        /////unsigned int* flag = (unsigned int*) get_variable_data_by_name(path, meta_group, flag_name);
    ///
    ///                        ////// Read in data from average file
    ///                        ////read_misr(argv[2],misr_aod,misr_cc);
    ///                        ////printf("From file, %f %d\n",misr_aod[22][72],misr_cc[22][72]);

    ///                        //printf("\nUpdating average MISR AOD with new data\n\n"); 
    ///                        //// loop over data and  
    ///                        //for(i=0;i<max_index;i++)
    ///                        for (row_loop = 0; row_loop <aod_height; row_loop++) {
    ///                            for (col_loop = 0; col_loop < aod_width; col_loop++) {
    ///                                index = row_loop * aod_width + col_loop;
    ///                                // Check QC flag
    ///                                //
    ///                                // if(QC flag is met)
    ///                                // {
    ///                                //printf("Comparing filenames\n");
    ///                                //if(strcmp(file_name,"/data/MISR_v23/l5ftl01.larc.nasa.gov/misrl2l3/MISR/MIL2ASAE.003/2000.03.14/MISR_AM1_AS_AEROSOL_P009_O001272_F13_0023.nc")==0)
    ///                                //{
    ///                                //    //printf("About to print aod\n");
    ///                                //    printf("%d aod = %f\n",i,aod_raw[i]);
    ///                                //} 
    ///                                if(aod[index]>-900.0)
    ///                                ////if(aod_raw[index]>-900.0)
    ///                                //if(aod_raw[i]>-900.0)
    ///                                {
    ///                                    lat_index = (int)lat[index]+90;
    ///                                    lon_index = (int)lon[index]+180;
    ///                                    if((lat[index]>-200) & (lon[index]>-300))
    ///                                    {
    ///                                        // Update AOD values
    ///                                        if(misr_aod[lat_index][lon_index]<-100)
    ///                                        {
    ///                                            misr_aod[lat_index][lon_index]=aod[index];
    ///                                            ////misr_aod[lat_index][lon_index]=aod_raw[index];
    ///                                            //misr_aod_unc[lat_index][lon_index]=unc_raw[index];
    ///                                            misr_cc[lat_index][lon_index]=1;
    ///                                        }
    ///                                        else
    ///                                        {        
    ///                                            misr_aod[lat_index][lon_index] = (float)((misr_aod[lat_index][lon_index]*
    ///                                                (float)misr_cc[lat_index][lon_index])+aod[index])/((float)(misr_cc[lat_index][lon_index]+1));
    ///                                            ////misr_aod[lat_index][lon_index] = (float)((misr_aod[lat_index][lon_index]*
    ///                                            ////    (float)misr_cc[lat_index][lon_index])+aod_raw[index])/((float)(misr_cc[lat_index][lon_index]+1));
    ///                                            //misr_aod_unc[lat_index][lon_index] = (float)(sqrt(pow(misr_aod_unc[lat_index][lon_index]*
    ///                                            //    misr_cc[lat_index][lon_index],2.0)+pow(unc_raw[index],2.0)))/((float)(misr_cc[lat_index][lon_index]+1));
    ///                                            //misr_aod_unc[lat_index][lon_index] = (float)((misr_aod_unc[lat_index][lon_index]*
    ///                                            //    (float)misr_cc[lat_index][lon_index])+unc_raw[index])/((float)(misr_cc[lat_index][lon_index]+1));
    ///                                            misr_cc[lat_index][lon_index]++;
    ///                                        }
    ///                                        //// Update counts
    ///                                        //if(misr_cc[lat_index][lon_index]==-9)
    ///                                        //{
    ///                                        //    misr_cc[lat_index][lon_index]=1;
    ///                                        //}
    ///                                        //else
    ///                                        //{        
    ///                                        //    misr_cc[lat_index][lon_index]++;
    ///                                        //}
    ///                                        //printf("%d %d %f %d\n",lat_index,lon_index,misr_aod[lat_index][lon_index],misr_cc[lat_index][lon_index]);
    ///                                    }
    ///                                }
    ///                                // }
    ///                            }
    ///                        }
    ///                        //printf("Here?\n");
    ///                        //printf("%f %d\n",misr_aod[74][188],misr_cc[74][188]);
    ///                        // Free up data
    ///                        free(file_name);
    ///                        ////free(aod_raw);
    ///                        free(aod);
    ///                        //free(unc_raw);
    ///                        free(lat);
    ///                        free(lon);
    ///                        free(aod_dims);
    ///                        //free(unc_dims);
    ///                        free(lat_dims);
    ///                        free(lon_dims);
    ///                    }
    ///                }
    ///                free(temp_str);
    ///                closedir(dr_2);
    ///            }
    ///        }
    ///        free(temp_dir);
    ///    }
    ///}
    ///closedir(dr);

    ///printf("Done reading\n");
    /////printf("%f %d\n",misr_aod[37][82],misr_cc[37][82]);
    ///// Open output file
//  ///  fout = fopen(average_file,"w");
   
    ///
    ///// Generate the output file name
    ///if( (fout=fopen(outfile_name,"w")) != NULL)
    ///{
    ///   // fprintf(fout,"LAT   LON    AVG_AOD  AVG_UNC   AOD_COUNT\n"); 
    ///    for(i=0;i<NLAT;i++)
    ///    {
    ///        for(j=0;j<NLON;j++)
    ///        {
    ///            lat_index = i-90;
    ///            lon_index = j-180;
    ///            fprintf(fout,"%d %d %f %d\n",lat_index,lon_index,misr_aod[i][j],misr_cc[i][j]);
    ///            ////fprintf(fout,"%d %d %f %f %d\n",lat_index,lon_index,misr_aod[i][j],misr_aod_unc[i][j],misr_cc[i][j]);
    ///        }
    ///    } 
    ///}
    ///printf("Done writing to file %s\n",outfile_name);
    ///fclose(fout);
    ///free(outfile_name);
    ///temp_str=NULL;
    ///temp_dir=NULL;
    /////printf("Args = %s %s\n",argv[1],argv[2]);
    return 0;
   
    ////////// Get the data 
    ////////float* aod_raw = (float*) get_variable_data_by_name(argv[1], data_group,aod_raw_name);
    ///////////float* aod_unc_raw = (float*) get_variable_data_by_name(path, data_group,  \
    ///////////                 aod_raw_unc_name);
    ////////float* lat = (float*) get_variable_data_by_name(argv[1], meta_group, lat_name);
    ////////float* lon = (float*) get_variable_data_by_name(argv[1], meta_group, lon_name);
    /////////////unsigned int* flag = (unsigned int*) get_variable_data_by_name(path, meta_group, flag_name);
    //////////
    //////////for(i=0;i<argc;i++)
    //////////{
    //////////    printf("%s\n",argv[i]);
    //////////} 

    ////////printf("Processing %s\n",argv[1]);
    ////////
    ////////
    ////////col=4032;
    ////////row=496;  
    ////////max_index=col*row;

    ////////// Read in data from average file
    ////////read_misr(argv[2],misr_aod,misr_cc);
    //////////printf("From file, %f %d\n",misr_aod[22][72],misr_cc[22][72]);

    //////////printf("\nUpdating average MISR AOD with new data\n\n"); 
    //////////// loop over data and  
    ////////for(i=0;i<max_index;i++)
    ////////{
    ////////    // Check QC flag
    ////////    //
    ////////    // if(QC flag is met)
    ////////    // {
    ////////   
    ////////    if(aod_raw[i]>-900.0)
    ////////    {
    ////////        lat_index = (int)lat[i]+90;
    ////////        lon_index = (int)lon[i]+180;
    ////////        if((lat[i]>-200) & (lon[i]>-300))
    ////////        {
    ////////            // Update AOD values
    ////////            if(misr_aod[lat_index][lon_index]<-100)
    ////////            {
    ////////                misr_aod[lat_index][lon_index]=aod_raw[i];
    ////////                misr_cc[lat_index][lon_index]=1;
    ////////                //if((lat_index==22) & (lon_index==72))
    ////////                //{
    ////////                //    printf("First instance, %f %d\n",misr_aod[lat_index][lon_index],misr_cc[lat_index][lon_index]);
    ////////                //}
    ////////            }
    ////////            else
    ////////            {        
    ////////                misr_aod[lat_index][lon_index] = (float)((misr_aod[lat_index][lon_index]*
    ////////                    (float)misr_cc[lat_index][lon_index])+aod_raw[i])/((float)(misr_cc[lat_index][lon_index]+1));
    ////////                misr_cc[lat_index][lon_index]++;
    ////////                //if((lat_index==22) & (lon_index==72))
    ////////                //{
    ////////                //    printf("In new data, aod=%f\n",aod_raw[i]);
    ////////                //    printf("In new data, avg_aod=%f\n",misr_aod[lat_index][lon_index]);
    ////////                //}
    ////////            }
    ////////            //// Update counts
    ////////            //if(misr_cc[lat_index][lon_index]==-9)
    ////////            //{
    ////////            //    misr_cc[lat_index][lon_index]=1;
    ////////            //}
    ////////            //else
    ////////            //{        
    ////////            //    misr_cc[lat_index][lon_index]++;
    ////////            //}
    ////////            //printf("%d %d %f %d\n",lat_index,lon_index,misr_aod[lat_index][lon_index],misr_cc[lat_index][lon_index]);
    ////////        }
    ////////    }
    ////////    // }
    ////////} 
//  ////////  fout = fopen(average_file,"w");
   
    //////////fprintf(fout,"LAT   LON    AVG_AOD     AOD_COUNT\n"); 
    ////////
    ////// Generate the output file name
    ////if( (fout=fopen(argv[2],"w")) != NULL)
    //////if( (fout=fopen(outfile_name,"w")) != NULL)
    ////{
    ////    for(i=0;i<NLAT;i++)
    ////    {
    ////        for(j=0;j<NLON;j++)
    ////        {
    ////            lat_index = i-90;
    ////            lon_index = j-180;
    ////            fprintf(fout,"%d %d %f %d\n",lat_index,lon_index,misr_aod[i][j],misr_cc[i][j]);
    ////            //if(misr_cc[i][j]>-9)
    ////            //{
    ////            //    //t_avg = (float)misr_aod[i][j]/misr_cc[i][j];
    ////            //    lat_index = i-90;
    ////            //    lon_index = j-180;
    ////            //    //printf("%d %d %f %f %d\n",i,j,misr_aod[i][j],t_avg,misr_cc[i][j]);
    ////            //    fprintf(fout,"%d %d %f %d\n",lat_index,lon_index,t_avg,misr_cc[i][j]);
    ////            //}
    ////            //else
    ////            //{
    ////            //    fprintf(fout,"%d %d %f %d\n",lat_index,lon_index,misr_aod[i][j],misr_cc[i][j]);
    ////            //}
    ////        }
    ////    } 
    ////    fclose(fout);
    ////}
    ////else
    ////{
    ////    printf("Error opening output file\n");
    ////}

    //printf("Done writing\n");
    //for (row = 0; row < 496; row++) {
    //    //for (col = 0; col < 4032; col++) {
    //        col = 0;
    //        int index = row * 4032 + col;
    //        if (data[index] != -9999) printf("%d [%d, %d]-> %5.3f\n", index, index/4032, index%4032, data[index]);
    //    //}
    //    printf("\n");
    //}
    //
    //printf("data at 4032*160+450 = %f\n", data[4032*160+450]);
    
    ////// Free up data
    ////free(aod_raw);
    //////free(aod_unc_raw);
    ////free(lat);
    ////free(lon);
    ////return 0;
}
