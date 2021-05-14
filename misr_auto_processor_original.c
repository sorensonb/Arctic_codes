/*----------------------------------------------------------------------------/
 *
 *
 * SYNTAX:
 *  $ gcc -o misr_exec netcdf_helper_newest.c -I/usr/include/ -lnetcdf
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
#include "/usr/include/netcdf.h"

#define NLON  360
#define NLAT  180

#define DEBUG 1
    
void handle_error(const char* err_name, int err_code) {
    printf("Error: %s = %d\n", err_name, err_code);
    exit(1);
}

void* get_variable_data_by_name(const char *path, const char* grouppath, const char *name) {

    // Open the file by name
    //int nc_open (const char * path, int omode, int * ncidp) 

    int ncid = 0;
    int err_open = nc_open(path, NC_NOWRITE, &ncid);
    if (err_open != NC_NOERR) handle_error("err_open", err_open);

    // Determine group hierarchy based on full variable path name
    int grp_ncid;
    int err_group = nc_inq_grp_full_ncid(ncid, grouppath, &grp_ncid);
    if (err_group != NC_NOERR) handle_error("err_group", err_group);

    int nvars = 0;
    int varids[100] = {0};
    int err_nvars = nc_inq_varids(grp_ncid, &nvars, varids);
    if (err_nvars != NC_NOERR) handle_error("err_nvars", err_nvars);
     
    //printf("%s contains %d variables\n", grouppath, nvars);
    
    if (DEBUG) {
       int var_i;
       for (var_i = 0; var_i < nvars; var_i++) {
           char name[NC_MAX_NAME + 1];
           nc_inq_varname(grp_ncid, varids[var_i], name);
           printf("%3d\t%3d\t%50s\n", var_i, varids[var_i], name);
       }
    }

    // Retrieve the varid by the variable name
    //int nc_inq_varid(int ncid, const char * name, int * varidp)

    int varid = 0;
    int err_var = nc_inq_varid(grp_ncid, name, &varid);
    if (err_var != NC_NOERR) handle_error("err_var", err_var);

    // Retrieve the metadata for the variable to learn
    // ... dimensions (need all of the dimension ids)
    //int nc_inq_varndims(int ncid, int varid, int * ndimsp) 
    //int nc_inq_vardimid(int ncid, int varid, int * dimids)

    int ndims = 1;
    int err_dims = nc_inq_varndims(grp_ncid, varid, &ndims);
    if (err_dims != NC_NOERR) handle_error("err_dims", err_dims);

    int *dimids = malloc(sizeof(int) * ndims);
    int err_dimi = nc_inq_vardimid(grp_ncid, varid, dimids);
    if (err_dimi != NC_NOERR) handle_error("err_dimi", err_dimi);

    // ... length of each dimension
    // ... total data count
    size_t total_data_count = 1;
    int dim_i;
    for (dim_i = 0; dim_i < ndims; dim_i++) {
        //int nc_inq_dim(int ncid, int dimid, char *name, int *len);
        size_t len = 1;
        int err_len = nc_inq_dim(grp_ncid, dimids[dim_i], NULL, &len);
        if (len == 0) printf("Length was ZERO!\n");
        if (err_len != NC_NOERR) handle_error("err_len", err_len);
        total_data_count *= len;
    }

    // ... type
    //int nc_inq_vartype(int ncid, int varid, nc_type *typep)

    int type = 0;
    int err_typ = nc_inq_vartype(grp_ncid, varid, &type);
    if (err_typ != NC_NOERR) handle_error("err_typ", err_typ);

    /* From netcdf.h 
     #define NC_NAT         0         ??? not a type ???
     #define NC_BYTE         1         byte        >>> unsigned char
     #define NC_CHAR         2         char
     #define NC_SHORT         3         short
     #define NC_INT         4         int
     #define NC_LONG         NC_INT     int
     #define NC_FLOAT         5         float
     #define NC_DOUBLE         6         double
     #define NC_UBYTE         7         unsigned byte
     #define NC_USHORT         8         unsigned short
     #define NC_UINT         9         unsigned int
     #define NC_INT64         10         long long
     #define NC_UINT64         11         unsigned long long
     #define NC_STRING         12        char**
     */

    // Get size based on the returned type
    int type_size = 0;
    switch (type) {
        case NC_BYTE :
            type_size = sizeof(char);
            break;
        case NC_CHAR :
            type_size = sizeof(char);
            break;
        case NC_SHORT :
            type_size = sizeof(short);
            break;
        case NC_INT :            // & NC_LONG (same value in header (4), see above)
            type_size = sizeof(int);
            break;
        case NC_FLOAT :
            type_size = sizeof(float);
            break;
        case NC_DOUBLE :
            type_size = sizeof(double);
            break;
        case NC_UBYTE :
            type_size = sizeof(unsigned char);
            break;
        case NC_USHORT :
            type_size = sizeof(unsigned short);
            break;
        case NC_UINT :
            type_size = sizeof(unsigned int);
            break;
        case NC_INT64 :
            type_size = sizeof(long long);
            break;
        case NC_UINT64 :
            type_size = sizeof(unsigned long long);
            break;
        case NC_STRING :
            type_size = sizeof(char*);
            break;
    }

    // .. end-edness
    //int int nc_inq_var_endian(int ncid, int varid, int *endian)

    /*        >>>> Should already be handeled within the library.
    int endian = 0;
    int err_end = nc_inq_var_endian(ncid, varid, &endian);

     #define NC_ENDIAN_LITTLE
     #define NC_ENDIAN_BIG
     #define NC_ENDIAN_NATIVE
     */
    
    // .. fill value
    // TODO
    // Types may be an issue
    /*
    float fill_value = 1;
    int err_scale = nc_get_att(ncid, varid, _FillValue, (void*) scale_factor);
     */
     
    // .. scale factor
    /* TODO
    float scale_factor = 1;
    int err_scale = nc_get_att(ncid, varid, "scale_factor", (void*) scale_factor);
     */
     
    // .. offset
    /* TODO
    float offset_factor = 0;
    int err_offset = nc_get_att(ncid, varid, "offset", (void*) offset);
     */
    
    if (DEBUG) {
        printf("%d %d\n", type_size, total_data_count);
    }
    // Retrieve the data using the varid
    //int nc_get_var(int ncid, int varid, void* ip)
    void* data = malloc(type_size * total_data_count);
    int err_data = nc_get_var(grp_ncid, varid, data);
    if (err_data != NC_NOERR) handle_error("err_data", err_data);
    
    // close the file
    // int nc_close(int ncid)
    int err_close = nc_close(ncid);
    if (err_close != NC_NOERR) handle_error("err_close", err_close);
    
    // .. do needed conversions to 'make it right'
    // TODO
        
    return (void*) data;
}

void read_misr(char *filename, float misr_aod[NLAT][NLON], int misr_cc[NLAT][NLON])
{
    FILE *fp;
    int i,j,tlat,tlon,tempcc;
    float tempaod;
   
    // If misr average data file is not already created, fill array with 
    // missing values. 
    if( (fp = fopen(filename,"r")) == NULL)
    {
        printf("File %s does not exist. Initializing array\n",filename);
        for(i=0;i<NLAT;i++)
        {
            for(j=0;j<NLON;j++)
            {
                misr_aod[i][j]=-999.;
                misr_cc[i][j]=-9;
            }
        }
    }
    // If file does exist, read its contents into the array
    else
    {
        printf("Reading data from file %s\n",filename);
        for(i=0;i<NLAT;i++)
        {
            for(j=0;j<NLON;j++)
            {
                fscanf(fp,"%d %d %f %d",&tlat,&tlon,&misr_aod[i][j],&misr_cc[i][j]);
            }
        }
        fclose(fp);
    }
}

int main(int argc, char *argv[]) {
   
    // Check command line arguments
    if(argc != 3)
    {
        //printf("SYNTAX: ./misr_exec data_file_name average_file\n");
        printf("SYNTAX: ./misr_exec year month\n");
        return 0;
    }

    //char path[]  = "/data/MISR_v23/l5ftl01.larc.nasa.gov/misrl2l3/MISR/MIL2ASAE.003/2006.05.05/MISR_AM1_AS_AEROSOL_P014_O033935_F13_0023.nc";
    char data_path[]  = "/data/MISR_v23/l5ftl01.larc.nasa.gov/misrl2l3/MISR/MIL2ASAE.003/";
    //char average_file[]  = "/home/bsorenson/HighLatitudeStudy/MISR/Updated_Data/30to90/misr_avg_aod_newData.txt";
    char meta_group[] = "/4.4_KM_PRODUCTS/";
    char data_group[] = "/4.4_KM_PRODUCTS/AUXILIARY/";
    char aod_raw_name[]  = "Aerosol_Optical_Depth_Raw";
    char aod_raw_unc_name[]  = "Aerosol_Optical_Depth_Uncertainty_Raw";
    char lat_name[]  = "Latitude";
    char lon_name[]  = "Longitude";
    char flag_name[] = "Aerosol_Retrieval_Screening_Flags";
    char avg_file_start[] = "/home/bsorenson/HighLatitudeStudy/MISR/Updated_Data/30to90/misr_avg_aod_newData";
    char *temp_str,*dir_year,*dir_month,*temp_dir,*file_name,*outfile_name;
    char *dir_day;

    int i, j, col, row, max_index,slen,pathlen,double_len,s2len;
    int lat_index,lon_index;
    float tlat,tlon,t_avg;

    //char testmonth[] = "01";     // would be passed from outside
    //char testyear[] = "2004";    // would be passed from outside

    struct dirent *de;
    struct dirent *de_2;

    float misr_aod[NLAT][NLON];     // average AOD array
    int misr_cc[NLAT][NLON];      // count array

    FILE *fout;

    //  = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    // Generate the output file name   
    outfile_name = (char*)malloc((strlen(avg_file_start)+11)*sizeof(char)+1);
    strcpy(outfile_name,avg_file_start);
    strcat(outfile_name,"_");
    strcat(outfile_name,argv[1]);
    strcat(outfile_name,argv[2]);
    strcat(outfile_name,".txt");

    col=4032;
    row=496;  
    max_index=col*row;
   
    // Initialize arrays 
    for(i=0;i<NLAT;i++)
    {
        for(j=0;j<NLON;j++)
        {
            misr_aod[i][j]=-999.;
            misr_cc[i][j]=-9;
        }
    }

    // Open the main data storage directory
    pathlen = strlen(data_path);
    DIR *dr = opendir(data_path);
    if(dr==NULL)
    {
        printf("Could not open directory\n");
        return 0;   
    }
    // Loop over all directories 
    while((de=readdir(dr)) != NULL)
    {
        printf("Current directory %s\n",de->d_name);
        slen = strlen(de->d_name);
        // Ignore '.' and '..'
        if(slen>2)
        {
            // Extract year and month from the directory name
            temp_dir = (char*)malloc(slen*sizeof(char)+1);
            strcpy(temp_dir,de->d_name);
            dir_year = strtok(temp_dir,".");
            dir_month = strtok(NULL,".");
            dir_day = strtok(NULL,"\n");
            if((strcmp(dir_year,argv[1])==0))
            {
                if((strcmp(dir_month,argv[2])==0) && (strcmp(dir_day,"14")==0))
                {
                    //max_index = pathlen+slen;
                    temp_str = (char*)malloc((pathlen+slen)*sizeof(char)+2);
                    strcpy(temp_str,data_path);
                    strcat(temp_str,de->d_name);
                    strcat(temp_str,"/");
                    printf("%s\n",temp_str);

                    // Search this directory for files
                    DIR *dr_2 = opendir(temp_str);
                    if(dr_2==NULL)
                    {
                        printf("Could not open directory %s\n",temp_str);
                        return 0;   
                    }
                    // Loop over all directories 
                    while((de_2=readdir(dr_2)) != NULL)
                    {
                        s2len = strlen(de_2->d_name);
                        //printf("Made it here: %s\n",de_2->d_name);
                        if(s2len>2)
                        {
                            double_len = strlen(temp_str);
                            //max_index = double_len+s2len;
                            file_name = (char*)malloc((double_len+s2len)*sizeof(char)+2);
                            strcpy(file_name,temp_str);
                            strcat(file_name,de_2->d_name);
                            printf("%s\n",file_name);

                            // Read data from this filename
                            // Get the data 
                            float* aod_raw = (float*) get_variable_data_by_name(file_name, data_group,aod_raw_name);
                            ///float* aod_unc_raw = (float*) get_variable_data_by_name(path, data_group,  \
                            ///                 aod_raw_unc_name);
                            float* lat = (float*) get_variable_data_by_name(file_name, meta_group, lat_name);
                            float* lon = (float*) get_variable_data_by_name(file_name, meta_group, lon_name);
                            //printf("after lon\n");
                            /////unsigned int* flag = (unsigned int*) get_variable_data_by_name(path, meta_group, flag_name);
    
                            ////// Read in data from average file
                            ////read_misr(argv[2],misr_aod,misr_cc);
                            ////printf("From file, %f %d\n",misr_aod[22][72],misr_cc[22][72]);

                            //printf("\nUpdating average MISR AOD with new data\n\n"); 
                            //// loop over data and  
                            for(i=0;i<max_index;i++)
                            {
                                // Check QC flag
                                //
                                // if(QC flag is met)
                                // {
                                //printf("Comparing filenames\n");
                                if(strcmp(file_name,"/data/MISR_v23/l5ftl01.larc.nasa.gov/misrl2l3/MISR/MIL2ASAE.003/2000.03.14/MISR_AM1_AS_AEROSOL_P009_O001272_F13_0023.nc")==0)
                                {
                                    //printf("About to print aod\n");
                                    printf("%d aod = %f\n",i,aod_raw[i]);
                                } 
                                if(aod_raw[i]>-900.0)
                                {
                                    lat_index = (int)lat[i]+90;
                                    lon_index = (int)lon[i]+180;
                                    if((lat[i]>-200) & (lon[i]>-300))
                                    {
                                        // Update AOD values
                                        if(misr_aod[lat_index][lon_index]<-100)
                                        {
                                            misr_aod[lat_index][lon_index]=aod_raw[i];
                                            misr_cc[lat_index][lon_index]=1;
                                        }
                                        else
                                        {        
                                            misr_aod[lat_index][lon_index] = (float)((misr_aod[lat_index][lon_index]*
                                                (float)misr_cc[lat_index][lon_index])+aod_raw[i])/((float)(misr_cc[lat_index][lon_index]+1));
                                            misr_cc[lat_index][lon_index]++;
                                        }
                                        //// Update counts
                                        //if(misr_cc[lat_index][lon_index]==-9)
                                        //{
                                        //    misr_cc[lat_index][lon_index]=1;
                                        //}
                                        //else
                                        //{        
                                        //    misr_cc[lat_index][lon_index]++;
                                        //}
                                        //printf("%d %d %f %d\n",lat_index,lon_index,misr_aod[lat_index][lon_index],misr_cc[lat_index][lon_index]);
                                    }
                                }
                                // }
                            } 
                            //printf("Here?\n");
                            //printf("%f %d\n",misr_aod[74][188],misr_cc[74][188]);
                            // Free up data
                            free(aod_raw);
                            free(file_name);
                            //free(aod_unc_raw);
                            free(lat);
                            free(lon);
                        }
                    }
                    free(temp_str);
                    closedir(dr_2);
                }
            }
            free(temp_dir);
        }
    }
    closedir(dr);

    //printf("Done reading\n");
    //printf("%f %d\n",misr_aod[37][82],misr_cc[37][82]);
    // Open output file
//    fout = fopen(average_file,"w");
   
    //fprintf(fout,"LAT   LON    AVG_AOD     AOD_COUNT\n"); 
    
    ///// Generate the output file name
    ///if( (fout=fopen(outfile_name,"w")) != NULL)
    ///{
    ///    for(i=0;i<NLAT;i++)
    ///    {
    ///        for(j=0;j<NLON;j++)
    ///        {
    ///            lat_index = i-90;
    ///            lon_index = j-180;
    ///            fprintf(fout,"%d %d %f %d\n",lat_index,lon_index,misr_aod[i][j],misr_cc[i][j]);
    ///            //if(misr_cc[i][j]>-9)
    ///            //{
    ///            //    //t_avg = (float)misr_aod[i][j]/misr_cc[i][j];
    ///            //    lat_index = i-90;
    ///            //    lon_index = j-180;
    ///            //    //printf("%d %d %f %f %d\n",i,j,misr_aod[i][j],t_avg,misr_cc[i][j]);
    ///            //    fprintf(fout,"%d %d %f %d\n",lat_index,lon_index,t_avg,misr_cc[i][j]);
    ///            //}
    ///            //else
    ///            //{
    ///            //    fprintf(fout,"%d %d %f %d\n",lat_index,lon_index,misr_aod[i][j],misr_cc[i][j]);
    ///            //}
    ///        }
    ///    } 
    ///}
    ///fclose(fout);
    ///free(outfile_name);
    temp_str=NULL;
    temp_dir=NULL;
    //printf("Args = %s %s\n",argv[1],argv[2]);
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
