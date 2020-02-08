;
; NAME:
;   Average_AerosolIndexCalculator.pro
;
; PURPOSE:
;   Calcualte the average AI value for each 1x1 degree lat/lon grid in the 
;   northern hemisphere from 48 degrees N to 90 degrees N. The averages
;   are printed to an output file at the end of the code.
;
; SYNTAX:
;  `
;   IDL> Average_AerosolIndexCalculator date first
;
;       date:  A string containing at least 'CCYYMMDD' must be passed
;       first: integer (0,1) 
;           0 => Append to the end of the file and do not print a header line
;           1 => Overwrite the previous file and write a header line
;
; EXAMPLE:
;   $ idl
;   IDL> .run Average_AerosolIndexCalculator.pro
;   % Compiled module: Calc_Average.
;   IDL> Calc_Average, '20140813', 0
; 
;   This calculates the average AI values for all files for August 13, 2014 
;
; MODIFICATIONS:
;   Blake Sorenson <blake.sorenson@und.edu>     - 2018/09/18:
;     Written (modification of AerosolIndexPlotting_Procedure_August112014t0046.pro by rcontreras)
;
;



;creating a procedure in order to plot figure 1 from the paper
;similar to the code of creating a single Aerosol Index plot however will now create several plots in order to create a single .EPS file
;attempting to save it as a ps file
Pro Average_AerosolIndexCalculator, date, first


;define the directory to the aerosol index files
homepath    = '/home/bsorenson/HighLatitudeStudy/OMI/Updated_Data/30to90/flag_test/'
pathlength  = '/home/shared/OMAERUV/OMAERUV_Parameters_20190212_download/Aerosol_Index/'
pathlength1 = '/home/shared/OMAERUV/OMAERUV_Parameters_20190212_download/Latitude/'
pathlength2 = '/home/shared/OMAERUV/OMAERUV_Parameters_20190212_download/Longitude/'
pathlength3 = '/home/shared/OMAERUV/OMAERUV_Parameters_20190212_download/X_Track_Quality_Flags/'
pathlengthM = '/home/shared/OMAERUV/OMAERUV_Parameters_20190212_download/Measurement_Quality_Flags/'
pathlengthP = '/home/shared/OMAERUV/OMAERUV_Parameters_20190212_download/Pixel_Quality_Flags/'
pathlengthG = '/home/shared/OMAERUV/OMAERUV_Parameters_20190212_download/Ground_Pixel_Quality_Flags/'


np = 1440
nl = 720

sdate = strtrim(string(date),1)
    
; Create arrays to hold the latitude and longitude ranges
lat_ranges = [30.0:90.0:1.0]
lon_ranges = [-180.0:180.0:1.0]

; Loop over latitude and longtiude ranges and average the AI data
num_lats = n_elements(lat_ranges)
num_lons = n_elements(lon_ranges)

avg_ai = dblarr(num_lats,num_lons)
std_ai = dblarr(num_lats,num_lons)
num_ai = intarr(num_lats,num_lons)


;search for the number of files which will be averaging for the current time frame
;date=strmid(dateandtime,0,8)
ailist             = file_search(pathlength+'UVAerosolIndex_'+sdate+'*.bin')
latlist            = file_search(pathlength1+'Latitude_'+ sdate+'*.bin')
lonlist            = file_search(pathlength2+'Longitude_'+ sdate+'*.bin')
xtracklist         = file_search(pathlength3+'XTrackQualityFlags_'+sdate+'*.bin')
msmnt_qc_list      = file_search(pathlengthM+'MeasurementQualityFlags_'+sdate+'*.bin')
pixel_qc_list      = file_search(pathlengthP+'PixelQualityFlags_'+sdate+'*.bin')
grnd_pixel_qc_list = file_search(pathlengthG+'GroundPixelQualityFlags_'+sdate+'*.bin')
num_files = n_elements(ailist)

; Find the output file name depending on the length of the date.
; If the date is in the format of CCYYMM, then write the data to a monthly
; averaged output file. If the format is CCYYMMDD, then write the data to
; a daily averaged output file.
if(strlen(sdate) eq 8) then outname = 'omi_daily_average_AI_xtrack_row023_newData.txt'
if(strlen(sdate) eq 6) then outname = 'omi_monthly_average_AI_xtrack_row023_newData.txt'

;outname = 'testidlout.txt'

; Open the output file
if( first eq 1) then begin
    OPENW, 1, homepath+outname
    printf, 1, 'Date,LatxLon,Avg,#_obs'
endif else begin
    OPENW, 1, homepath+outname, /APPEND
endelse
; -----------------------------------------------------------------------------
;
; Begin file loop
;
; -----------------------------------------------------------------------------
for fi = 0, num_files-1 do begin 
    
    AI = READ_BINARY(ailist[fi],DATA_TYPE = 4)
    
    ;read in the binary files keeping track of columns and rows
    CBA = SIZE(AI, /DIMENSIONS)
    CBA2 = CBA/60
    
    
    AI = FLTARR(60,CBA2)
    print, ailist[fi]
    OPENR,2,ailist[fi]
    READU,2,AI
    close,2
    
    LAT = FLTARR(60,CBA2)
    OPENR,2,latlist[fi]
    READU,2,LAT
    close,2
    
    
    LON = FLTARR(60,CBA2)
    OPENR,2,lonlist[fi]
    READU,2,LON
    close,2
    
    
    ;read in the cross track flags parameters
    Xtrack = BYTARR(60,CBA2)
    OPENR,2,xtracklist[fi]
    READU,2,XTrack
    close,2

    ;read in the measurement qc flags parameters
    MsmntQC = BYTARR(60,CBA2)
    OPENR,2,msmnt_qc_list[fi]
    READU,2,MsmntQC
    close,2
   
    ;read in the pixel qc flags parameters
    PixelQC = BYTARR(60,CBA2)
    OPENR,2,pixel_qc_list[fi]
    READU,2,PixelQC
    close,2
   
    ;read in the pixel qc flags parameters
    GrndPixelQC = BYTARR(60,CBA2)
    OPENR,2,grnd_pixel_qc_list[fi]
    READU,2,GrndPixelQC
    close,2
   
    ; Reshape all arrays to only use the first 23 rows
    AI  = AI[0:23,*]
    LAT = LAT[0:23,*]
    LON = LON[0:23,*]
    XTrack  = XTrack[0:23,*]
    MsmntQC = MsmntQC[0:23,*]
    PixelQC = PixelQC[0:23,*]
    GrndPixelQC = GrndPixelQC[0:23,*]
 
    
    ; = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    ;
    ; Old averaging code from original program.
    ;
    ; = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    ;;;UVAI = fltarr(np,nl)
    ;;;count = fltarr(np,nl)
    ;;;count[*,*] = 0
    ;;;
    ;;;for ii = 0,59 do begin
    ;;;    for ii2 = 0,cba2[0]-1 do begin
    ;;;        if ai[ii,ii2] gt -20 then begin
    ;;;        
    ;;;            ;only enter this loop if our cross track quality flags meet certain requirements
    ;;;            if xtrack[ii,ii2] eq 0 or (xtrack[ii,ii2] and 4) eq 4 then begin
    ;;;                index1 = floor(Lat(ii,ii2)*4 + 360.)
    ;;;                index2 = floor(Lon(ii,ii2)*4 + 720.)
    ;;;                
    ;;;                if(index1 lt 0) then index1 = 0
    ;;;                if(index1 ge 719) then index1 = 719
    ;;;                if(index2 lt 0) then index2 = 0
    ;;;                if(index2 ge 1439) then index2 = 1439
    ;;;                
    ;;;                UVAI(index2, index1) = (UVAI(index2,index1)*count(index2,index1) + ai(ii,ii2))/(count(index2,index1)+1l)
    ;;;                ;if((ii eq 59) and (ii2 eq 1309)) then print,UVAI(index2,index1)
    ;;;                count(index2, index1) = count(index2,index1) + 1l
    ;;;            
    ;;;            endif
    ;;;        endif
    ;;;    endfor
    ;;;endfor
   
    ;;;aa = where(count eq 0,ctaa)
    ;;;if ctaa gt 0 then begin 
    ;;;    UVAI[aa] = !VALUES.F_NAN
    ;;;endif

    for i =0,num_lats-2 do begin
        for j =0,num_lons-2 do begin
            ; Find the AI data within the current 1x1 degree lat/lon box
            ;ai_in_coord = AI(where( ((LAT gt lat_ranges(i)) && (LAT lt lat_ranges(i+1))) && ((LON gt lat_ranges(j)) && (LON lt lon_ranges(j+1))) ))

            ; Make sure XTrack data are good
            ;;;; Now find where xtrack are good 
            ;;;ai_in_coord = AI(where(Xcheck eq 0 or (Xcheck and 4) eq 4))
            ;;;print, ai_in_coord, xtrack(where(AI eq ai_in_coord))
            ai_indices = where(((LAT gt lat_ranges[i]) and (LAT lt lat_ranges[i+1]) ) and ((LON gt lon_ranges[j]) and (LON lt lon_ranges[j+1]) ) )
            ai_in_coord = AI[ai_indices] ;, ai_in_coord

            Xcheck = xtrack[ai_indices]
            ;; Ignore missing data
            ai_in_range = ai_in_coord[where(abs(ai_in_coord) lt 20.0)]
            Xcheck = xcheck[where(abs(ai_in_coord) lt 20.0)]
            ;Xcheck = xtrack[ where(ai_in_coord eq ai_in_range)]
            ;Xcheck = xtrack[ where(AI eq ai_in_range and (xtrack eq 0))]
            ;Xcheck = xtrack[ where(AI eq ai_in_range and (xtrack eq 0 or (xtrack and 4) eq 4 ))]

            ; Only use AI values with an XTrackQC flag value of 0, meaning the
            ; values are not affected by the row anomaly.
            ; NOTE: Should XTrack values of 4 (meaning the values are affected,
            ;       but corrected) be used also?
            good_xcheck_indices = where(Xcheck eq 0)
            ai_in_range2 = ai_in_range[good_xcheck_indices]

            ;;if(n_elements(Xcheck) eq 1) then begin
            ;;    printf, 1, n_elements(ai_in_range2), n_elements(Xcheck), Xcheck[good_xcheck_indices]
            ;;endif else begin 
            ;;    printf, 1, n_elements(ai_in_range2), n_elements(Xcheck), Xcheck[0],' yay'
            ;;endelse

            ;; If there are any good data, add the average and standard devaition
            ;; of the data to the arrays
            ;print, lat_ranges(i),lon_ranges(j),ai_in_range
            if(n_elements(ai_in_range2) gt 1) then begin
                num_ai(i,j) = num_ai(i,j)+n_elements(ai_in_range2) 
                ; Keep a running total of average*num_obs which will be divided
                ; by the total number of obs after all data have been gathered
                avg_ai(i,j) = avg_ai(i,j)+total(ai_in_range2)
                ;avg_ai(i,j) = avg_ai(i,j)+(mean(ai_in_range)*num_ai(i,j))
;                std_ai(i,j) = stddev(ai_in_range) 
;                templatxlon = string(fix(lat_ranges(i)))+'x'+strtrim(string(fix(lon_ranges(j))),1)
;;                printf, 1, date, templatxlon, avg_ai(i,j), std_ai(i,j), num_ai(i,j)
            ; If there are no good data, add missing values to the arrays
            ;;;endif else begin
            ;;;    avg_ai(i,j) = -999.9 
            ;;;    std_ai(i,j) = -999.9 
            ;;;    num_ai(i,j) = 0
            endif ;endelse
        endfor
    endfor
    
endfor ; file loop
print, 'Finished with file loop'
print, 'Calculating averages'

; Calculate the total average for each lat/lon grid point
for i=0,num_lats-2 do begin
    for j=0,num_lons-2 do begin
        ; Calcualte the average AI values for the entire time period
        templatxlon = string(fix(lat_ranges(i)))+'x'+strtrim(string(fix(lon_ranges(j))),1)
        ;print, templatxlon, avg_ai(i,j), float(num_ai(i,j)), avg_ai(i,j)/float(num_ai(i,j))
        avg_ai(i,j) = avg_ai(i,j)/float(num_ai(i,j))

        ;if(num_ai(i,j) eq 0) then avg_ai(i,j) = -999.9
        ;    templatxlon = string(fix(lat_ranges(i)))+'x'+strtrim(string(fix(lon_ranges(j))),1)
        if(num_ai(i,j) ne 0) then printf, 1, date, templatxlon, avg_ai(i,j), num_ai(i,j)
        ;printf, 1, date, templatxlon, avg_ai(i,j), std_ai(i,j), num_ai(i,j)
    endfor    
endfor

print, 'Closing output file'
close, 1

END
