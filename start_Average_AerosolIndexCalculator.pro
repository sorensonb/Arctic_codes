; 
; NAME:
;   start_Average_AerosolIndexCalculator.pro
; 
; PURPOSE:
;   Automate the calculation of monthly-averaged aerosol index values. 
; 
; CALLS:
;   Average_AerosolIndexCalculator.pro
; 
; MODIFICATIONS:
;   Blake Sorenson <blake.sorenson@und.edu>     - 2018/11/2
;     Updated comments (written earlier this semester)
;



Pro Run_AI, first

; Compile Average_AerosolIndexCalculator

start_dates = [200410,200501,200601,200701,200801,200901,201001,201101,201201,201301,201402,201501,201601,201701,201801,201901]
end_dates   = [200412,200512,200612,200712,200812,200912,201012,201112,201212,201312,201412,201512,201612,201712,201812,201902]

start_length = n_elements(start_dates)
end_length   = n_elements(end_dates)

;sdate = start_date
;edate = end_date
first_num=1


if(first eq 1) then first_num=1

; Since data for 2004 has already been processed, start the loop at 1 to begin
; processing data at 2005
for i=0,start_length-1 do begin
    sdate=start_dates[i]
    edate=end_dates[i]
    for j=sdate,edate do begin
        Average_AerosolIndexCalculator,j,first_num
        if(first_num eq 1) then first_num=0
    endfor
endfor

end
