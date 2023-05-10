# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Short version notes
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

1  Latitude
2  Longitude
3  Original Value
4  Altered Value
5  Climo Value (difference)
6  Time (of day)
7  Szea
8  Vzea
9  Raza
10 Albedo 1 (used for binning)
11 Albedo 2 (not used for binning)
12 Normalized Radiance 1
13 Normalized Radiance 2
14 Cloud (flag?)
15 X-track flag (with removed row anomoly)
16 Ground flag
17 Row number (base 1)

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Short version notes
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

out_files-ltc3/     = same as ltc3, but with GPQF flag added at end.
                      Contains only rows that are always good 
                      (no row anomaly rows at all).
                      No files for April or May
    201708161504 contains 7644 lines
                 last two columns are GPQF flag and row number
out_files-ltc3/new/ = identical to old, but with fewer lines.
    201708161504 contains 7644 lines
                      Contains only rows that are always good 
                      (no row anomaly rows at all).
                 last column is the GPQF flag
out_files-ltc3/old/ = same as ltc3, but with no GPQF flag value.
                      longer file length than ltc3 or new. 
                      Contain data for April and May.
                      VALUES ARE NOT THE SAME AS LTC3, NEW, OR LTC4
    201708161504 contains 7945 lines
                 last column is the GPQF flag
out_files-ltc4      = same as ltc3/new, but with row number added at end.
                      201708161504 is IDENTICAL to ltc3 version

Looked at local versions of files, compared to raindrop
201807041812, housed in ltc3 locally
- ltc4: IDENTICAL
- ltc3

201507121137, stored locally in ltc3
- ltc3    : identical to local version in ltc3
- ltc3_old: contains more lines than local version, pert values are different
- ltc3_new: does not contain all rows, identical to ltc4
- ltc4:     identical to local version in ltc3

200804210805, not housed in ltc3 locally 
- ltc3    : NO DATA
- ltc3_old: contains all rows, much higher values than other types
- ltc3_new: does not contain all rows, identical to ltc4
- ltc4:     does not contain all rows, identical to ltc3_new

200607242155
- ltc3    : NO DATA
- ltc3_old: contains all rows, much higher values than other types
- ltc3_new: does not contain all rows, identical to ltc4
- ltc4:     does not contain all rows, identical to ltc3_new

# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
# Email communications
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

10:22 AM, 26 April 2021

Hi Blake,

I have not had much time as I wanted to have for this in the last week, but 
that is not to say that I had none.

I did not realize that you didn't have the data that was being processed. 
The process failed and I took that as completion until I checked while 
working on the line removal. I have finished the processing of the latter 
half (2012-2019) of the data we have, trained using the same months (April, 
May, June) from all years (with a +- 5 days for training periods).

My first attempt using python to train and remove the rows as we had talked
about did not go well and I am rethinking the approach slightly and may 
likely do it in c using the (now) processed data.

If you are going to create a set of plots for this data, I would appreciate 
read access to the directory.

Shawn

Forgot to mention the processed data are in:

raindrop.atmos.und.edu/Research/OMI/out_files-monthly/

-----------------------------------------

11:54 AM, 26 April 2021

Hi Shawn,

No worries! I see the data in that directory now, so thanks for working on 
that! 

My qual exam studying is taking over more and more of my time now, so I 
won't have much time to work on OMI stuff, but I will slowly chip away at a 
quick plotting script for these processed data files. When I have a script 
working and have some images generated, I will let you know where to find 
them!

That's unfortunate that your first attempt didn't go well. I looked at the 
OMI data in the HDF5 files and tried to see which index was causing the 
problem, and using both IDL and Python, it looks like row index 51 (row 52) 
is the problem row, which was surprising since I remember us saying that row
53 was the problem row earlier. In your first Python attempt, were you 
removing row 53 or 52?

Blake

-----------------------------------------

11:36 AM, 13 May 2021

Hi Blake,

In case you didn't already see it, the OMI processing has completed.

The earlier and later periods are trained on 'their' half of the data. 2012 
is a bit of an overlap. I have the logs of how they were done and when you 
need the specific information, I can get it. They still need to be plotted.

In general, each output month was trained on the same months as input (with 
a five day buffer added to training data) from the entire (half) period, 
including itself.

Data are in:
/Research/OMI/out_files-monthly/

Shawn

-----------------------------------------

10:15 AM, 1 July 2021

Hi Blake,

Been too long for this, but data is now being processed and output using the
cleaned input tracks. I am currently running as follows: For any given 
"reconstruction" month, the model is trained on the same months from the 
previous, same, and next year. The exception being 2020, which is only 
trained on the previous and same year. The longer training periods were 
causing biases in my opinion, as well as significant processing problems.

Just a quick question, do you start your numbering at 1? meaning tracks 1-60?
I had to subtract 1 from your values to remove the appropriate tracks.

The new data are on Raindrop in:

/Research/OMI/out_files-monthly_test/

As before, my initial output is in #####_CLI.txt. I am moving them over 
periodically as the processing finishes to the expected format.

Additionally, I am including a modified version of your excellent plotting 
code. The command line and basic output is the same, but I have thrown an 
informative histogram over the colorbar. Please feel free to re-use. The 
output is 'odata' for original, 'fdata' for filter, and 'pdata' for 
processed. Since the bad tracks are removed with a '-1' xtrack flag, they 
are not included in any of the figures.

Shawn

-----------------------------------------

3:43 PM, 8 July 2021

Hi Shawn,

I'm looking at your OMI files in /Research/OMI/out_files-monthly_test, and 
it looks like the earliest files here are from 2006. I remember there being 
2005 files in your previous dataset (in out_files-monthly.20210518), so do 
you know if you processed 2005 in your latest run?

Thanks!

Blake

-----------------------------------------

5:12 PM, 8 July 2021

Hi Blake,

I will double check. But, if memory serves. 2005 was used to train, but was
not used to reconstruct (since no 2004) wa available and I was using a +/-1
year window.

I am processing again, using all years to make climatology. Jianglong 
doesn't want a trend removal instroduced by the temporal window. The 
problem is computing the climo, I need to modify the code a bit to only 
compute the climo once. For this run, it will be 2006-2020 (15years). 
I will check to make sure 2005 is good, and if so, I will include it 
(16 years)

Shawn

-----------------------------------------

1:08 PM, 19 July 2021

Hi Blake,

I am running the OMIprocess with the single climo (for each month). It's 
taken a little longer than expected, but is about halfway through now, it 
was started on Saturday AM.

So far, April and May for all years are available. June should finish some 
time tonight.

The output files are in (please ignore the 'monthly' in the title, that's 
wrong);

/Research/OMI/out_files-monthly_test_ltc/

JZ mentioned that you had a new track removal file? I will need to run the 
climo again, but you are welcome to send that over when you have a chance.

Shawn

-----------------------------------------

2:36 PM, 29 July 2021

Hi Blake,

OMI process has finished with the new climo parameters. The data are 
available in:

/Research/OMI/out_files-ltc

and the figures in:

/Research/OMI/out_files-ltc/figs

I have also included the source code that generated all files with this 
email.

Shawn

-----------------------------------------

11:41 AM, 11 August 2021
Me to shawn:

Hi Shawn,

It took me a while, but I finally figured out what was causing the bad 
data for July and August 2019. You can see in the plot of the climatology 
and trends for each of the 3 datasets (original, JZ, and your data) that 
there's a similar ring of positive trend around about the North Pole in 
the uncorrected data and your perturbed data, which looks like row 
contamination. I re-ran my remote-ocean monthly climatology code a large 
handful of times, removing each sensor row one by one, and I finally 
determined that row 43 was the "problem child" for July and August 2019. 
You can see this in the drift check plot I also attached, in which the only 
difference between the two plots is that I excluded row 43 in July and August
2019 in the bottom plot, where the giant spike at the end is removed. 
I updated the bad row database to include row 43 as a bad row for all days 
in July and August of 2019, and the updated file is attached.

Could you process the OMI data again using this updated bad row file? I am
hopeful that this is the last time you should need to run your OMI data 
process because the climatology and trends in your data for August look 
great outside of this anomalous ring. 

Thanks!

Blake

-----------------------------------------

2:27 PM, 13 August 2021

Hi Blake,

New data are now available in:

/Research/OMI/out_files-ltc2/

and figures (some still being generated) in:

/Research/OMI/out_files-ltc2/figs/

I used a method to "patch" the data with the new rows. Could you double 
check those dates that you noticed a problem to make sure it is not still an 
issue?

Shawn


-----------------------------------------

9:54 AM, 16 August 2021

Hi Shawn,

I checked the output file from the climatology process that I just ran 
across the new data in /Research/OMI/out_files-ltc2/, and it is identical 
to the climatology I made with the previous dataset. To make sure that it 
wasn't just something wrong with my code, I compared the last two versions 
of some of your OMI files (/Research/OMI/out_files-ltc/ and 
/Research/OMI/out_files-ltc2/):
201905081110
201907040430
201907182111
201908011129
201908211738
and the old and new versions of these five files are also identical. Is it 
possible that the old bad row file was used in your last run?

Blake

-----------------------------------------

12:42 PM, 1 September 2021

Hey Blake,

OMI run is finished. New data are in: /Research/OMI/out_files-ltc3/

Do you want me to generate the figures?

Shawn

-----------------------------------------

4:36 PM, 1 September 2021

Hi Shawn,

I just finished one of my analyses on your most recent output, and it looks 
like there's a bad line getting through the analysis. I compared the most 
recent output with the previous version, and they appear to agree perfectly 
between 2005 and 2019. However, the data for 2020 have the same jumps that 
we saw in the data a couple of months ago that were caused by bad rows 
sneaking through the analysis. After comparing the data plot with the bad 
row plot I made with the 2020 row anomaly file (both attached), the jumps 
line up perfectly with the identified bad rows. What's weird is that row 53 
listed as the problem row in the 2020 row anomaly file, but row 53 shouldn't
have been included in your analysis because 53 is listed in the "good rows 
only" row anomaly file I sent you yesterday (or the other day... whenever 
that was). I am going to head out in a bit, but I will look at the raw data 
tomorrow and see if maybe another row is going bad at exactly the same times
as row 53 in 2020 but isn't getting picked up by my code.

I'll keep you updated on what I find, but for now, I would hold off on 
starting up your next process until we figure out what's causing the error in
this data.

Blake

-----------------------------------------

2:15 PM, 2 September 2021

Hi Shawn,

Just got back from my 3-hr marathon of classes, so I've got some time 
before Mike's retirement party at 3:30 if you still want to meet! 

I just plotted up some of the output files from your last run, and I can 
now confirm that row 53 is getting through the row removal process. I 
attached 3 figures: one for 20200416, one for 20200723, and one for 
20200816. The first two show the bad data in row 53 over the Arctic, but 
the last file is weird because row 53 is not shown. My guess is that this is 
because the actual XTrack flag might actually be set correctly for this date.
The other weird things is that row 8 is removed in the 20200723 figure, which
is when the row 8 data are actually missing in the HDF5 files, but it is not
removed for the other two figures. Neither row 8 nor row 53 should be 
making it through because they're both flagged as bad rows in the row file, 
so I'm not really sure what's going on.

I will also be driving (riding, I guess) for most of the day tomorrow, so 
if you aren't free to meet today, we can just plan to meet next Tuesday! 

Blake

-----------------------------------------

7:11 PM, 8 September 2021

They are now located in

/Research/OMI/out_files-ltc3/new
and
/Research/OMI/out_files-ltc3/old

Shawn


-----------------------------------------

12:34 PM, 9 September 2021

Hey Blake,

I just noticed that the output is missing the "row" column, but the other 
processing we talked about looks good. I apologize for my mistake. I am 
testing the change to output the row, do you need this urgently?

Shawn

-----------------------------------------

12:37 PM, 9 September 2021

Hi Shawn,

Awesome! I'm glad the other processing looks good! That's okay, I don't need
the data urgently. I have other analysis code that I'm working on, so I can 
keep working on that code while the data get reprocessed.

Thank you!

Blake
 
-----------------------------------------

6:20 PM, 14 September 2021

HI Blake,

Data is finished and ready in:

/Research/OMI/out_files-ltc4/

I reduced the width of the columns, otherwise the format is the same with the
addition of the row as the track (row) as the last column.

Shawn

-----------------------------------------

10:46 AM, 15 September 2021
 
Hi Shawn,

Great, thank you! 

By the way, I processed the new data, and it looks like it fixed the trends 
around the 82N latitude we saw in the previous dataset, so the fixed worked!

-----------------------------------------

12:19 PM, 26 November 2021

Hi Blake,

Here is the reference to it.  Turns out the June-Sept are in ltc3 and 
April-May in ltc4. My intent was to put them all in ltc4 but I only changed 
the settings for the first process of the three. My apologies again for the 
inconvenience. If you need the prior ltc3 runs, they are in old/

I have all the logs if you have specific questions.

Shawn

-----------------------------------------

1:11 PM, 29 November 2021

Hi Blake,

With regard to your more recent email regarding the two periods, I think 
there is some misunderstanding. The June-Sept are (were) in ltc3 and the 
April and May in ltc4. 

The 'over-written' files that you mentioned were created in the absence of 
the respective file. I suspect it has to do with your file open command, 
maybe the mode parameter.

In any event, I don't think there is a need to reprocess. I am simply going 
to copy the files from ltc3 to ltc4, where they should have been in the 
first place.

Sorry for my initial mistake in the output directories. I have removed the 
0 size files and am currently copying the April-May to ltc4

Shawn

-----------------------------------------

11:50 AM, 14 January 2023

Hi Blake,

The relevant code is the in the OMIprocess.py file

szea_bins=np.linspace(0,90,37)       # 37 bin edges >> 36 bins >>  
                                       90 / 36 =  2.5 degree bins
vzea_bins=np.linspace(0,90,37)       # 37 bin edges >> 36 bins >>  
                                       90 / 36 =  2.5 degree bins
raza_bins=np.linspace(-1,181,92)     # 92 bin edges >> 91 bins >> 
                                       182 / 91 =  2.0 degree bins
salb_bins=np.linspace(0.0,1.0,21)    # 21 bin edges >> 20 bins >> 
                                       1.0 / 20 = 0.05 albedo bins
flag_bins=np.linspace(0,7,8)          #  8 bin edges >>  7 bins

The only dimensions on which the climo is trained is the solar zenith, 
viewing zenith, and relative azimuth angles, with the solar albedo and 
ground type, as shown above. The output files include the original point 
data and the climo value that was removed from it.

The output fields for a particular orbit have the following columns:

Latitude
Longitude
Original Value
Altered Value
Climo Value (difference)
Time (of day)
Szea
Vzea
Raza
Albedo 1 (used for binning)
Albedo 2 (not used for binning)
Normalized Radiance 1
Normalized Radiance 2
Cloud (flag?)
X-track flag (with removed row anomoly)
Ground flag
Row number (base 1)

So, it should be possible to get the point data if they want the original 
location from the data. The climo value is not associated with a particular 
location. There is no time dependence or regression of historical values. 
The one caveat to that is that we train a given month only with the same 
months from other years (with an additional 5 days before and after the 
chosen month, eg. see: fast_extract1.sh, etc.). There is also no 5-d 
histogram, but one could probably be generated from the climo value if we 
slice it effectively for display. This has not been done yet.

The climo.bin file is just the bin_mean variable (with a shape of 
(len(szea_bins), len(vzea_bins), len(raza_bins), len(salb_bins), 
len(flag_bins))) that has been "pickled" using the standard python module. 
If you need help reading it let me know. The file is specific for each run, 
so would need to be regenereated.

I hope that helps.

Best regards,
Shawn
