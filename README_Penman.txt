Supplemental Material for

Incorporation of the Penman-Monteith potential evapotranspiration method into a Palmer Drought Severity Index tool

Darren Ficklin, Sally Letsinger, Hamed Gholizadeh, and Justin Maxwell 

(Indiana University, Bloomington, Indiana)


Introduction
This .zip file contains new and amended MATLAB code to add the calculation of the Penman-Monteith potential evapotranspiration method to the Jacobi et al., 2013, PDSI tool, as well as example data files. Descriptions of the new and amended MATLAB scripts can be found below.  All necessary data files (AWC, Albedo, Elevation, and micrometerological variables) are provided for a time series of 35 years (1979-2013) at one latitude/longitude coordinate to allow users to execute a full run of the Penman-Monteith calculations in the PDSI tool. 

The original README.txt file is reproduced at the bottom of this file for easy reference for codes that were not modified in this update of the tool. 

*** = modified from the original Jacobi et al., 2013, code.

1. ***PDSI_Tool_Launcher.m - This file opens the graphical user interface (GUI) and the PDSI_Tool_Launcher.fig, thus launching the tool. This file has been modified from the original Jacobi et al., 2013, code.

2. Penman_PET.m - The Penman_PET function calculates the potential evapotranspiration using the Penman-Monteith method.

3. kdown2net_rad.m - A subroutine of Penman_PET.m that calculates net radiation from inputs of month, latitude, downwelling shortwave radiation (W/m2), albedo, and temperature.

4. Days_each_month.m - A function called by Penman_PET.m that returns the number of days in each month for a given year. 

5. *** Hamon_PET.m - Modified from the original Jacobi et al. (2013) to allow the number of simulation years to be greater than 12.

6. *** Z_Index.m - Modified from the original Jacobi et al. (2013) to allow the number of simulation years to be greater than 12.


Example data files from 49.44153595 -95.14763641 (just north of Lake of the Woods, Minnesota, USA) from 1979-2013: 

1. awc.txt - This file contains the Available Water Content (AWC) value for Lake of the Woods, Minnesota.

1.1 Column "AWC", AWC value for Lake of the Woods, Minnesota, inches

2. albedo.txt - This file contains the average annual albedo value for Lake of the Woods, Minnesota. 

2.1 Column "Albedo", Albedo value for Lake of the Woods, Minnesota, dimensionless, fraction. 

3. elevation.txt - This file contains the elevation value for Lake of the Woods, Minnesota. 

3.1 Column "Elevation", Elevation value for Lake of the Woods, Minnesota, meters

4. PDSI_PM_input_example.txt - This file contains climate data for Lake of the Woods, Minnesota, as well as the latitude and year for each observation.  

4.1 Column "Latitude", latitude of climate division centroid, degrees

4.2 Column "Year", year of the temperature and precipitation observation, years, data is monthly so each year value repeats 12 times

4.3 Column "Temperature_min", monthly average minimum temperature observation, degrees Fahrenheit (converted to C or K in the code)

4.4 Column "Temperature_max", monthly average maximum temperature observation, degrees Fahrenheit (converted to C or K in the code)

4.5 Column "Precipitation", total monthly precipitation, inches

4.6 Column "Solar radiation", monthly average downwelling shortwave radiation, W/m2

4.7 Column "Wind speed", monthly average wind speed velocity, m/s

4.8 Column "Specific humidity", monthly average specific humidity, kg/kg (converted to relative humidity in the code; see code comments for instructions on how to use relative humidity directly as an input.)

******The calibration time period for the above time period is 1950-1999. This can be changed to any time period by modifyingline 379 to the beginning year of calibration, line 383 to the ending year of calibration, line 397 to the beginning year of calibration and line 401 to the ending year of calibration. ***************  


=============================================
ORIGINAL README FILE FROM JACOBI ET AL., 2013
=============================================

Auxiliary Material for

A Tool For Calculating the Palmer Drought Indices

John Jacobi, Debra Perrone, Leslie Lyons Duncan, George Hornberger

(Vanderbilt Institute for Energy and Environment, Vanderbilt University, Nashville, Tennessee)

Water Resources Research, 2013

Introduction
This .zip file contains the MATLAB code for the PDSI tool, as well data files and a user manual. The user manual also contains detailed instructions on the use of the tool. Detailed descriptions of the different MATLAB scripts can be found be below, in the user manual, and in the comments of the MATLAB scripts themselves.  Two of the data files (GA_02_AWC and GA_02_TempAndPrecip) are provided as sample data files to help the user begin using the code. The NCDC_soilcnst data file is for those users who wish to replicate the PDSI calculations exactly as they appear in the NCDC code. 

1. BackTrack.m - This function backtracks through previous PX1 and PX2 values. Backtracking occurs in two instances: (1) after the PPe reaches 100 and (2) when the PPe is zero. In either case, the backtracking function works by backtracking through PX1 and PX2 until reaching a month where PPe = 0. Either PX1 or PX2 is assigned to X as the backtracking progresses.

2. Between0s.m - Established dry and wet spells are indicated by PPe = 0, and abatement of such spells are signified by 0 < PPe < 100. Frequently a possible abatement of an established spell is interrupted by a return to dry or wet conditions, without the spell having ended, and this is indicated by a return of PPe = 0. In such instances – i.e., when non-zero, non-one hundred PPe values occur between values of PPe = 0 – this function is called. 

3. Count_Loc.m - This is an indexing function that counts the number of locations for which the PDSI is being calculated and logs the beginning and end of each location’s data record.

4. Function_Ud.m - Called when there is an established wet spell, this function calculates the PPe value for a given month and then either checks for backtracking opportunities or continues through the code.

5. Function_Uw.m - Called during an established drought, this function calculates the PPe value for a given month and then either checks for backtracking opportunities or continues through the code.

6. GA_02_AWC.txt - This file contains the Available Water Content (AWC) value for the second Georgia climate division. A list of AWC values for U.S. climate divisions can be found in the NCDC_soilcnst.txt file. This file is provided as an example and can be used to test the tool upon download.

6.1 Column "AWC", AWC value for second Georgia climate division, inches

7. GA_02_TempAndPrecip.txt - This file contains temperature and precipitation data for the second Georgia climate division, as well as the latitude and year for each observation.  This file is provided as an example and can be used to test the tool upon download.

7.1 Column "Latitude", latitude of climate division centroid, degrees

7.2 Column "Year", year of the temperature and precipitation observation, years, data is monthly so each year value repeats 12 times

7.3 Colum "Temperature", monthly average temperature observation, degrees Fahrenheit

7.4 Column "Precipitation", monthly precipitation observation, inches

8. Hamon_PET.m - Hamon_PET.m is a function file that calculates the monthly potential evapotranspiration (PET) using the Hamon method. Latitude is required as an input to calculate the average monthly solar insolation. The average monthly solar insolation and temperature are then used to calculate the monthly potential evapotranspiration.

9. Main.m - This function calculates PX1 and PX2 and calls the backtracking loop. If the absolute value of PX1 or PX2 exceeds 1, that value becomes the new PX3.

10. NCDC_soilcnst.txt - This file contains the climate division information needed to accurately replicate the NCDC calculations of PDSI. The first column is the climate division number, the second column is the AWC value for that division, the third column is the B value (which is the NCDC equivalent of Thornthwaite’s PET exponent, a, the fourth column is the climate division heat index, and the fifth column is the negative tangent of the latitude at the centroid of the division.

10.1 Column "Climate division", the US climate division number, climate division numbers

10.2 Column "AWC", the AWC value for the climate division, inches

10.3 Column "B", NCDC B value, unitless, this value is the equivalent of Thornthwaite's PET exponent, a

10.4 Column "Heat Index", NCDC heat index value, unitless

10.5 Column "Negative tangent of latitude", the negative tangent of the latitude at the centroid of the division, unitless

11. PDSI_Central.m - This function is the central component of the PDSI calculations. All sub-functions are called from PDSI_Central.m, and once calculations are complete, the PDSI and PHDI values are finalized and tabulated in this function.

12. PDSI_Tool_Launcher.fig - This file contains the supporting image for the GUI.

13. PDSI_Tool_Launcher.m - Running this file opens the graphical user interface (GUI), thus launching the tool.

14. Thornthwaite_PET.m - The Thornthwaite_PET function calculates the potential evapotranspiration using Thornthwaite's method.

15. User_Manual.txt - An abridged version of the user manual for users without the software required to open PDFs.

16. User_Manual.pdf - A full user manual for the tool. Directions on tool operation and data setup and arrangement can be found in this document.

17. WaterBalance.m - This function calculates the Thornthwaite water balance using inputs from the PET function and user-loaded precipitation data.

18. Z_Index.m - This function calculates Palmer's Z index using inputs from the water balance
function.
