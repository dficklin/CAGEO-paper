%% © 2014 Indiana University %%

function [PET] = Penman_PET(T_F_all,T_F_minmax,rsds_all,u_2_all,Humidity_data_all,Alb_Data,Elevation_Data,lat_list,count_loc,beg_row_col,Days_in_month)
%Penman_PET.m is a function that calculates potential evapotranspiration (PET) using the 
%Penman-Monteith method (Allen RG, Pereira LS, Raes D, Smith M, 1998, Crop 
%evapotranspiration-guidelines for computing crop water requiremenTa. FAO Irrigation and Drainage Paper 56)
%%%
%%%Input(s):
%%%       T_F_all: Mean air temperature in degrees Farenheit. This data should be input as a column vector.
%%%       T_F_minmax: Minimum and maximum air temperature in degrees Farenheit. This data should be input as a n by 2 matrix(i.e.,[min max]).
%%%       Note: If T_F_minmax is not available, then use [].
%%%       rsds_all: Average daily downward shortwave radiation at surface (W m-2). This data should be input as a column vector.
%%%       u_2_all: Wind speed at 2m height m/s. This data should be input as a column vector.
%%%       Humidity_data_all: Specific humidity (kg/kg). This data should be input as a column vector.
%%% Note: If "relative humidity (%)" is loaded instead of specific humdity, comment lines 82 and 92, and uncomment line 93.
%%%       Alb_Data (should be between [0-1]): Albedo. This data should be a constant, or a column vector.
%%%       Elevation_Data: Surface elevation in meters. This data should be a constant, or a column vector.
%%%       lat_list: Latitude vector. This data should be a column vector.
%%%       count_loc: Number of locations.
%%%       beg_row_col: Vector of the row numbers, which contains the column where the temperature data for 
%%%       each of the different locations begins.

%%%Output:
%%%       PET: Penman–Monteith potential evapotranspiration (inch/month)

Fahrenheit2Cel = @(T_F)(T_F-32)*5/9;%%%Farenheit to Celsius
if length(T_F_minmax)==0 %If mean temperature is used:
    Ta_all=Fahrenheit2Cel(T_F_all);%%%Air temperature in Celsius. 
    % Ta_all=T_F_all;%%%Uncomment this line if input temperature is in Celsius
else  %If min and max temperature are used:
    Ta_min_all=Fahrenheit2Cel(T_F_minmax(:,1));%Min Air temperature in Celsius
    Ta_max_all=Fahrenheit2Cel(T_F_minmax(:,2));%Max Air temperature in Celsius
    Ta_all = (Ta_min_all + Ta_max_all)/2;%%%Mean temperature in Celsius.
    %     Ta_min_all=T_F_minmax(:,1);%%%Uncomment this line if input temperature is in Celsius
    %     Ta_max_all=T_F_minmax(:,2);%%%Uncomment this line if input temperature is in Celsius
    %     Ta_all = (T_F_minmax(:,1) + T_F_minmax(:,2))/2;%%%Uncomment this line if input temperature is in Celsius
    
end
beg_row_col_modified= [beg_row_col;length(T_F_all)+1];
PET=[];
Net_all=[];
e_s_all=[];
e_a_all=[];
rh_all=[];
tmp_all=[];

Pressure_all=[];
for i=1:count_loc %Number of locations
    PET_mm=[];
    Ta=Ta_all(beg_row_col_modified(i):beg_row_col_modified(i+1)-1);
    if length(T_F_minmax)>0
        Ta_min=Ta_min_all(beg_row_col_modified(i):beg_row_col_modified(i+1)-1);
        Ta_max=Ta_max_all(beg_row_col_modified(i):beg_row_col_modified(i+1)-1);
    end
    
    rsds=rsds_all(beg_row_col_modified(i):beg_row_col_modified(i+1)-1);
    u_2=u_2_all(beg_row_col_modified(i):beg_row_col_modified(i+1)-1);
    u_2(u_2<=0.5)=0.5;%Limit very low wind speeds for the calculations (m/s)
    
    G=0;%Soil heat flux density (MJ m-2/day). It is assumed to be negligible, so G=0;
    if length(T_F_minmax)==0
        e_s=.611.*exp(17.3*Ta./(Ta+237.3));%%% Saturation vapor pressure in kPa
    else %If min/max temperature are used, then e_s is calculated for min and max remperature and then averaged, 
         %as per Allen et al. 1998, FAO Chapter 3, Eqn 12
        e_s_min=.611.*exp(17.3*Ta_min./(Ta_min+237.3));% Saturation vapor pressure in kPa for Ta_min; Allen et al. 1998, FAO Chapter 3, Eqn 11
        e_s_max=.611.*exp(17.3*Ta_max./(Ta_max+237.3));% Saturation vapor pressure in kPa for Ta_max; Allen et al. 1998, FAO Chapter 3, Eqn 11
        e_s=(e_s_min+e_s_max)/2;% Mean saturation vapor pressure in kPa; Allen et al. 1998, FAO Chapter 3, Eqn 12
    end
    
    if length(Elevation_Data)==1
        %If constant elevation is used:
        elevation=Elevation_Data*ones(length(rsds),1);
    else
        %If elevation file is loaded:
        elevation=Elevation_Data(i)*ones(length(rsds),1);
    end
    Pressure_sea_level=101.3; %Sea level Air pressure in kPa (=1013 mbar)
    Pressure = Pressure_sea_level*((293-0.0065*elevation)/293).^5.26; %Pressure in kPa; Allen et al. 1998, FAO Chapter 3, Eqn 7
    Psy_cnst = (0.665e-3) * Pressure;    %kPa/C; Allen et al. 1998, FAO Chapter 3, Eqn 8
    
    %%%If relative humidity is loaded instead of specific humdity,
    %%%comment lines 82 and 92, and uncomment line 93.
    Spec_Hum=Humidity_data_all(beg_row_col_modified(i):beg_row_col_modified(i+1)-1);
    %Relative Humidity in (%) estimated using the ideal gas law, and the
    %approximations of q~=mixing ratio and RH~=q/qs (specific humidity/saturated specific
    %humidity). American Meteorological Society glossary: http://glossary.ametsoc.org/wiki/Specific_humidity;
    %American Meteorological Society glossary: http://glossary.ametsoc.org/wiki/Relative_humidity;
    %Stull, R.B., 1995, Meteorology Today for Scientists and Engineers,West
    %Publishing Company, San Francisco, p. 86-87, Eqns 5.4 and 5.6.
    %Another relationship could employ a slightly different
    %approximation of q (specific humidity):
    %rh=(((Spec_Humid*(Pressure-0.378*e_s))/(e_s*0.622))*100
    rh=100*(Spec_Hum.*Pressure)./(e_s*0.622);  
%   rh=Humidity_data_all(beg_row_col_modified(i):beg_row_col_modified(i+1)-1);
    
    %Actual vapor pressure (kPa); %Allen et al. 1998, FAO Chapter 3, Eqn 17
    e_a=rh.*e_s/100 ;
    
    %Slope vapor pressure curve in kPa/C; Allen et al. 1998, FAO Chapter 3, Eqn 13
    Delta = 4098*((.6108.*exp(17.27*Ta./(Ta+237.3))))./(Ta+237.3).^2;  
    
    kdown = rsds * 0.0864;    %W/m2 to MJ/m2day
    % net radiation from downward shortwave radiation (rsds):
    if length(Alb_Data)==1
        %If constant Albedo is used:
        Albedo=Alb_Data*ones(length(rsds),1);
    else
        %If Albedo file is loaded:
        Albedo=Alb_Data(i)*ones(length(rsds),1);%%%If one Albedo is used for one point over time.
    end
    
    [ netrad ] = kdown2net_rad(Ta,kdown,e_a,Albedo,lat_list(i));
    Net_all=[Net_all;netrad];
    PET_mm=((.408.*Delta.*(netrad-G)+900.*Psy_cnst.*u_2.*(e_s-e_a)./(Ta+273))./(Delta+(Psy_cnst.*(1+0.34.*u_2))));
    
    %Finding negative temperature indices:
    PET_mm(Ta<0)=0;%For negative (mean) temperatures in Celsius: PET=0.
    PET_mm=reshape(PET_mm,12,length(PET_mm)/12);%Reshaping the PET vector into m by 12 matrix, where m is the number of years.
    PET_mm=PET_mm';
    PET_mm_mon=PET_mm.*Days_in_month;%Convert mm/day to mm/month.
    PET_in=PET_mm_mon./25.4;%Convert mm/month to inches/month.
    PET=[PET;PET_in];%PET in inches/month
    
    e_s_all=[e_s_all;e_s];
    e_a_all=[e_a_all;e_a];
    rh_all=[rh_all;rh];
    tmp_all=[tmp_all;Ta];
    Pressure_all=[Pressure_all;Pressure];
end



