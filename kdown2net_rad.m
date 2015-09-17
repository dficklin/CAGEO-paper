%% © 2014 Indiana University %%

function [ netrad ] = kdown2net_rad(TCa_mean,kdown,e_a,Albedo,lat)

%kdown2net_rad.m is a function that calculates net radiation from downward shortwave radiation. 
%%%Input:
%%%       TCa_mean is the air temperature in degrees Celsius. This data should be input as a column vector.
%%%       kdown is the average daily downward shortwave radiation at surface (MJ/m2day). This data should be input as a column vector.
%%%       e_a is actual vapor pressure (kPa). This data should be input as a column vector. 
%%%       albedo (should be between [0-1]). This data should be input as a column vector. 
%%%       lat is the latitude.

%%%Output:
%%%       netrad: net radiation (column vector) in MJ/m2day.


%%%Julian day of the middle of each month:
% days_mo is the number of days in each of the 12 months.
% days_mo = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
% % J is the Julian day of the middle of month i. 
% julian = zeros(1,12);
% for i = 1:12
%     if i == 1
%         julian(i) = days_mo(i)/2;
%     else
%         julian(i) = days_mo(i)/2 + sum(days_mo(1,1:i-1));
%     end
% end

%%%Most equations below are sourced from: 
%%%Allen RG, Pereira LS, Raes D and Smith M, 1998. Crop evapotranspiration: 
%%%Guidelines for computing crop requirements. Irrigation and Drainage Paper No. 56, 
%%%FAO, Rome, Italy.

%%%Pre-allocated julian day for speed:
julian=[15.5000   45.0000   74.5000  105.0000  135.5000  166.0000  196.5000  227.5000  258.0000  288.5000  319.0000  349.5000];
Num_years=length(kdown)/12;

julian_repmat=repmat(julian',Num_years,1);

% lat_rad is the latitudes in radians for a given location.
lat_rad = lat*2*pi/360; %Convert latitude in degrees to radians
    
% TLAT is the negative tangent of the latitude (where the latitude is
% in radians).
TLAT = -tan(lat_rad);

stbo = 4.903e-9;  %Stephan Boltzmann constant, MJ/m2K4day
IO = 0.082;     %Solar constant; MJ/m2min
TKa_mean = TCa_mean+273.15;%Celsius to Kelvin.

% Calculate PHI, the solar declination angle on Julian day (J) in radians:
PHI = .4093*sin(((2*pi/365)*julian_repmat)-1.39); %Allen et al. 1998, FAO Chapter 3, Eqn 24.

%sunset hour angle; radians, Allen et al. 1998, FAO Chapter 3, Eqn 25 
omega_sunset = acos((TLAT).*tan(PHI));  

%Allen et al. 1998, FAO Chapter 3, Eqn 28 (partial)
cosz=(omega_sunset.*(sin(lat_rad).*sin(PHI)))+(cos(lat_rad).*cos(PHI).*sin(omega_sunset)); 

%inverse relative distance Earth-Sun, radians, Allen et al. 1998, FAO Chapter 3, Eqn 23
d_r=1+(0.033*cos(2*pi*julian_repmat/365));    

%Extraterrestrial radiation or radiation at top of atmosphere
%MJ/m2day; Allen et al. 1998, FAO Chapter 3, Eqn 28 (partial)
atmosrad = ((24*60)/pi)*IO*d_r.*cosz;   

%Clear sky radiation approximation, using coefficients suggested by Allen et al. 1998
%MJ/m2day; Allen et al. 1998, FAO Chapter 3, Eqn 36
Rso = (0.25+0.5)*atmosrad;

%Solar radiation or net incoming shortwave radiation, MJ/m2day; Allen et al. 1998, FAO Chapter 3, Eqn 38
solrad = kdown.*(1-Albedo); 

%Solar radiation fraction; dimensionless, Allen et al. 1998, FAO Chapter 3, Eqn 39 (partial)
%rs [kdown]/rs0 [clear-sky radiation] in Eqn 39
sun = kdown./Rso;  
%add limit to sun of <=1, so if sun > 1, then 1:
sun(find(sun>1))=1;

%Terrestrial radiation, or net longwave radiation, MJ/m2day, Allen et al. 1998, FAO Chapter 3, Eqn 39
terrad = (stbo*TKa_mean.^4).*(0.34-(0.14*(e_a.^0.5))).*((1.35*(sun))-0.350);

%Net radiation = incoming net shortwave - outgoing net longwave; Allen et al. 1998, FAO Chapter 3, Eqn 40
netrad = solrad - terrad; 
   
end

