%% © 2012 Vanderbilt University %%

function PET = Hamon_PET(T_F,lat_list,count_loc,beg_row_mat,lat_mat)

% Hamon_PET.m is a function file that calculates the monthly potential 
% evapotranspiration (PET) using the Hamon method (cf. Hamon, W.R., 1963;
% International Association of Scientific Hydrology, Pub. 63). Latitude is 
% required as an input to calculate the average monthly solar insolation.
% The average  monthly solar insolation and temperature are then used to  
% calculate the monthly potential evapotranspiration. 
%
% NOTE:
% T_F is the temperature in degrees F. This data should be input as a
% column vector.
% lat_list is the list of latitudes in the input data. There should be one
% latitude for each location, and this informatio should be input as a
% column vector.
% count_loc is the number of different locations in the input data.
% beg_row_mat is a vector of the row numbers where the temperature data for
% each of the different locations begins when the data is listed in a  
% matrix such that years represent rows and columns represent months.
% lat_mat is the record of latitudes in degrees for all locations arranged
% in a column vector, where a latitude is associated with each YEAR of the
% record for all locations.

%% ASSIGN VARIABLES 

% temp is a matrix of the temperature data, for all locations, arranged 
% such that the rows represent the years in the record and columns 
% represent the months of the year.
temp = (reshape(T_F,12,length(T_F)/12))';

% days_mo is the number of days in each of the 12 months.
days_mo = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

% J is the Julian day of the middle of month i. Pre-allocate J for speed.
J = zeros(1,12);

for i = 1:12
    if i == 1
        J(i) = days_mo(i)/2;
    else
        J(i) = days_mo(i)/2 + sum(days_mo(1,1:i-1));
    end
end

PET = [];

for j = 1:count_loc % j is the counter for each of the different locations.
    %% CALCULATE THE NUMBER OF DAYLIGHT HOURS
    
    % lat_rad is a vector of the latitudes in radians for the different 
    % locations in the total data record.
    lat_rad = (lat_list(j)*2*pi)/360; 

    % TLAT is the negative tangent of the latitiude (where the latitude is 
    % in radians).
    TLAT = -tan(lat_rad);
    
    % Calculate PHI, the solar declination angle on Julian day (J).
    PHI = .4093.*sin(((2.*(pi())./365).*J)-1.405);
    % Calculate w (omega), the sunset hour angle on Julian day (J).
    w = acos(TLAT.*PHI);
    % Calculate D, the number of daylight hours (day length) on Julian  
    % day (J).
    D = 24.*w./pi;
    % Calculate DL, the day length normalizer.
    DL = D./12;
    
    b = beg_row_mat(j); % b is the row number of the temperature and
                        % precipitation data records (i.e., matrices) 
                        % where the data for location j begins.
    
    if j == count_loc
        e = length(lat_mat); % e is the row number of the temperature and
                             % precipitation data records (i.e., matrices)
                             % where the data for location j ends.
    else
        e = beg_row_mat(j + 1) - 1;
    end
    
    T = temp(b:e,1:12); % T is the temperature record, arranged as a 
                        % matrix, for location j.
    
    y = size(T,1); % y is the number of years in the temperature record for
                   % location j.
    %% CALCULATE PET
    
    % PET_mm is the potential evapotranspiration (in millimeters) for 
    % location j. Pre-allocate for speed.
    PET_mm = zeros((e - b + 1),12);
    
    for k = 1:y; % k is the counter for each year of the data on record 
                 % for each of the different locations.
        for i = 1:12; % i is the counter for each of the 12 months in a 
                      % year.
            if temp(k,i) <= 32
                % Haith and Shoemaker (cf. Haith and Shoemaker, 1987;  
                % Journal of the American Water Resources Association; 
                % Vol. 23, No. 3, June 1987) set PET equal to zero on days 
                % for which the mean monthly temperature is less than or  
                % equal to 32 degrees Fahrenheit.
                PET_mm(k,i) = 0;
            else
                PET_mm(k,i) = (2.1*(DL(i)^2)*0.6108*exp((17.27* ...
                5*(T(k,i) - 32)/9)/(237.3 + 5*(T(k,i) - 32)/9))/ ...
                (5*(T(k,i) - 32)/9 + 273.3))*days_mo(i);
            end
        end
    end
    
    % Convert PET (for all locations) from millimeters to inches.
    PET_in = (PET_mm./25.4);

    % Catalogue PET for all locations.
    PET = [PET; PET_in];
    
end

end