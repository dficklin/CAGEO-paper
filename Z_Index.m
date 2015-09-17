%% © 2012 Vanderbilt University %%

function [Z_all] = ...
          Z_Index(P,PET,ET,PR,R,RO,PRO,L,PL,beg_row_mat,end_row_mat, ...
          count_loc,beg_calyr_mat,end_calyr_mat,beg_row_col,end_row_col)

% This function calculates Palmer's Z index using inputs from the water 
% balance function.

% NOTE: P, PET, ET, PR, R, RO, S, L, and PL SHOULD BE READ IN AS COLUMNS 
% AND SHOULD BE GIVEN IN INCHES.
% beg_row_mat is a vector of the row numbers where the temperature data for
% each of the different locations begins when the data is listed in a  
% matrix such that years represent rows and columns represent months.
% lat_col is the record of latitudes in degrees for all locations, 
% arranged in a column vector, where a latitude is associated with each
% observation (i.e., month) of the record.
% lat_mat_rec is the record of latitudes in degrees for all locations
% arranged in a column vector, where a latitude is associated with each
% YEAR of the record.

Z_all = [];

for j = 1:count_loc
    
    b = beg_row_mat(j); % b is the row number of the reshaped temperature
                        % and precipitation data records (i.e., matrices) 
                        % where the data for location j begins.
                        
    e = end_row_mat(j); % e is the row number of the reshaped temperature
                        % and precipitation data records (i.e., matrices)
                        % where the data for location j ends.                   
    
    g = beg_row_col(j); % g is the row number of the temperature and
                        % precipitation data records (i.e., column vectors)
                        % where the data for location j begins.
                        
    d = end_row_col(j); % d is the row number of the temperature and 
                        % precipitation data records (i.e, column vectors)
                        % where the data for location j ends.
    
    y_col = (d - g + 1)/12;  % y_col is the number of years in the data 
                             % record for location j arranged in a column.
                             
    % The potential values - PET, ET, PR, PL - and the actual values - R, 
    % RO, S, L, P - are reshaped for location j such that the rows of the 
    % matrix represent years and the columns of the matrix represent 
    % months.
    PET_Reshape = PET(b:e,1:12);
    ET_Reshape = (reshape(ET(g:d),12,y_col))';
    PR_Reshape = (reshape(PR(g:d),12,y_col))';
    PL_Reshape = (reshape(PL(g:d),12,y_col))';
    R_Reshape = (reshape(R(g:d),12,y_col))';
    RO_Reshape = (reshape(RO(g:d),12,y_col))';
    PRO_Reshape = (reshape(PRO(g:d),12,y_col))';
    L_Reshape = (reshape(L(g:d),12,y_col))';
    P_Reshape = (reshape(P(g:d),12,y_col))';
        
    % ALPHA, BETA, GAMMA, DELTA CALCULATIONS
    % A calibration period is used to calculate alpha, beta, gamma, and 
    % and delta, four coeffients dependent on the climate of the area being
    % examined. The NCDC and CPC use the calibration period January 1931
    % through December 1990 (cf. Karl, 1986; Journal of Climate and Applied 
    % Meteorology, Vol. 25, No. 1, January 1986).
    cy = end_calyr_mat(j) - beg_calyr_mat(j) + 1; % cy is the number of 
                                                  % years in the 
                                                  % calibration period for 
                                                  % location j.
    for i = 1:12
        % The coefficients are caclulated using average values for each 
        % month. The averages are calculated using the calibration period.
        P_bar(i) = sum(P_Reshape(beg_calyr_mat(j):end_calyr_mat(j), ...
                   i),1)/cy;
        ET_bar(i) = sum(ET_Reshape(beg_calyr_mat(j):end_calyr_mat(j), ...
                    i),1)/cy;
        PET_bar(i) = sum(PET_Reshape(beg_calyr_mat(j):end_calyr_mat(j), ...
                     i),1)/cy;
        R_bar(i) = sum(R_Reshape(beg_calyr_mat(j):end_calyr_mat(j), ...
                   i),1)/cy;
        PR_bar(i) = sum(PR_Reshape(beg_calyr_mat(j):end_calyr_mat(j), ...
                    i),1)/cy;
        L_bar(i) = sum(L_Reshape(beg_calyr_mat(j):end_calyr_mat(j), ...
                   i),1)/cy;
        PL_bar(i) = sum(PL_Reshape(beg_calyr_mat(j):end_calyr_mat(j), ...
                    i),1)/cy;
        RO_bar(i) = sum(RO_Reshape(beg_calyr_mat(j):end_calyr_mat(j), ...
                    i),1)/cy;
        PRO_bar(i) = sum(PRO_Reshape(beg_calyr_mat(j):end_calyr_mat(j), ...
                   i),1)/cy;
        % Calculate alpha.
        if PET_bar(i) == 0
            if ET_bar(i) == 0
                alpha(i) = 1;
            else
                alpha(i) = 0;
                fprintf('CHECK DATA: PET is less than ET.')
            end
        else
            alpha(i) = ET_bar(i)/PET_bar(i);
        end
        % Calculate beta.
        if PR_bar(i) == 0
            if R_bar(i) == 0
                beta(i) = 1;
            else
                beta(i) = 0;
                fprintf('CHECK DATA: PR is less than R.')
            end
        else
            beta(i) = R_bar(i)/PR_bar(i);
        end
        % Calculate gamma.
        if PRO_bar(i) == 0
            if RO_bar(i) == 0
                gamma(i) = 1;
            else
                gamma(i) = 0;
                fprintf('CHECK DATA: PRO is less than RO.')
            end
        else
            gamma(i) = RO_bar(i)/PRO_bar(i);
        end
        % Calcuate delta.
        if PL_bar(i) == 0
            if L_bar(i) == 0
                delta(i) = 1;
            else
                delta(i) = 0;
                fprintf('CHECK DATA: PL is less than L.')
            end
        else
            delta(i) = L_bar(i)/PL_bar(i);
        end
        
    end
    
    % CALIBRATED CAFEC, K, AND d CALCULATION
    % NOTE: 
    % The Z index is calculated with a calibrated K (weighting factor) but
    % a full record d (difference between actual precipitation and CAFEC -
    % climatically appropriate for exisiting conditions - precipitation).
    % CAFEC precipitation is calculated analagously to a simple water
    % balance, where precipitation is equal to evaporation plus runoff 
    % (and groundwater recharge) plus or minus any change in soil moisture 
    % storage. 
    for k = 1:cy
        if k == 1
            v(k) = beg_calyr_mat(j); % v is the row number of the reshaped
                                     % PET, PR, PRO, PL, and P matrices
                                     % that corresponds to the calibration
                                     % year k.
        else
            v(k) = v(k-1) + 1;
        end
        for i = 1:12
            % CAFEC_hat is calculated for month i of year k of the 
            % calibration period.
            CAFEC_hat(k,i) = alpha(i)*PET_Reshape(v(k),i) + beta(i)* ...
                             PR_Reshape(v(k),i) + gamma(i)* ... 
                             PRO_Reshape(v(k),i) - delta(i)* ...
                             PL_Reshape(v(k),i);
            % Calculate d_hat, the difference between actual precipitation
            % and CAFEC precipitation for month i of year k of the 
            % calibration period.
            d_hat(k,i) = P_Reshape(v(k),i) - CAFEC_hat(k,i);
        end
    end
    
    for i = 1:12
        % NOTE: D_hat, T_hat, K_hat, and z_hat are all calibrated
        % variables - i.e., they are calculated only for the calibration
        % period.
        % Calculate D_hat, the average of the absolute values of d_hat for
        % month i.
        D_hat(i) = mean(abs(d_hat(:,i)));
        % Calculate T_hat, a measure of the ratio of "moisture demand" to
        % "moisture supply" for month i and location j.
        T_hat(i) = (PET_bar(i) + R_bar(i) + RO_bar(i)) / ...
                   (P_bar(i) + L_bar(i));
        % Calculate K_hat, the denominator of the K equation for month i.
        K_hat(i) = 1.5 * log10((T_hat(i)+2.8)/D_hat(i)) + .50;
        % Calculate z_hat, the numerator of the K equation for month i.
        z_hat_m(i) = (D_hat(i)*K_hat(i));
    end
    
    z_hat = sum(z_hat_m);
    
    for i = 1:12
        % Calculate the weighting factor, K, using the calibrated 
        % variables K_hat and z_hat. The purpose of the weighting factors
        % is to adjust the  departures from normal precipitation d 
        % (calculated below), such that they are comparable among different 
        % locations and for different months. The K tends to be large in 
        % arid regions and small in humid regions (cf. Alley, 1984; Journal 
        % of Climate and Applied Meteorology, Vol. 23, No. 7, July 1984).
        K(i) = ((17.67)*K_hat(i))/(z_hat);
    end
    
    % FULL RECORD CAFEC AND d CALCULATION
%     for n = 1:length(P_Reshape)
       for n = 1:size(P_Reshape,1)
        for i = 1:12
            % Calculate the CAFEC precipitation for month i and year n of
            % the full record.
            CAFEC(n,i) = alpha(i)*PET_Reshape(n,i) + beta(i)*...
                         PR_Reshape(n,i) + gamma(i)*PRO_Reshape(n,i) - ...
                         delta(i)*PL_Reshape(n,i);
            % Calculate d_hat, difference between actual precipitation and
            % CAFEC precipitatio for month i and year n of the full record.
            d(n,i) = P_Reshape(n,i) - CAFEC(n,i);
            % Calculate the Z index or the "moisture anomaly index" for 
            % month i and year n of the full record.
            z(n,i) = K(i) * d(n,i);
        end
       end
    z_reshape = reshape((z'),size(z,1)*12,1);
    Z_all = [Z_all; z_reshape];
end

end