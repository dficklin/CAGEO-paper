%% © 2012 Vanderbilt University %%

function varargout = PDSI_Tool_Launcher(varargin)
% PDSI_Tool_Launcher MATLAB code for PDSI_Tool_Launcher.fig
%      PDSI_Tool_Launcher, by itself, creates a new PDSI_Tool_Launcher or raises the existing
%      singleton*.
%
%      H = PDSI_Tool_Launcher returns the handle to a new PDSI_Tool_Launcher or the handle to
%      the existing singleton*.
%
%      PDSI_Tool_Launcher('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PDSI_Tool_Launcher.M with the given input 
%      arguments.
%
%      PDSI_Tool_Launcher('Property','Value',...) creates a new PDSI_Tool_Launcher or raises 
%      the existing singleton*.  Starting from the left, property value 
%      pairs are applied to the GUI before PDSI_Tool_Launcher_OpeningFcn gets called. 
%      An unrecognized property name or invalid value makes property 
%      application stop. All inputs are passed to PDSI_Tool_Launcher_OpeningFcn via 
%      varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PDSI_Tool_Launcher

% Last Modified by GUIDE 13-Jun-2014 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PDSI_Tool_Launcher_OpeningFcn, ...
                   'gui_OutputFcn',  @PDSI_Tool_Launcher_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PDSI_Tool_Launcher is made visible.
function PDSI_Tool_Launcher_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PDSI_Tool_Launcher (see VARARGIN)

% Choose default command line output for PDSI_Tool_Launcher
handles.output = hObject;

set(handles.pet_buttongroup,'SelectionChangeFcn', ...
    @pet_buttongroup_SelectionChangeFcn);
set(handles.calyrs_buttongroup, 'SelectionChangeFcn', ...
    @calyrs_buttongroup_SelectionChangeFcn);
set(handles.radiobutton_thorn,'Value',1)
set(handles.radiobutton_noaa,'Value',1)
set(handles.help_text,'Visible','off')


set(handles.edit_const_alb,'string','0.23')
set(handles.edit_alb_file,'string','')
set(handles.edit_alb_file,'Enable','off')
set(handles.radiobutton1_const_Alb,'Value',1)
set(handles.pushbutton_browse,'Enable','off')
set(handles.radiobutton_alb_file,'Value',0)
set(handles.uipanel_albedo,'SelectionChangeFcn',@uipanel_albedo_SelectionChangeFcn);
set(handles.pushbutton_Elevation,'Enable','off')
set(handles.pushbutton_OK,'Enable','off')

set(handles.radiobutton1_const_Alb,'Enable','off')
set(handles.radiobutton_alb_file,'Enable','off')
set(handles.edit_const_alb,'Enable','off')





W1 = {'Welcome to the MATLAB Palmer Drought Indices Calculator.'};
W2 = {'Please load your data, select your desired PET calculation method'};
W3 = {'and calibration period, and run the program.'};
W4 = {'Click the "Help" button for directions and answers to any questions.'};
dash = {'--------'};
space = {''};
set(handles.welcome_text,'String',[space;W1;dash;W2;W3;dash;W4;]);
H4 = {'------- Climate Inputs -------'};
H5 = {'- The US Historical Climate Network (USHCN) has historical temperature and precipitation data available for US stations at http://cdiac.ornl.gov/epubs/ndp/ushcn/access.html.'};
H6 = {'- Data must be consecutive (i.e., no gaps) and combined into one .txt file. Do not include column headers.'};
H7 = {'- Data should be chronologically organized into eight columns: Column 1 is the latitude in degrees, Column 2 is the year, Columns 3 and 4 are min/max temperature (Farenheit),'};
H8 = {'  Column 5 is precipitation (inches), Column 6 is the average daily downward shortwave radiation at surface (W m-2), Column 7 is wind speed (m/s), and Column 8 is specific humdity (kg/kg).'};
H9 = {'- Note: Column 8 can be relative humidity (See the comments in the penman.m code).'};  
H10 = {'- Note: Hamon and Thornswaite methods can be used if Columns 6 to 8 are note available.'};
H11 = {'- Note: The code can be used if mean temperature is available instead of min/max temperature.'};
H12 = {'- Note that all observations must have a latitude and a year associated with it.'};
H13 = {'- If there are data for multiple locations, make sure that the data are arranged in the same order.'}; 
H14 = {'-------- Available Water Capacity (AWC) Input --------'};
H15 = {'- Different locations have different field capacities. An AWC value for each location must be loaded into the program in inches.'};
H16 = {'- AWC data should be organized in a column in one .txt file. The AWCs should be organized in the same location order as the temperature and precipitation data.'};
H17 = {'- Note that only one AWC value is needed per location (i.e. the number of AWC values must equal the number of stations).'};
H18 = {'-------- Albedo --------'};
H19 = {'- Different locations have different Albedo (%). Albedo is used in Penman method to calculate net radiation. If Albedo is not available, user can input a constant Albedo.'};
H20 = {'- Note that if Albedo file is loaded only one Albedo value is needed per location (i.e. the number of Albedo values must equal the number of stations).'};
H21 = {'-------- Elevation --------'};
H22 = {'- Different locations have different Elevation (m). If Elevation is not available, then sea level is used.'};
H23 = {'- Note that if Elevation file is loaded only one Elevation value is needed per location (i.e. the number of Elevation values must equal the number of stations).'};
H24 = {'------- Calibration Period -------'};
H25 = {'- The calibration period is used to calculate the "Climatologically Appropriate for Existing Conditions" (CAFEC) precipitation for a location.'};
H26 = {'- The CAFEC is used to calculate weighting factors used in the Z-Index calculation (see Palmer, 1965).'};
H27 = {'- NOAA uses the period January 1931 to December 1990 as its calibration period (see Karl, 1986), and this option is provided.'};
H28 = {'- In the absence of long data records or to use a more comprehensive timespan, the option is also provided to use the full record as the calibration period.'};
H29 = {'------- Operation -------'};
H30 = {'- To operate the GUI, load your input data using the "Load Data" buttons, select your desired PET calculation method and calibration period, and click "Run."'};
H31 = {'------- Output Details -------'};
H32 = {'- Units for the output of the PET calculations are in inches. Main drought index outputs are PET, Z-Index, PDSI, and PHDI (see Palmer, 1965).'};
H33 = {'- Results are output into a text file (Palmer.txt) that is saved in the working folder of the current directory. The text file should be opened in Notepad or imported into Excel.'};
H34 = {''};
H35 = {'* For additional help please refer to README_Penman.txt.'};
H36 = {'TO EXIT THIS HELP MENU, PRESS THE "HELP" BUTTON IN THE LOWER RIGHT CORNER'};
set(handles.help_text,'String',[H4;H5;H6;H7;H8;H9;H10;H11; ...
    H12;H13;H14;H15;H16;H17;H18;H19;H20;H21;H22; ...
    H23;H24;H25;H26;H27;H28;H29;H30;H31;H32;H33;H34;H35;H36]);

global RB_Cal RB_PET alb file_load Elevation_file_load 

% Set button values to defaults in case user doesn't press them.
RB_PET = 2;
RB_Cal = 1;
alb=1;
file_load=0;
Elevation_file_load=0;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PDSI_Tool_Launcher wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PDSI_Tool_Launcher_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_data.
function pushbutton_data_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global DATA lat_col yr_col T_F T_F_minmax P rsds u_2 Spec_Hum mean_temp

D = uiimport('-file');
fprintf('%s\n%s %s.\n','Climate data successfully loaded.');
E = fieldnames(D);
DATA = getfield(D,char(E));

% Extract latitude, temperature, and precipitation data from large data
% file.
lat_col = DATA(1:end,1);
yr_col = DATA(1:end,2);
if   (isempty(mean_temp)==1) | (mean_temp==1)
    T_F = DATA(1:end,3);
    T_F_minmax=[];
    P = DATA(1:end,4);
    if size(DATA,2)==7
        set(handles.radiobutton_penman,'Enable','on')
        rsds= DATA(1:end,5);
        u_2= DATA(1:end,6);
        Spec_Hum= DATA(1:end,7);%%%note: The Specific humidity unit is kg/kg here, not g/kg
    else
        disp('Note: Not enough input data to run Penman! You can only run Thornthwaite or Hamon.')
        set(handles.radiobutton_penman,'Enable','off')
        rsds=[];
        u_2=[];
        Spec_Hum=[];
    end
else
    T_F_minmax = DATA(1:end,3:4);
    T_F=(T_F_minmax(:,1)+T_F_minmax(:,2))/2;
    P = DATA(1:end,5);
    if size(DATA,2)==8
        set(handles.radiobutton_penman,'Enable','on')
        rsds= DATA(1:end,6);
        u_2= DATA(1:end,7);
        Spec_Hum= DATA(1:end,8);%%%note: The Specific humidity unit is kg/kg, not g/kg
    else
        disp('Note: Not enough input data to run Penman! You can only run Thornthwaite or Hamon!')
        set(handles.radiobutton_penman,'Enable','off')
        rsds=[];
        u_2=[];
        Spec_Hum=[];
    end
end
% --- Executes on button press in pushbutton_awc.
function pushbutton_awc_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_awc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global AWC

A = uiimport('-file');
fprintf('%s\n%s %s.\n','AWC data successfully loaded.');
W = fieldnames(A);
AWC_DATA = getfield(A,char(W));

% Extract AWC data from larger data file.
AWC = AWC_DATA(1:end,1);

% --- Executes on button press in togglebutton_help.
function togglebutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_help

button_state = get(hObject,'Value');
if button_state == get(hObject,'Max')
	set(handles.help_text,'Visible','on')
elseif button_state == get(hObject,'Min')
	set(handles.help_text,'Visible','off')
end

function pet_buttongroup_SelectionChangeFcn(hObject, eventdata)

global RB_PET alb

%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject);

switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'radiobutton_hamon'
        %execute this code when hamon_radiobutton is selected
        RB_PET = 1;
        set(handles.radiobutton1_const_Alb,'Enable','off')
        set(handles.radiobutton_alb_file,'Enable','off')
        set(handles.edit_const_alb,'Enable','off')
        set(handles.pushbutton_browse,'Enable','off')
        set(handles. pushbutton_Elevation,'Enable','off')
        set(handles.pushbutton_OK,'Enable','off')
        set(handles.pushbutton_run,'Enable','on')

    case 'radiobutton_thorn'
        %execute this code when thorn_radiobutton is selected
        RB_PET = 2;
        set(handles.radiobutton1_const_Alb,'Enable','off')
        set(handles.radiobutton_alb_file,'Enable','off')
        set(handles.edit_const_alb,'Enable','off')
        set(handles.pushbutton_browse,'Enable','off')
        set(handles. pushbutton_Elevation,'Enable','off')
        set(handles.pushbutton_OK,'Enable','off')
        set(handles.pushbutton_run,'Enable','on')


    case 'radiobutton_penman'
        %execute this code when penman_radiobutton is selected
        RB_PET = 3;
        %         set(handles.pushbutton_penman,'Enable','on')
        set(handles.radiobutton1_const_Alb,'Enable','on')
        set(handles.radiobutton_alb_file,'Enable','on')
        set(handles.edit_const_alb,'Enable','on')
        set(handles.pushbutton_OK,'Enable','on')
        set(handles.pushbutton_Elevation,'Enable','on')
        if alb==1
            set(handles.pushbutton_browse,'Enable','off')
        else
            set(handles.pushbutton_browse,'Enable','on')
            set(handles.edit_const_alb,'Enable','off')
        end
        set(handles.pushbutton_run,'Enable','off')
        
    otherwise
        % Code for when there is no match.
        
end
%updates the handles structure
guidata(hObject, handles);

function calyrs_buttongroup_SelectionChangeFcn(hObject, eventdata)

global RB_Cal

%retrieve GUI data, i.e. the handles structure
handles = guidata(hObject); 
 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'radiobutton_noaa'
      %execute this code when hamon_radiobutton is selected
      RB_Cal = 1;
    case 'radiobutton_full'
      %execute this code when thorn_radiobutton is selected
      RB_Cal = 2;
    otherwise
       % Code for when there is no match.
 
end
%updates the handles structure
guidata(hObject, handles);

% --- Executes on button press in pushbutton_run.
function pushbutton_run_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global RB_Cal RB_PET AWC lat_col yr_col T_F T_F_minmax P rsds u_2 Spec_Hum Alb_Data Elevation_file_load Elevation_Data


[beg_row_col,end_row_col,beg_row_mat,end_row_mat,count_loc, ...
           lat_list,lat_mat] = Count_Loc(lat_col);
if RB_Cal == 1
    % Use the calibration period that the NCDC and CPC use - January 1931
    % through December 1990 (cf. Karl, 1986; Journal of Climate and Applied
    % Meteorology, Vol. 25, No. 1, January 1986).
    
    % NOTE:
    % THE ORDERING OF LOCATIONS MUST BE CONSISTENT ACROSS INPUTS!
    % yr_col is the record of years for all locations, arranged in a column
    % vector, where a year is associated with each observation (i.e., 
    % month) of the record for all locations.
    % beg_row_col is a vector of the row numbers where the temperature data
    % for each of the different locations begins.
    % beg_row_mat is a vector of the row numbers where the temperature data 
    % for each of the different locations begins when the data is listed in 
    % a matrix such that years represent rows and columns represent months.
    % count_loc is the number of different locations in the total data 
    % record.
    
    % yr_mat_rec is the record of years for all locations, arranged in a
    % matrix, where a year is associated with each MONTH of the record for 
    % all locations.
    yr_mat_rec = (reshape(yr_col,12,length(yr_col)/12))';
    
    % yr is the record of years for all locations, arranged in a column 
    % vector, where each year of the record is listed only once (i.e., one 
    % year for each year of observations) for all locations.
    yr = yr_mat_rec(:,1);
    
    % count1931 is a counter that tracks the rows in the total data record
    % (i.e. column vector) that holds the January 1931 temperature and
    % precipitation data for each of the different locations.
    count1931 = 1;
    
    % count1990 is a counter that tracks the rows in the total data record
    % (i.e., column vector) that holds the December 1990 temperature and
    % precipitation data for each of the different locations.
    count1990 = 1;
    
    % beg_yr_col is a vector of the row numbers that hold the January 1931
    % temperature and precipitation data for each of the different 
    % locations.
    beg_yr_col = [];
    
    % end_yr_col is a vector of the row numbers that the December 1990
    % temperature and precipitation data for each of the different 
    % locations.
    end_yr_col = [];
    
    for j = 1:count_loc
        for m = 1:(end_row_col(j) - beg_row_col(j) + 1)
            if m == 1
                v(m) = beg_row_col(j);
                if yr_col(v(m)) == 1950
                    count1931(m) = count1931(m);
                    count1990(m) = count1990(m);
                    beg_yr_col = [beg_yr_col; count1931(m)];
                elseif yr_col(v(m)) == 1999
                    count1931(m) = count1931(m);
                    count1990(m) = count1990(m);
                    end_yr_col = [end_yr_col; count1990(m)];
                else
                    count1931(m) = count1931(m);
                    count1990(m) = count1990(m);
                end
                continue
            elseif m == (end_row_col(j) - beg_row_col(j) + 1)
                v(m) = end_row_col(j);
            else
                v(m) = v(m-1) + 1;
            end
            if yr_col(v(m)) == 1950
                count1931(m) = count1931(m - 1) + 1;
                count1990(m) = count1990(m - 1) + 1;
                beg_yr_col = [beg_yr_col; count1931(m)];
            elseif yr_col(v(m)) == 1999
                count1931(m) = count1931(m - 1) + 1;
                count1990(m) = count1990(m - 1) + 1;
                end_yr_col = [end_yr_col; count1990(m)];
            else
                count1931(m) = count1931(m - 1) + 1;
                count1990(m) = count1990(m - 1) + 1;
            end
        end
        beg_calyr_col(j) = min(beg_yr_col);
        beg_yr_col = [];
        beg_calyr_mat(j) = (beg_calyr_col(j) - 1)/12 + 1;
        end_calyr_col(j) = max(end_yr_col);
        end_yr_col = [];
        end_calyr_mat(j) = end_calyr_col(j)/12;
    end
    
else
    % Use the entire period of record as the calibration period (cf. Karl, 
    % 1986; Journal of Climate and Applied Meteorology, Vol. 25, No. 1, 
    % January 1986).
    
    for j = 1:count_loc
        beg_calyr_col(j) = 1; 
        end_calyr_col(j) = (end_row_col(j) - beg_row_col(j) + 1);
        beg_calyr_mat(j) = (beg_calyr_col(j) - 1)/12 + 1;
        end_calyr_mat(j) = end_calyr_col(j)/12;
    end
        
end

if RB_PET == 1
    PET = Hamon_PET(T_F,lat_list,count_loc,beg_row_mat,lat_mat);
elseif RB_PET == 2
    PET = Thornthwaite_PET(T_F,lat_list,count_loc,beg_row_col, ...
          lat_col);
else
    yr_uniqu=unique(yr_col);
    [ Days_in_month ] = Days_each_month( yr_uniqu );%%%Number of days in each month.
    if Elevation_file_load==1
        Elevation=Elevation_Data;
    else
        Elevation=0;
    end
    PET = Penman_PET(T_F,T_F_minmax,rsds,u_2,Spec_Hum,Alb_Data,Elevation,lat_list,count_loc,beg_row_col,Days_in_month);
    %%%       Note: If T_F_minmax is not available, then use [].
end
[ET,PR,R,RO,PRO,L,PL] = WaterBalance(AWC,PET,P,beg_row_col, ...
                        end_row_col,count_loc);

[Z_all] = Z_Index(P,PET,ET,PR,R,RO,PRO,L,PL,beg_row_mat,end_row_mat, ...
          count_loc,beg_calyr_mat,end_calyr_mat,beg_row_col,end_row_col);

[table] = PDSI_Central(Z_all,count_loc,beg_row_col,end_row_col, ...
          lat_col,yr_col,PET);

% Open the text file to which the table of values is to be written with 
% write permission.
fid = fopen('Palmer.txt','w');
fprintf(fid, '%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\r\n','Latitude','Year','PET','Z-Index','PPe','PX1','PX2','PX3','X','PDSI','PHDI');
fprintf(fid, '%10.6f %10.0f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\r\n',(table'));
fclose(fid);

% If the tool successfully executes the Palmer calculations, then print 
% text, including the directory and folder to which the output text file 
% was written, to the command window.
currentFolder = pwd;
screen_text_l1 = 'Success! The output file - "Palmer.txt" -';
screen_text_l2 = 'is located in';
fprintf('%s\n%s %s.\n',screen_text_l1,screen_text_l2,currentFolder)




function edit_const_alb_Callback(hObject, eventdata, handles)
% hObject    handle to edit_const_alb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_const_alb as text
%        str2double(get(hObject,'String')) returns contents of edit_const_alb as a double


% --- Executes during object creation, after setting all properties.
function edit_const_alb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_const_alb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_alb_file_Callback(hObject, eventdata, handles)
% hObject    handle to edit_alb_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_alb_file as text
%        str2double(get(hObject,'String')) returns contents of edit_alb_file as a double


% --- Executes during object creation, after setting all properties.
function edit_alb_file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_alb_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_browse.
function pushbutton_browse_Callback(hObject, eventdata, handles)
global  Data file_load

% hObject    handle to pushbutton_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.txt' , 'Select Albedo Text File');
Current_Dir{ 1 } = [ PathName FileName ];
set(handles.edit_alb_file,'string',Current_Dir{ 1 })
Data=load(Current_Dir{ 1 });
file_load=1;

% --- Executes on button press in pushbutton_OK.
function pushbutton_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global alb Data Alb_Data file_load

if alb==1
    Alb_Data=str2double(get(handles.edit_const_alb,'string'));
    h = msgbox(['Albedo Constant Value: ',num2str(Alb_Data)]);
    set(handles.pushbutton_run,'Enable','on')
elseif (alb==2) && (file_load==1)
    Alb_Data=Data;
    h = msgbox('Albedo File Successfully Loaded!');
    set(handles.pushbutton_run,'Enable','on')
else
    h = msgbox('Please Upload the Albedo File!');
    set(handles.pushbutton_run,'Enable','off')
end


% --- Executes when selected object is changed in uipanel_albedo.
function uipanel_albedo_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_albedo 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global alb
handles = guidata(hObject); 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'radiobutton1_const_Alb'
      %execute this code when hamon_radiobutton is selected
      alb = 1;
      set(handles.edit_const_alb,'Enable','on')
      set(handles.pushbutton_browse,'Enable','off')
      set(handles.pushbutton_run,'Enable','off')

    case 'radiobutton_alb_file'
      %execute this code when thorn_radiobutton is selected
      alb = 2;
      set(handles.pushbutton_browse,'Enable','on')
      set(handles.edit_const_alb,'Enable','off')
      set(handles.pushbutton_run,'Enable','off')


    otherwise
       % Code for when there is no match.
 
end


% --- Executes during object creation, after setting all properties.
function uipanel_albedo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_albedo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in uipanel7.
function uipanel7_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel7 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global mean_temp
handles = guidata(hObject); 
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'radiobutton_Mean_tmp'
      %execute this code when hamon_radiobutton is selected
      mean_temp = 1;
    case 'radiobutton_minmax_tmp'
      %execute this code when thorn_radiobutton is selected
      mean_temp = 0;

    otherwise
       % Code for when there is no match.
 
end


% --- Executes on button press in pushbutton_Elevation.
function pushbutton_Elevation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global  Elevation_Data Elevation_file_load

% hObject    handle to pushbutton_browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uigetfile('*.txt' , 'Select Elevation Text File');
Current_Dir{ 1 } = [ PathName FileName ];
Elevation_Data=load(Current_Dir{ 1 });
Elevation_file_load=1;
h = msgbox('Elevation File Successfully Loaded!');



