% This file contains all the variables that should be modified at beginning of every cruise
%-----------------------------
 
struct_levels_to_print(0);

warning('off');


%-----------------------------
CRUISE = 'AMT29';

%-----------------------------
% Variables to be changed during cruise according to specific setups and user needs
%
% Dates
% Set date interval to be processed (format is 'yyyymmdd')
% (this will have to be changed each time the scripts are run)
inidate = '20191121';
enddate = '20191122';

% Hour of the day for which Wapped files are searched
% (day is not processed if a file for the specified hour is not found)
% Format is '0HH'
WAPhour = '001';

% Underway subdirectory where to find special wapped data
% Leave with simple / if no special case
% UWAY_WAP_SUBDIR = 'withACS167/'; 
% UWAY_WAP_SUBDIR = 'withAC9277/'; 
UWAY_WAP_SUBDIR = '/'; 

% Parameters specific for Underway plotting/processing
% (this will change depending on specific section fo the cruise)
% Setup to automatically change based on UWAY_WAP_SUBDIR
%
% Implemented instruments to selct from are 
% {'ctd','acs','bb3','cstar','acs2','ac9','clam'}
if UWAY_WAP_SUBDIR == 'withACS167/'
    dh8_instruments = {'ctd','acs','acs2','cstar','bb3','clam'};
    % Ports must corresponds to same ports as in dh8_instruments
    dh8_ports = {2,7,5,6,1,4}; 
    % Serial numbers are mainly needed for acs and ac9 config files, leave blank for other instruments
    dh8_serialnumber = {[],122,167,[],1173,[]}; 
elseif UWAY_WAP_SUBDIR == 'withAC9277/'
    dh8_instruments = {'ctd','acs','ac9','bb3','clam'};
    % Ports must corresponds to same ports as in dh8_instruments
    dh8_ports = {2,7,5,1,4}; 
    % Serial numbers are mainly needed for acs and ac9 config files, leave blank for other instruments
    dh8_serialnumber = {[],122,277,1173,[]}; 
elseif UWAY_WAP_SUBDIR == '/'
    dh8_instruments = {'ctd','acs','cstar','bb3','clam'};
    % Ports must corresponds to same ports as in dh8_instruments
    dh8_ports = {2,7,6,1,4}; 
    % Serial numbers are mainly needed for acs and ac9 config files, leave blank for other instruments
    dh8_serialnumber = {[],122,[],1173,[]}; 
endif
%-----------------------------

%-----------------------------
% Paths
%MAIN_PATH = '/data/datasets/cruise_data/active/AMT29/Public_Read_Only_Copy/DY110_Public/Optics_group/';
MAIN_PATH = '/cruise_data/AMT29/Public_Read_Only_Copy/DY110_Public/Optics_group/';
%MAIN_PATH = [MAIN_PATH,'/Data/',CRUISE,'/'];     % Root directory for current AMT cruise
PATH_DATA = [MAIN_PATH,'Data/'];        % Directory with all raw and wapped data
PATH_SOURCE = [MAIN_PATH,'Source/'];% Directory with all source code
OUT_PROC = [MAIN_PATH,'Processed/'];    % Output directory for processed oct and mat files
OUT_FIGS = [MAIN_PATH,'Figures/'];      % Output directory for figures

addpath([PATH_SOURCE, 'Octave_functions']);
addpath([PATH_SOURCE]);
%-----------------------------
% Each directory will contain a series of subdirectories for each instrument
% (e.g. Underway, Optics_rig, BB3_ctd etc. etc.)
OPTIC_DIR = 'Optics_rig/';
UWAY_DIR = 'Underway/';
BB3_DIR = 'BB3_ctd/';
CTD_DIR = 'Ship_CTD/';
% Specific data subdirectories
DATA_WAPPED = 'WAP_extracted/';
DATA_RAW = 'Raw/';
DATA_FLOW = 'Flow/';
DATA_WITH_BB3 = 'with_BB3/';
%-----------------------------
% Ship's system directories
PATH_SHIP = [PATH_DATA,'Ship_uway/'];
PATH_TS = [PATH_SHIP,'SURFMETV3/']; % Directory with ship underway ctd
PATH_GPS = [PATH_SHIP,'GPS/']; % Directory with ship gps
ship_uway_fname = '*Surf-DY-SM_DY1*'; % Name of ship underway system files
%-----------------------------

%-----------------------------
% Parameters specific for Optics rig plotting/processing
%
% Wether cdt is saved as ASCII format (false for AMT26; true for AMT27)
ctdASCII = true;
% Limits for temperature and salinity profiles
Tlim = [0 20];
Slim = [33 35];
% Limits for absorption and attenuation profiles
alim = [0.1 0.4];
clim = [0.05 0.6];
chlac9lim = [0 5];
%-----------------------------

% Processors to be used by parcellfun in run_step1par.m
NProc = 4;
% Name of processed file to be saved
fproc_name = ['optics_' lower(CRUISE) '_'];
% Limits for time-series plots
acs_raw_lim = [-0.03 0.1]; % acs
flow_lim = [20 45];        % flow rate
bb3_lim = [50 140];       % backscattering
SST_lim = [15 20];         % CTD temperature
SSS_lim = [35 36.5];        % CTD salinity
% Limits for step2 time series
acs_lim = [-0.05 0.3];
bb_opt_lim = [70 150];
cstar_lim = [0.75 0.85];
spectra_alim = [0.03];
spectra_clim = [1];
chl_lim = [0.01 5];

%-----------------------------
% Parameters specific for BB3 plotting/processing
%
% Limits for bb3 profiles
bb3lim = [50 300];
%-----------------------------

%-----------------------------
% Parameters specific for underway transect
%
latlim = 54;
trans_SST = [01 30];
trans_SSS = [33 38];
trans_chl = [0.01 5];
trans_cp650 = [0.01 0.3];
%-----------------------------

%-----------------------------
% useful functions
movavg = inline("filter(1/mavgwd*ones(1, mavgwd), 1, x)", "x", "mavgwd"); % this is equivalent to a standard moving average with a window size of mavgwd
%-----------------------------
