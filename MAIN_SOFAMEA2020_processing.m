% -------------------------------------------------------------------------
% Initialisation
% -------------------------------------------------------------------------
clearvars;
close all;
clc;

% Set folders
toolboxFolder = 'C:\Users\moissene\Documents\Professionnel\routines\github\GaitEventsAlgorithms';
btkFolder     = 'C:\Users\moissene\Documents\Professionnel\routines\btk';
patientFolder = 'C:\Users\moissene\Documents\Professionnel\publications\communications\2020\SOFAMEA2020\Patients';
subjectFolder = 'C:\Users\moissene\Documents\Professionnel\publications\communications\2020\SOFAMEA2020\Healthy subjects';
exportFolder  = 'C:\Users\moissene\Documents\Professionnel\publications\communications\2020\SOFAMEA2020\Outputs';
addpath(toolboxFolder);
addpath(btkFolder);

% -------------------------------------------------------------------------
% Set conditions
% -------------------------------------------------------------------------
conf  = {'R_conf1' 'R_conf2' 'R_conf3' 'R_conf4' ...
         'R_confCGM10' 'R_confCGM24a' 'R_confCGM24b' 'R_confCGM24c' ...
         'R_confLEARa' 'R_confLEARb' 'R_confLEARc'};
acase = {'case1' 'case2' 'case3'};

for iconf = 1:length(conf)
for icase = 1:length(acase)
% -------------------------------------------------------------------------
% Process data
% -------------------------------------------------------------------------
% List available C3D files
cd(patientFolder)
c3dList = dir('*.c3d');

% Estimate gait events for each C3D file
for file = 1:length(c3dList)
    % Load 3D trajectories of markers and related frame rate and frame number
    cd(patientFolder);
    c3dfile = c3dList(file).name;
    btkData = btkReadAcquisition(c3dfile);
    btkClearEvents(btkData);
    Markers = btkGetMarkers(btkData);
    f = btkGetPointFrequency(btkData);
    n = btkGetPointFrameNumber(btkData);

    % Define set of markers    
    % Configuration 1: 1 marker for hind and forefoot
    if strcmp(conf{iconf},'R_conf1')
        hindFootMarkersName = {'RHEE'};
        foreFootMarkersName = {'RDMT5'};
    % Configuration 2: all markers for hind and forefoot
    elseif strcmp(conf{iconf},'R_conf2')
        hindFootMarkersName = {'RHEE','RTPR','RSITA'};
        foreFootMarkersName = {'RDMT1','RDMT5','RTOE','RHLX'};
    % Configuration 3: all markers for hindfoot, all markers + middle foot
    % markers for forefoot
    elseif strcmp(conf{iconf},'R_conf3')
        hindFootMarkersName = {'RHEE','RTPR','RSITA'};
        foreFootMarkersName = {'RDMT1','RDMT5','RTOE','RHLX','RPMT1','RPMT5','RCUN'};
    % Configuration 4: all markers + middle foot markers for hindfoot, all
    % markers for forefoot
    elseif strcmp(conf{iconf},'R_conf4')
        hindFootMarkersName = {'RHEE','RTPR','RSITA','RPMT1','RPMT5','RCUN'};
        foreFootMarkersName = {'RDMT1','RDMT5','RTOE','RHLX'};
    % -- Right CGM 1.0
    elseif strcmp(conf{iconf},'R_confCGM10')
        hindFootMarkersName = {'RHEE'};
        foreFootMarkersName = {'RTOE'};
    % -- Right CGM 2.4 (middle foot markers not used)
    elseif strcmp(conf{iconf},'R_confCGM24a')
        hindFootMarkersName = {'RHEE'};
        foreFootMarkersName = {'RDMT1','RDMT5'};
    % -- Right CGM 2.4 (middle foot markers used as forefoot marker)
    elseif strcmp(conf{iconf},'R_confCGM24b')
        hindFootMarkersName = {'RHEE'};
        foreFootMarkersName = {'RDMT1','RDMT5','RCUN'};
    % -- Right CGM 2.4 (middle foot markers used as hindfoot marker)
    elseif strcmp(conf{iconf},'R_confCGM24c')
        hindFootMarkersName = {'RHEE','RCUN'};
        foreFootMarkersName = {'RDMT1','RDMT5'};
    % -- Right Leardini (middle foot markers not used)
    elseif strcmp(conf{iconf},'R_confLEARa')
        hindFootMarkersName = {'RHEE','RTPR'};
        foreFootMarkersName = {'RDMT1','RDMT5','RTOE','RHLX'};
    % -- Right Leardini (middle foot markers used as forefoot marker)
    elseif strcmp(conf{iconf},'R_confLEARb')
        hindFootMarkersName = {'RHEE','RTPR'};
        foreFootMarkersName = {'RDMT1','RDMT5','RTOE','RHLX','RPMT1','RPMT5','RCUN'};
    % -- Right Leardini (middle foot markers used as hindfoot marker)
    elseif strcmp(conf{iconf},'R_confLEARc')
        hindFootMarkersName = {'RHEE','RTPR','RPMT1','RPMT5','RCUN'};
        foreFootMarkersName = {'RDMT1','RDMT5','RTOE','RHLX'};
    end
    
    % Filter 3D trajectories of markers
    % Butterworth filter, zero-phase filter, 2nd order, low-pass, cut-off 10 Hz    
    [B,A] = butter(2,10/(f/2),'low');
    Hmarkers = [];
    for i = 1:length(hindFootMarkersName)
        Hmarkers(:,:,i) = filtfilt(B, A, Markers.(hindFootMarkersName{i}));
    end
    Fmarkers = [];
    for i = 1:length(foreFootMarkersName)
        Fmarkers(:,:,i) = filtfilt(B, A, Markers.(foreFootMarkersName{i}));
    end

    % Define gait and vertical axes
    tdiff = Hmarkers(end,:,1) - Hmarkers(1,:,1);
    [~,gaitAxis] = max(abs(tdiff));
    verticalAxis = 3; % Z axis vertical

    % Estimate gait events using Ghoussayni algorithm
    [FS,FO] = algorithm_Ghoussayni(Hmarkers,Fmarkers,acase{icase},gaitAxis,verticalAxis,n,f);
    
    % Get gait events measured using forceplates and identified the related
    % estimated gait events
    [mFS,mFO,eFS,eFO] = setMeasuredEstimatedGaitEvents(btkData,FS,FO,n,f);
    
    if isempty(mFS) == 0 && isempty(mFO) == 0 && isempty(eFS) == 0 && ...
       isempty(eFO) == 0
        % Store results in the C3D file
        btkAppendEvent(btkData,'eFS',(eFS+btkGetFirstFrame(btkData))/f,'Right');
        btkAppendEvent(btkData,'eFO',(eFO+btkGetFirstFrame(btkData))/f,'Right');
        btkAppendEvent(btkData,'mFS',(mFS+btkGetFirstFrame(btkData))/f,'Right');
        btkAppendEvent(btkData,'mFO',(mFO+btkGetFirstFrame(btkData))/f,'Right');

        % Store errors (frames)
        info.format = 'Integer';
        info.values = eFS-mFS;
        btkAppendMetaData(btkData,'EventDetectionErrors','FS',info);
        info.format = 'Integer';
        info.values = eFO-mFO;
        btkAppendMetaData(btkData,'EventDetectionErrors','FO',info);

        % Update C3D file
        cd(exportFolder);
        btkWriteAcquisition(btkData,[conf{iconf},'_',acase{icase},'_',c3dfile]);
    end
end
end
end