% -------------------------------------------------------------------------
% Initialisation
% -------------------------------------------------------------------------
clearvars;
close all;
clc;

% Set folders
toolboxFolder = 'C:\Users\moissene\Documents\Professionnel\routines\github\GaitEventsAlgorithms';
btkFolder     = 'C:\Users\moissene\Documents\Professionnel\routines\btk';
patientFolder = 'C:\Users\moissene\Documents\Professionnel\publications\articles\1- en cours\Freslier - Gait events\C3D Files';
subjectFolder = 'C:\Users\moissene\Documents\Professionnel\routines\github\GaitEventsAlgorithms\NormSubjects';
exportFolder  = 'C:\Users\moissene\Documents\Professionnel\publications\communications\2020\SOFAMEA2020\Data';
addpath(toolboxFolder);
addpath(btkFolder);

% -------------------------------------------------------------------------
% Set conditions
% -------------------------------------------------------------------------
conf = 'R_confLEARb';
acase = 'case3';

% -------------------------------------------------------------------------
% Asymptomatic subjects
% -------------------------------------------------------------------------
% List available C3D files
cd(subjectFolder)
c3dList = dir('*.c3d');

% Estimate gait events for each C3D file
for file = 1%:length(c3dList)
    % Load 3D trajectories of markers and related frame rate and frame number
    cd(subjectFolder);
    c3dfile = c3dList(file).name;
    btkData = btkReadAcquisition(c3dfile);
    btkClearEvents(btkData);
    Markers = btkGetMarkers(btkData);
    f = btkGetPointFrequency(btkData);
    n = btkGetPointFrameNumber(btkData);

    % Define set of markers    
    % Configuration 1: 1 marker for hind and forefoot
    if strcmp(conf,'R_conf1')
        hindFootMarkersName = {'RHEE'};
        foreFootMarkersName = {'RDMT5'};
    % Configuration 2: all markers for hind and forefoot
    elseif strcmp(conf,'R_conf2')
        hindFootMarkersName = {'RHEE','RTPR','RSITA'};
        foreFootMarkersName = {'RDMT1','RDMT5','RTOE','RHLX'};
    % Configuration 3: all markers for hindfoot, all markers + middle foot
    % markers for forefoot
    elseif strcmp(conf,'R_conf3')
        hindFootMarkersName = {'RHEE','RTPR','RSITA'};
        foreFootMarkersName = {'RDMT1','RDMT5','RTOE','RHLX','RPMT1','RPMT5','RCUN'};
    % Configuration 4: all markers + middle foot markers for hindfoot, all
    % markers for forefoot
    elseif strcmp(conf,'R_conf4')
        hindFootMarkersName = {'RHEE','RTPR','RSITA','RPMT1','RPMT5','RCUN'};
        foreFootMarkersName = {'RDMT1','RDMT5','RTOE','RHLX'};
    % -- Right CGM 1.0
    elseif strcmp(conf,'R_confCGM10')
        hindFootMarkersName = {'RHEE'};
        foreFootMarkersName = {'RDMT5'};
    % -- Right CGM 2.4 (middle foot markers not used)
    elseif strcmp(conf,'R_confCGM24a')
        hindFootMarkersName = {'RHEE'};
        foreFootMarkersName = {'RDMT1','RDMT5'};
    % -- Right CGM 2.4 (middle foot markers used as forefoot marker)
    elseif strcmp(conf,'R_confCGM24b')
        hindFootMarkersName = {'RHEE'};
        foreFootMarkersName = {'RDMT1','RDMT5','RCUN'};
    % -- Right CGM 2.4 (middle foot markers used as hindfoot marker)
    elseif strcmp(conf,'R_confCGM24c')
        hindFootMarkersName = {'RHEE','RCUN'};
        foreFootMarkersName = {'RDMT1','RDMT5'};
    % -- Right Leardini (middle foot markers not used)
    elseif strcmp(conf,'R_confLEARa')
        hindFootMarkersName = {'RHEE','RPMT5'};
        foreFootMarkersName = {'RDMT1','RDMT5','RTOE','RHLX'};
    % -- Right Leardini (middle foot markers used as forefoot marker)
    elseif strcmp(conf,'R_confLEARb')
        hindFootMarkersName = {'RHEE','RPMT5'};
        foreFootMarkersName = {'RDMT1','RDMT5','RTOE','RHLX','RPMT1','RPMT5','RCUN'};
    % -- Right Leardini (middle foot markers used as hindfoot marker)
    elseif strcmp(conf,'R_confLEARc')
        hindFootMarkersName = {'RHEE','RPMT5','RPMT1','RPMT5','RCUN'};
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
    [FS,FO] = algorithm_Ghoussayni(Hmarkers,Fmarkers,acase,gaitAxis,verticalAxis,n,f);
    
    % Get gait events measured using forceplates and identified the related
    % estimated gait events
    [mFS,mFO,eFS,eFO] = setMeasuredEstimatedGaitEvents(btkData,FS,FO,n,f);
    
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
    btkWriteAcquisition(btkData,[conf,'_',acase,'_',c3dfile]);
end