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
addpath(toolboxFolder);
addpath(btkFolder);

% -------------------------------------------------------------------------
% Asymptomatic subjects
% -------------------------------------------------------------------------
% List available C3D files
cd(subjectFolder)
c3dList = dir('*.c3d');

% Estimate gait events for each C3D file
for file = 2%1:length(c3dList)
    % Load 3D trajectories of markers and related frame rate and frame number
    cd(subjectFolder);
    c3dfile = c3dList(file).name;
    btkData = btkReadAcquisition(c3dfile);
    btkClearEvents(btkData);
    Rmarkers = btkGetMarkers(btkData);
    f = btkGetPointFrequency(btkData);
    n = btkGetPointFrameNumber(btkData);

    % Filter 3D trajectories of markers
    % Butterworth filter, zero-phase filter, 2nd order, low-pass, cut-off 10 Hz
    rightFootMarkersName = {'RDMT5','RHEE'};
    Fmarkers = [];
    for i = 1:length(rightFootMarkersName)
        [B,A] = butter(2,10/(f/2),'low');
        Fmarkers.(rightFootMarkersName{i}) = filtfilt(B, A, Rmarkers.(rightFootMarkersName{i}));
    end

    % Define gait and vertical axes
    tdiff = Fmarkers.(rightFootMarkersName{i})(end,:) - Fmarkers.(rightFootMarkersName{i})(1,:);
    [~,gaitAxis] = max(abs(tdiff));
    verticalAxis = 3; % Z axis vertical

    % Estimate gait events using Ghoussayni algorithm
    [FS,FO] = algorithm_Ghoussayni(Fmarkers,rightFootMarkersName,gaitAxis,verticalAxis,n,f);
    
    % Get gait events measured using forceplates and identified the related
    % estimated gait events
    [mFS,mFO,eFS,eFO] = setMeasuredEstimatedGaitEvents(btkData,FS,FO,n);
    
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
    btkWriteAcquisition(btkData,c3dfile);
end