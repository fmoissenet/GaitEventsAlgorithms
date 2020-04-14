% Turn warnings off
warning('off','MATLAB:interp1:NaNinY')
clear all

%% Set folders
% toolboxFolder = 'C:\Users\Ganglabor\Documents\automatisiertenEvents\Matlab\git_repo\GaitEventsAlgorithms';
% addpath(toolboxFolder);
% btkFolder     = 'C:\Users\FreslierM\Documents\MATLAB\btk';
% addpath(btkFolder);
subjectFolder = 'C:\Users\Ganglabor\Documents\automatisiertenEvents\Matlab\Daten\data_hemi';
% subjectFolder = 'C:\Users\Ganglabor\Documents\automatisiertenEvents\Matlab\Daten\NormSubjects';
exportFolder  = 'C:\Users\Ganglabor\Documents\automatisiertenEvents\Matlab\output\hemi';

%% hold the description of the data
T = readtable([subjectFolder '\hemi_CP_sofamea.xlsx']);
% T = readtable([subjectFolder '\norm_description.xlsx']);
left = T.left;
right = T.right;
affectedSide = T.LesionSide;
subjects = T.Subject_Study;
clear T

% -------------------------------------------------------------------------
%% Set Marker names
% -------------------------------------------------------------------------
FootMarkersName = { 'HEE' 'TPR' 'SITA' 'ANK'...
                    'PMT1' 'PMT5' 'CUN' ...
                    'DMT1' 'DMT5' 'TOE' 'HLX'};
R_FootMarkersName = { 'RHEE' 'RTPR' 'RSITA' 'RANK'...
                      'RPMT1' 'RPMT5' 'RCUN' ...
                      'RDMT1' 'RDMT5' 'RTOE' 'RHLX'};
L_FootMarkersName = { 'LHEE' 'LTPR' 'LSITA' 'LANK'...
                      'LPMT1' 'LPMT5' 'LCUN' ...
                      'LDMT1' 'LDMT5' 'LTOE' 'LHLX'};
FS_stat_hemi.nameMarkers = FootMarkersName;
FO_stat_hemi.nameMarkers = FootMarkersName;

%% subjects
for subj= 1:size(subjects,1)
    dataStr = subjects{subj,1}; % subject name
    switch(affectedSide{subj,1})
        case 'Left'
            infoSide = left{subj,1};
            textSplit_aff = regexp(infoSide,'Tr \d*| FP \d','match');
            trialCell = regexp(textSplit_aff{1,1},'\d*','match');
            trial_aff = str2num(trialCell{1,1});

            infoSide = right{subj,1};
            textSplit_unaff = regexp(infoSide,'Tr \d*| FP \d','match');
            trialCell = regexp(textSplit_unaff{1,1},'\d*','match');
            trial_unaff = str2num(trialCell{1,1});
        case 'Right'
            infoSide = right{subj,1};
            textSplit_aff = regexp(infoSide,'Tr \d*| FP \d','match');
            trialCell = regexp(textSplit_aff{1,1},'\d*','match');
            trial_aff = str2num(trialCell{1,1});

            infoSide = left{subj,1};
            textSplit_unaff = regexp(infoSide,'Tr \d*| FP \d','match');
            trialCell = regexp(textSplit_unaff{1,1},'\d*','match');
            trial_unaff = str2num(trialCell{1,1});
    end
    if ~isempty(dataStr)
        hemi = {'affected' 'unaffected'};
        for h = 1:length(hemi)
            FS_stat_hemi.(hemi{1,h}).nameSubjects{subj,1} = dataStr;
            FO_stat_hemi.(hemi{1,h}).nameSubjects{subj,1} = dataStr;
            % informations of the side
            switch(hemi{1,h})
                case 'affected'
                    switch(affectedSide{subj,1})
                        case 'Left'
                            side = 'L';
                        case 'Right'
                            side = 'R';
                    end
                case 'unaffected'
                    switch(affectedSide{subj,1})
                        case 'Left'
                            side = 'R';
                        case 'Right'
                            side = 'L';
                    end
            end

            % path of the c3d
            switch(hemi{1,h})
                case 'affected'
                    if trial_aff == trial_unaff
                        c3dPath = [subjectFolder '\' dataStr '.c3d'];
                    else
                        c3dPath = [subjectFolder '\' dataStr '_affected.c3d'];
                    end
                    fpCell = regexp(textSplit_aff{1,2},'\d*','match');
                case 'unaffected'
                    if trial_aff == trial_unaff
                        c3dPath = [subjectFolder '\' dataStr '.c3d'];
                    else
                        c3dPath = [subjectFolder '\' dataStr '_unaffected.c3d'];
                    end
                    fpCell = regexp(textSplit_unaff{1,2},'\d*','match');
            end

%             clear infoSide textSplit_aff textSplit_unaff trialCell trial_aff trial_unaff
            
            %% Load 3D trajectories of markers and related frame rate and frame number
            btkData = btkReadAcquisition(c3dPath);
            Markers = btkGetMarkers(btkData);
            f = btkGetPointFrequency(btkData);
            n = btkGetPointFrameNumber(btkData);
            
            %% Define gait and vertical axes
            markers_corrected = f_rotCoordinateSystem(Markers);
            gaitAxis = 1; % X axis
            verticalAxis = 3; % Z axis vertical
            
            %% Filter 3D trajectories of markers
            % Butterworth filter, zero-phase filter, 2nd order, low-pass, cut-off 10 Hz
            [B,A] = butter(2,10/(f/2),'low');
            nonfiltered_markers = [];
            filtered_markers = [];
            switch(side)
                case 'L'
                    for i = 1:length(L_FootMarkersName)
                        nonfiltered_markers(:,:,i) = markers_corrected.(L_FootMarkersName{1,i});
                        filtered_markers(:,:,i) = filtfilt(B, A, markers_corrected.(L_FootMarkersName{1,i}));
                    end
                case 'R'
                    for i = 1:length(R_FootMarkersName)
                        nonfiltered_markers(:,:,i) = markers_corrected.(R_FootMarkersName{1,i});
                        filtered_markers(:,:,i) = filtfilt(B, A, markers_corrected.(R_FootMarkersName{1,i}));
                    end
            end
            pelvicMk.filtLASI = filtfilt(B, A, markers_corrected.LASI);
            pelvicMk.filtRASI = filtfilt(B, A, markers_corrected.RASI);
            pelvicMk.filtLPSI = filtfilt(B, A, markers_corrected.LPSI);
            pelvicMk.filtRPSI = filtfilt(B, A, markers_corrected.RPSI);
            pelvicMk.filtSACR = filtfilt(B, A, markers_corrected.SACR);
            
            %% Get gait events measured using forceplates
            [mFS,mFO] = getMeasuredGaitEvents(btkData,n,str2num(fpCell{1,1}));
% figure
% hold on
            for iMk = 1:length(FootMarkersName)
                %% Ghoussayni
                [FS_Gh,FO_Gh] = Ghoussayni(filtered_markers(:,:,iMk),gaitAxis,verticalAxis,n,f);
                % Get the estimated gait events related to the forceplate
                [eFS,eFO] = getEstimatedGaitEvents(btkData,FS_Gh,FO_Gh,f,side);
                % !!! the FS and FO are given in frames! some measurements are done with 300Hz and not 150Hz!
                % give all the results in 150Hz : divided by 2 = 150/300
                if isempty(mFS) == 0 && isempty(mFO) == 0 && isempty(eFS) == 0 && isempty(eFO) == 0
                    if f==300
                        FS_stat_hemi.(hemi{1,h}).Ghoussayni.absError(subj,iMk) = floor((eFS-mFS)/2);
                        FO_stat_hemi.(hemi{1,h}).Ghoussayni.absError(subj,iMk) = floor((eFO-mFO)/2);
                    else
                        FS_stat_hemi.(hemi{1,h}).Ghoussayni.absError(subj,iMk) = eFS-mFS;
                        FO_stat_hemi.(hemi{1,h}).Ghoussayni.absError(subj,iMk) = eFO-mFO;
                    end
                end

                %% Zeni
                [FS_Zeni,FO_Zeni] = Zeni(filtered_markers(:,:,iMk),pelvicMk,gaitAxis,f);%,verticalAxis,n);
                % Get the estimated gait events related to the forceplate
                [eFS,eFO] = getEstimatedGaitEvents(btkData,FS_Zeni,FO_Zeni,f,side);
                % !!! the FS and FO are given in frames! some measurements are done with 300Hz and not 150Hz!
                % give all the results in 150Hz : divided by 2 = 150/300
                if isempty(mFS) == 0 && isempty(mFO) == 0 && isempty(eFS) == 0 && isempty(eFO) == 0
                    if f==300
                        FS_stat_hemi.(hemi{1,h}).Zeni.absError(subj,iMk) = floor((eFS-mFS)/2);
                        FO_stat_hemi.(hemi{1,h}).Zeni.absError(subj,iMk) = floor((eFO-mFO)/2);
                    else
                        FS_stat_hemi.(hemi{1,h}).Zeni.absError(subj,iMk) = eFS-mFS;
                        FO_stat_hemi.(hemi{1,h}).Zeni.absError(subj,iMk) = eFO-mFO;
                    end
                end
                
                %% Desailly
                [FS_Des,FO_Des] = Desailly(nonfiltered_markers(:,:,1),nonfiltered_markers(:,:,iMk),gaitAxis,verticalAxis,f);
                % Get the estimated gait events related to the forceplate
                [eFS,eFO] = getEstimatedGaitEvents(btkData,FS_Des,FO_Des,f,side);
                % !!! the FS and FO are given in frames! some measurements are done with 300Hz and not 150Hz!
                % give all the results in 150Hz : divided by 2 = 150/300
                if isempty(mFS) == 0 && isempty(mFO) == 0 && isempty(eFS) == 0 && isempty(eFO) == 0
                    if f==300
                        FS_stat_hemi.(hemi{1,h}).Desailly.absError(subj,iMk) = floor((eFS-mFS)/2);
                        FO_stat_hemi.(hemi{1,h}).Desailly.absError(subj,iMk) = floor((eFO-mFO)/2);
                    else
                        FS_stat_hemi.(hemi{1,h}).Desailly.absError(subj,iMk) = eFS-mFS;
                        FO_stat_hemi.(hemi{1,h}).Desailly.absError(subj,iMk) = eFO-mFO;
                    end
                end
                
                %% Hsue
                [FS_Hsue,FO_Hsue] = Hsue(filtered_markers(:,:,iMk),pelvicMk,gaitAxis,f);
                % Get the estimated gait events related to the forceplate
                [eFS,eFO] = getEstimatedGaitEvents(btkData,FS_Hsue,FO_Hsue,f,side);
                % !!! the FS and FO are given in frames! some measurements are done with 300Hz and not 150Hz!
                % give all the results in 150Hz : divided by 2 = 150/300
                if isempty(mFS) == 0 && isempty(mFO) == 0 && isempty(eFS) == 0 && isempty(eFO) == 0
                    if f==300
                        FS_stat_hemi.(hemi{1,h}).Hsue.absError(subj,iMk) = floor((eFS-mFS)/2);
                        FO_stat_hemi.(hemi{1,h}).Hsue.absError(subj,iMk) = floor((eFO-mFO)/2);
                    else
                        FS_stat_hemi.(hemi{1,h}).Hsue.absError(subj,iMk) = eFS-mFS;
                        FO_stat_hemi.(hemi{1,h}).Hsue.absError(subj,iMk) = eFO-mFO;
                    end
                end

                %% Hreljac
                [FS_Hreljac,FO_Hreljac] = Hreljac(filtered_markers(:,:,iMk),pelvicMk,gaitAxis,verticalAxis,f);
                % Get the estimated gait events related to the forceplate
                [eFS,eFO] = getEstimatedGaitEvents(btkData,FS_Hreljac,FO_Hreljac,f,side);
                % !!! the FS and FO are given in frames! some measurements are done with 300Hz and not 150Hz!
                % give all the results in 150Hz : divided by 2 = 150/300
                if isempty(mFS) == 0 && isempty(mFO) == 0 && isempty(eFS) == 0 && isempty(eFO) == 0
                    if f==300
                        FS_stat_hemi.(hemi{1,h}).Hreljac.absError(subj,iMk) = (eFS-mFS)/2;
                        FO_stat_hemi.(hemi{1,h}).Hreljac.absError(subj,iMk) = (eFO-mFO)/2;
                    else
                        FS_stat_hemi.(hemi{1,h}).Hreljac.absError(subj,iMk) = eFS-mFS;
                        FO_stat_hemi.(hemi{1,h}).Hreljac.absError(subj,iMk) = eFO-mFO;
                    end
                end

            end % Markers
            
            btkDeleteAcquisition(btkData);
        end % hemi (affected,unaffected)
    end
end % subjects

FS_stat_hemi.affected.Ghoussayni.Mean = mean(abs(FS_stat_hemi.affected.Ghoussayni.absError));
FS_stat_hemi.affected.Ghoussayni.Std = std(abs(FS_stat_hemi.affected.Ghoussayni.absError));
FO_stat_hemi.affected.Ghoussayni.Mean = mean(abs(FO_stat_hemi.affected.Ghoussayni.absError));
FO_stat_hemi.affected.Ghoussayni.Std = std(abs(FO_stat_hemi.affected.Ghoussayni.absError));

FS_stat_hemi.affected.Zeni.Mean = mean(abs(FS_stat_hemi.affected.Zeni.absError));
FS_stat_hemi.affected.Zeni.Std = std(abs(FS_stat_hemi.affected.Zeni.absError));
FO_stat_hemi.affected.Zeni.Mean = mean(abs(FO_stat_hemi.affected.Zeni.absError));
FO_stat_hemi.affected.Zeni.Std = std(abs(FO_stat_hemi.affected.Zeni.absError));

FS_stat_hemi.affected.Desailly.Mean = mean(abs(FS_stat_hemi.affected.Desailly.absError));
FS_stat_hemi.affected.Desailly.Std = std(abs(FS_stat_hemi.affected.Desailly.absError));
FO_stat_hemi.affected.Desailly.Mean = mean(abs(FO_stat_hemi.affected.Desailly.absError));
FO_stat_hemi.affected.Desailly.Std = std(abs(FO_stat_hemi.affected.Desailly.absError));

FS_stat_hemi.affected.Hsue.Mean = mean(abs(FS_stat_hemi.affected.Hsue.absError));
FS_stat_hemi.affected.Hsue.Std = std(abs(FS_stat_hemi.affected.Hsue.absError));
FO_stat_hemi.affected.Hsue.Mean = mean(abs(FO_stat_hemi.affected.Hsue.absError));
FO_stat_hemi.affected.Hsue.Std = std(abs(FO_stat_hemi.affected.Hsue.absError));

FS_stat_hemi.affected.Hreljac.Mean = mean(abs(FS_stat_hemi.affected.Hreljac.absError));
FS_stat_hemi.affected.Hreljac.Std = std(abs(FS_stat_hemi.affected.Hreljac.absError));
FO_stat_hemi.affected.Hreljac.Mean = mean(abs(FO_stat_hemi.affected.Hreljac.absError));
FO_stat_hemi.affected.Hreljac.Std = std(abs(FO_stat_hemi.affected.Hreljac.absError));

FS_stat_hemi.unaffected.Ghoussayni.Mean = mean(abs(FS_stat_hemi.unaffected.Ghoussayni.absError));
FS_stat_hemi.unaffected.Ghoussayni.Std = std(abs(FS_stat_hemi.unaffected.Ghoussayni.absError));
FO_stat_hemi.unaffected.Ghoussayni.Mean = mean(abs(FO_stat_hemi.unaffected.Ghoussayni.absError));
FO_stat_hemi.unaffected.Ghoussayni.Std = std(abs(FO_stat_hemi.unaffected.Ghoussayni.absError));

FS_stat_hemi.unaffected.Zeni.Mean = mean(abs(FS_stat_hemi.unaffected.Zeni.absError));
FS_stat_hemi.unaffected.Zeni.Std = std(abs(FS_stat_hemi.unaffected.Zeni.absError));
FO_stat_hemi.unaffected.Zeni.Mean = mean(abs(FO_stat_hemi.unaffected.Zeni.absError));
FO_stat_hemi.unaffected.Zeni.Std = std(abs(FO_stat_hemi.unaffected.Zeni.absError));

FS_stat_hemi.unaffected.Desailly.Mean = mean(abs(FS_stat_hemi.unaffected.Desailly.absError));
FS_stat_hemi.unaffected.Desailly.Std = std(abs(FS_stat_hemi.unaffected.Desailly.absError));
FO_stat_hemi.unaffected.Desailly.Mean = mean(abs(FO_stat_hemi.unaffected.Desailly.absError));
FO_stat_hemi.unaffected.Desailly.Std = std(abs(FO_stat_hemi.unaffected.Desailly.absError));

FS_stat_hemi.unaffected.Hsue.Mean = mean(abs(FS_stat_hemi.unaffected.Hsue.absError));
FS_stat_hemi.unaffected.Hsue.Std = std(abs(FS_stat_hemi.unaffected.Hsue.absError));
FO_stat_hemi.unaffected.Hsue.Mean = mean(abs(FO_stat_hemi.unaffected.Hsue.absError));
FO_stat_hemi.unaffected.Hsue.Std = std(abs(FO_stat_hemi.unaffected.Hsue.absError));

FS_stat_hemi.unaffected.Hreljac.Mean = mean(abs(FS_stat_hemi.unaffected.Hreljac.absError));
FS_stat_hemi.unaffected.Hreljac.Std = std(abs(FS_stat_hemi.unaffected.Hreljac.absError));
FO_stat_hemi.unaffected.Hreljac.Mean = mean(abs(FO_stat_hemi.unaffected.Hreljac.absError));
FO_stat_hemi.unaffected.Hreljac.Std = std(abs(FO_stat_hemi.unaffected.Hreljac.absError));

save([exportFolder '\FS_stat_hemi.mat'],'FS_stat_hemi');
save([exportFolder '\FO_stat_hemi.mat'],'FO_stat_hemi');

