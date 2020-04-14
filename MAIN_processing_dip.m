% Turn warnings off
warning('off','MATLAB:interp1:NaNinY')
clear all

%% Set folders
% toolboxFolder = 'C:\Users\Ganglabor\Documents\automatisiertenEvents\Matlab\git_repo\GaitEventsAlgorithms';
% addpath(toolboxFolder);
% btkFolder     = 'C:\Users\FreslierM\Documents\MATLAB\btk';
% addpath(btkFolder);
subjectFolder = 'C:\Users\Ganglabor\Documents\automatisiertenEvents\Matlab\Daten\data_diparese';
% subjectFolder = 'C:\Users\Ganglabor\Documents\automatisiertenEvents\Matlab\Daten\NormSubjects';
exportFolder  = 'C:\Users\Ganglabor\Documents\automatisiertenEvents\Matlab\output\norm';

%% hold the description of the data
T = readtable([subjectFolder '\di_CP_sofamea.xlsx']);
% T = readtable([subjectFolder '\norm_description.xlsx']);
left = T.left;
right = T.right;
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
FS_stat_dip.nameMarkers = FootMarkersName;
FO_stat_dip.nameMarkers = FootMarkersName;

%% subjects
for subj= 1:size(subjects,1)
    dataStr = subjects{subj,1}; % subject name
    if ~isempty(dataStr)
        side = {'L' 'R'};
        for s = 1:length(side)
            FS_stat_dip.nameSubjects{(subj-1)*2+s,1} = [dataStr '_' side{1,s}];
            FO_stat_dip.nameSubjects{(subj-1)*2+s,1} = [dataStr '_' side{1,s}];
            % informations of the side from excel sheet
            switch(side{1,s})
                case 'L'
                    infoSide = left{subj,1};
                case 'R'
                    infoSide = right{subj,1};
            end
            textSplit = regexp(infoSide,'Tr \d*| FP \d','match');
            trialCell = regexp(textSplit{1,1},'\d*','match');
            trial = str2num(trialCell{1,1});
            fpCell = regexp(textSplit{1,2},'\d*','match');
            
            %% path of the c3d
            if trial <= 9
                c3dPath = [subjectFolder '\' dataStr '_0' trialCell{1,1} '.c3d'];
            else
                c3dPath = [subjectFolder '\' dataStr '_' trialCell{1,1} '.c3d'];
            end

            clear textSplit trialCell trial infoSide
            
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
            switch(side{1,s})
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
                [eFS,eFO] = getEstimatedGaitEvents(btkData,FS_Gh,FO_Gh,f,side{1,s});
                % !!! the FS and FO are given in frames! some measurements are done with 300Hz and not 150Hz!
                % give all the results in 150Hz : divided by 2 = 150/300
                if isempty(mFS) == 0 && isempty(mFO) == 0 && isempty(eFS) == 0 && isempty(eFO) == 0
                    if f==300
                        FS_stat_dip.Ghoussayni.absError((subj-1)*2+s,iMk) = floor((eFS-mFS)/2);
                        FO_stat_dip.Ghoussayni.absError((subj-1)*2+s,iMk) = floor((eFO-mFO)/2);
                    else
                        FS_stat_dip.Ghoussayni.absError((subj-1)*2+s,iMk) = eFS-mFS;
                        FO_stat_dip.Ghoussayni.absError((subj-1)*2+s,iMk) = eFO-mFO;
                    end
                end

                %% Zeni
                [FS_Zeni,FO_Zeni] = Zeni(filtered_markers(:,:,iMk),pelvicMk,gaitAxis,f);%,verticalAxis,n);
                % Get the estimated gait events related to the forceplate
                [eFS,eFO] = getEstimatedGaitEvents(btkData,FS_Zeni,FO_Zeni,f,side{1,s});
                % !!! the FS and FO are given in frames! some measurements are done with 300Hz and not 150Hz!
                % give all the results in 150Hz : divided by 2 = 150/300
                if isempty(mFS) == 0 && isempty(mFO) == 0 && isempty(eFS) == 0 && isempty(eFO) == 0
                    if f==300
                        FS_stat_dip.Zeni.absError((subj-1)*2+s,iMk) = floor((eFS-mFS)/2);
                        FO_stat_dip.Zeni.absError((subj-1)*2+s,iMk) = floor((eFO-mFO)/2);
                    else
                        FS_stat_dip.Zeni.absError((subj-1)*2+s,iMk) = eFS-mFS;
                        FO_stat_dip.Zeni.absError((subj-1)*2+s,iMk) = eFO-mFO;
                    end
                end
                
                %% Desailly
                [FS_Des,FO_Des] = Desailly(nonfiltered_markers(:,:,1),nonfiltered_markers(:,:,iMk),gaitAxis,verticalAxis,f);
                % Get the estimated gait events related to the forceplate
                [eFS,eFO] = getEstimatedGaitEvents(btkData,FS_Des,FO_Des,f,side{1,s});
                % !!! the FS and FO are given in frames! some measurements are done with 300Hz and not 150Hz!
                % give all the results in 150Hz : divided by 2 = 150/300
                if isempty(mFS) == 0 && isempty(mFO) == 0 && isempty(eFS) == 0 && isempty(eFO) == 0
                    if f==300
                        FS_stat_dip.Desailly.absError((subj-1)*2+s,iMk) = floor((eFS-mFS)/2);
                        FO_stat_dip.Desailly.absError((subj-1)*2+s,iMk) = floor((eFO-mFO)/2);
                    else
                        FS_stat_dip.Desailly.absError((subj-1)*2+s,iMk) = eFS-mFS;
                        FO_stat_dip.Desailly.absError((subj-1)*2+s,iMk) = eFO-mFO;
                    end
                end
                
                %% Hsue
                [FS_Hsue,FO_Hsue] = Hsue(filtered_markers(:,:,iMk),pelvicMk,gaitAxis,f);
                % Get the estimated gait events related to the forceplate
                [eFS,eFO] = getEstimatedGaitEvents(btkData,FS_Hsue,FO_Hsue,f,side{1,s});
                % !!! the FS and FO are given in frames! some measurements are done with 300Hz and not 150Hz!
                % give all the results in 150Hz : divided by 2 = 150/300
                if isempty(mFS) == 0 && isempty(mFO) == 0 && isempty(eFS) == 0 && isempty(eFO) == 0
                    if f==300
                        FS_stat_dip.Hsue.absError((subj-1)*2+s,iMk) = floor((eFS-mFS)/2);
                        FO_stat_dip.Hsue.absError((subj-1)*2+s,iMk) = floor((eFO-mFO)/2);
                    else
                        FS_stat_dip.Hsue.absError((subj-1)*2+s,iMk) = eFS-mFS;
                        FO_stat_dip.Hsue.absError((subj-1)*2+s,iMk) = eFO-mFO;
                    end
                end

                %% Hreljac
                [FS_Hreljac,FO_Hreljac] = Hreljac(filtered_markers(:,:,iMk),pelvicMk,gaitAxis,verticalAxis,f);
                % Get the estimated gait events related to the forceplate
                [eFS,eFO] = getEstimatedGaitEvents(btkData,FS_Hreljac,FO_Hreljac,f,side{1,s});
                % !!! the FS and FO are given in frames! some measurements are done with 300Hz and not 150Hz!
                % give all the results in 150Hz : divided by 2 = 150/300
                if isempty(mFS) == 0 && isempty(mFO) == 0 && isempty(eFS) == 0 && isempty(eFO) == 0
                    if f==300
                        FS_stat_dip.Hreljac.absError((subj-1)*2+s,iMk) = (eFS-mFS)/2;
                        FO_stat_dip.Hreljac.absError((subj-1)*2+s,iMk) = (eFO-mFO)/2;
                    else
                        FS_stat_dip.Hreljac.absError((subj-1)*2+s,iMk) = eFS-mFS;
                        FO_stat_dip.Hreljac.absError((subj-1)*2+s,iMk) = eFO-mFO;
                    end
                end

            end % Markers
            
            btkDeleteAcquisition(btkData);
        end % side (left,right)
    end
end % subjects

FS_stat_dip.Ghoussayni.Mean = mean(abs(FS_stat_dip.Ghoussayni.absError));
FS_stat_dip.Ghoussayni.Std = std(abs(FS_stat_dip.Ghoussayni.absError));
FO_stat_dip.Ghoussayni.Mean = mean(abs(FO_stat_dip.Ghoussayni.absError));
FO_stat_dip.Ghoussayni.Std = std(abs(FO_stat_dip.Ghoussayni.absError));

FS_stat_dip.Zeni.Mean = mean(abs(FS_stat_dip.Zeni.absError));
FS_stat_dip.Zeni.Std = std(abs(FS_stat_dip.Zeni.absError));
FO_stat_dip.Zeni.Mean = mean(abs(FO_stat_dip.Zeni.absError));
FO_stat_dip.Zeni.Std = std(abs(FO_stat_dip.Zeni.absError));

FS_stat_dip.Desailly.Mean = mean(abs(FS_stat_dip.Desailly.absError));
FS_stat_dip.Desailly.Std = std(abs(FS_stat_dip.Desailly.absError));
FO_stat_dip.Desailly.Mean = mean(abs(FO_stat_dip.Desailly.absError));
FO_stat_dip.Desailly.Std = std(abs(FO_stat_dip.Desailly.absError));

FS_stat_dip.Hsue.Mean = mean(abs(FS_stat_dip.Hsue.absError));
FS_stat_dip.Hsue.Std = std(abs(FS_stat_dip.Hsue.absError));
FO_stat_dip.Hsue.Mean = mean(abs(FO_stat_dip.Hsue.absError));
FO_stat_dip.Hsue.Std = std(abs(FO_stat_dip.Hsue.absError));

FS_stat_dip.Hreljac.Mean = mean(abs(FS_stat_dip.Hreljac.absError));
FS_stat_dip.Hreljac.Std = std(abs(FS_stat_dip.Hreljac.absError));
FO_stat_dip.Hreljac.Mean = mean(abs(FO_stat_dip.Hreljac.absError));
FO_stat_dip.Hreljac.Std = std(abs(FO_stat_dip.Hreljac.absError));

save([exportFolder '\FS_stat_dip.mat'],'FS_stat_dip');
save([exportFolder '\FO_stat_dip.mat'],'FO_stat_dip');

