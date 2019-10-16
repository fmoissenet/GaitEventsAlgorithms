% path where the data are recorded
Path = 'K:\Ganglabor_Daten\Forschung\UKBB\EventsAutomatisch\Data\';
% hold the description of the data
T = readtable('C:\Users\FreslierM\Documents\Marie\Forschung_Projekte_UKBB\automatisiertenEvents\data_description.xlsx');
left = T.left;
right = T.right;
numbers = T.Number;
clear T
results = struct;
indexResults = 0;
% global variables
fz_Threshold = 10;

for n=1:size(T,1)
    dataStr = numbers{n,1};
    if ~isempty(dataStr)
        % path of these data
%         dataPath = [Path dataStr];
        dataPath = Path;
        indexResults = indexResults + 1;

        %% informations of left side
        infoLeft = left{n,1};
        textSplit = regexp(infoLeft,'Tr \d*| FP \d','match');
        trialCell = regexp(textSplit{1,1},'\d*','match');
        trial = str2num(trialCell{1,1});
        fpCell = regexp(textSplit{1,2},'\d*','match');

        % path of the c3d
        if trial <= 9
%             c3dPath = [dataPath '\' dataStr '_0' trialCell{1,1} '.c3d'];
            c3dPath = [dataPath '\walkNorm_' dataStr 'a0' trialCell{1,1} '.c3d'];
        else
%             c3dPath = [dataPath '\' dataStr '_' trialCell{1,1} '.c3d'];
            c3dPath = [dataPath '\walkNorm_' dataStr 'a' trialCell{1,1} '.c3d'];
        end
        results.names{indexResults,1} = [dataStr '_left'];
        
        % open the c3d
        c3dHandle = btkReadAcquisition(c3dPath);
        
        ff = btkGetFirstFrame(c3dHandle);
        freqPoints = btkGetPointFrequency(c3dHandle);
        % forceplates
        freqAnalog = btkGetAnalogFrequency(c3dHandle);
        forces = btkGetAnalogs(c3dHandle);
        if isfield(forces,['Fz' fpCell{1,1}])
            fz = -1*forces.(['Fz' fpCell{1,1}]);
        else if isfield(forces,['Force_Fz' fpCell{1,1}])
            fz = -1*forces.(['Force_Fz' fpCell{1,1}]);
            else
                disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
                disp(['trial ' dataStr '_' trialCell{1,1} '.c3d has not good force names']);
                disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            end
        end
        
        %% determine FS and FO on this forceplate (Fz>fz_Threshold), into
        % seconds
        [FS_forceplate,TO_forceplate] = getEventsFromFP(fz,fz_Threshold,ff,freqAnalog,freqPoints);
        results.forceplate_FS(indexResults) = FS_forceplate;
        results.forceplate_TO(indexResults) = TO_forceplate;
        
        %% Markers
        markers = btkGetPoints(c3dHandle);
        Left_CoMAll_toSACR = CenterOfMassAllPoints(markers, 'Left', dataStr);
        %% determine FS and TO from algorithm
        left_FS = find_FS_CoM_vel_201811(c3dHandle,Left_CoMAll_toSACR,'Left');
        left_FO = find_FO_CoM_201811(c3dHandle,Left_CoMAll_toSACR,'Left');
        results.kinematics_FS(indexResults) = left_FS;
        results.kinematics_TO(indexResults) = left_FO;
        
        clear Left_CoMAll_toSACR left_FS left_FO
        %% ******************************************************
        indexResults = indexResults + 1;
        % informations of right side
        infoRight = right{n,1};
        textSplit = regexp(infoRight,'Tr \d*| FP \d','match');
        trialCell = regexp(textSplit{1,1},'\d*','match');
        trial = str2num(trialCell{1,1});
        fpCell = regexp(textSplit{1,2},'\d*','match');

        % path of the c3d
        if trial <= 9
            c3dPath = [dataPath '\' dataStr '_0' trialCell{1,1} '.c3d'];
        else
            c3dPath = [dataPath '\' dataStr '_' trialCell{1,1} '.c3d'];
        end
        results.names{indexResults,1} = [dataStr '_right'];
        
        % open the c3d
        c3dHandle = btkReadAcquisition(c3dPath);
        
        ff = btkGetFirstFrame(c3dHandle);
        freqPoints = btkGetPointFrequency(c3dHandle);
        % forceplates
        freqAnalog = btkGetAnalogFrequency(c3dHandle);
        forces = btkGetAnalogs(c3dHandle);
        if isfield(forces,['Fz' fpCell{1,1}])
            fz = -1*forces.(['Fz' fpCell{1,1}]);
        else if isfield(forces,['Force_Fz' fpCell{1,1}])
            fz = -1*forces.(['Force_Fz' fpCell{1,1}]);
            else
                disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
                disp(['trial ' dataStr '_' trialCell{1,1} '.c3d has not good force names']);
                disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            end
        end
        
        %% determine FS and FO on this forceplate (Fz>fz_Threshold), into
        % seconds
        [FS_forceplate,TO_forceplate] = getEventsFromFP(fz,fz_Threshold,ff,freqAnalog,freqPoints);
        results.forceplate_FS(indexResults) = FS_forceplate;
        results.forceplate_TO(indexResults) = TO_forceplate;
        
        %% Markers
        markers = btkGetPoints(c3dHandle);
        right_CoMAll_toSACR = CenterOfMassAllPoints(markers, 'Right', dataStr);
        %% determine FS and TO from algorithm
        right_FS = find_FS_CoM_vel_201811(c3dHandle,right_CoMAll_toSACR,'Right');
        right_FO = find_FO_CoM_201811(c3dHandle,right_CoMAll_toSACR,'Right');
        results.kinematics_FS(indexResults) = right_FS;
        results.kinematics_TO(indexResults) = right_FO;
        
        clear right_CoMAll_toSACR right_FS right_FO
    end
    
end

% diff_FS = (results.kinematics_FS(1,1:20)-results.forceplate_FS(1,1:20))';
% diff_TO = (results.kinematics_TO(1,1:20)-results.forceplate_TO(1,1:20))';
% 
% figure
% plot(diff_FS)
% hold on
% plot(diff_TO,'r')