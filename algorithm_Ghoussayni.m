clearvars;
close all;
clc;
cd('C:\Users\moissene\Documents\Professionnel\publications\articles\1- en cours\Freslier - Gait events');
addpath(genpath('.\btk'));

% Load 3D trajectories of markers and related frame rate and frame number
c3Dfile = './C3D Files/EVT127.c3d';
btkData = btkReadAcquisition(c3Dfile);
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

% Define gait axis and direction
tdiff = Fmarkers.(rightFootMarkersName{i})(end,:) - Fmarkers.(rightFootMarkersName{i})(1,:);
[~,gaitAxis] = max(abs(tdiff));
verticalAxis = 3; % Z axis vertical
if tdiff(gaitAxis) > 0
    gaitDirection = 1;
else
    gaitDirection = -1;
end

% Calculate the velocities of the markers along gait (V1) and vertical (V2) axes
for t = 1:n-1
    for i = 1:length(rightFootMarkersName)
        V1markers.(rightFootMarkersName{i})(t) = gaitDirection*...
                                                 (Fmarkers.(rightFootMarkersName{i})(t+1,gaitAxis)- ...
                                                  Fmarkers.(rightFootMarkersName{i})(t,gaitAxis))/ ...
                                                  (1/f);
        V2markers.(rightFootMarkersName{i})(t) = abs(Fmarkers.(rightFootMarkersName{i})(t+1,verticalAxis)- ...
                                                  Fmarkers.(rightFootMarkersName{i})(t,verticalAxis))/ ...
                                                  (1/f);
    end
end
                         
% Velocity thresholds (empirically set)
% 50 mm/s in the original article for barefoot gait
% 500 mm/s in Bruening et al., 2014
% Adjusted for vertical velocity to 200 mm/s
v1Threshold = 1000;
v2Threshold = 200;

% Detect events using the velocity threshold
Tevents.R_FS = [];
Tevents.R_FO = [];
for t = 1:n-1
    % Foot strike defined using heel marker
    if isempty(Tevents.R_FS) && isempty(Tevents.R_FO)
        if (V1markers.RHEE(t) <= v1Threshold && ...
            V2markers.RHEE(t) <= v2Threshold)
                Tevents.R_FS = [Tevents.R_FS t];
        end
    elseif ~isempty(Tevents.R_FS) && isempty(Tevents.R_FO)
    	% Do nothing: wait for a first FO (assume that we detect first a FS)
    elseif ~isempty(Tevents.R_FS) && ~isempty(Tevents.R_FO) && ...
           length(Tevents.R_FS) > length(Tevents.R_FO)
    	% Do nothing: wait for the next FO (assume that we detect first a FS)
    elseif ~isempty(Tevents.R_FS) && ~isempty(Tevents.R_FO) && ...
           length(Tevents.R_FS) == length(Tevents.R_FO)
        if (V1markers.RHEE(t) <= v1Threshold && ...
            V2markers.RHEE(t) <= v2Threshold)
                Tevents.R_FS = [Tevents.R_FS t];
        end
    end
    % Foot off defined using forefoot marker
    if isempty(Tevents.R_FS) && isempty(Tevents.R_FO)
        % Do nothing: wait for a first FS (assume that we detect first a FS)
    elseif ~isempty(Tevents.R_FS) && isempty(Tevents.R_FO)
        if (V1markers.RDMT5(t) >= v1Threshold && ...
            V2markers.RDMT5(t) >= v2Threshold)
                Tevents.R_FO = [Tevents.R_FO t];
        end
    elseif ~isempty(Tevents.R_FS) && ~isempty(Tevents.R_FO) && ...
           length(Tevents.R_FO) < length(Tevents.R_FS)
        if (V1markers.RDMT5(t) >= v1Threshold && ...
            V2markers.RDMT5(t) >= v2Threshold)
                Tevents.R_FO = [Tevents.R_FO t];
        end
    elseif ~isempty(Tevents.R_FS) && ~isempty(Tevents.R_FO) && ...
           length(Tevents.R_FO) == length(Tevents.R_FS)
        % Do nothing: wait for the nest FS (assume that we detect first a FS)
    end
end

% Compare velocity threshold-based results with forceplate data (up to 4 PFs)
temp = btkGetGroundReactionWrenches(btkData);
for i = 1:length(temp)
    for t = 1:length(temp(1).F(:,3))
        if temp(i).F(t,3) < 10 % 10 N threshold on Z axis (vertical)
            temp(i).F(t,:) = 0;
        end
    end
    Forces(:,i) = interpft(temp(i).F(:,3),n); % Z axis (vertical)
    temp1(i) = NaN;
    if min(find(Forces(:,i)>1e-4))
        temp1(i) = min(find(Forces(:,i)>1e-4));
    end
    if ~isnan(temp1(i))
        temp2(i) = min(find(Forces(temp1(i):end,i)<1e-4)+temp1(i));
    else
        temp2(i) = NaN;
    end
    twindow = 10; % find similar event in a 10 frame window
    for j = 1:length(Tevents.R_FO)
        % Foot strike detection
        if abs(temp1(i)-Tevents.R_FS(j)) <= twindow
            mFS = temp1(i); % measured events
            eFS = Tevents.R_FS(j); % estimated events
        end
        % Foot off detection
        if abs(temp2(i)-Tevents.R_FO(j)) <= twindow
            mFO = temp2(i); % measured events
            eFO = Tevents.R_FO(j); % estimated events
        end
    end
end

% Plot results
figure;
hold on;
ylim([0 4000]);
plot(V1markers.RHEE,'Color','blue');
plot(V2markers.RHEE,'Color','green');
plot(V1markers.RDMT5,'Color','blue','LineStyle','--');
plot(V2markers.RDMT5,'Color','green','LineStyle','--');
% for i = 1:length(Tevents.R_FO)
%     line([Tevents.R_FS(i) Tevents.R_FS(i)],[-1e6 1e6],'Color','red','LineStyle','-');
%     line([Tevents.R_FO(i) Tevents.R_FO(i)],[-1e6 1e6],'Color','red','LineStyle','--');
% end
line([mFS mFS],[-1e6 1e6],'Color','red','LineStyle','-');
line([eFS eFS],[-1e6 1e6],'Color','blue','LineStyle','-');
line([mFO mFO],[-1e6 1e6],'Color','red','LineStyle','--');
line([eFO eFO],[-1e6 1e6],'Color','blue','LineStyle','--');

% Store errors (frames)
errorFS = eFS - mFS
errorFO = eFO - mFO