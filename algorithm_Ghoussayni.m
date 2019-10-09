clearvars;
close all;
clc;
cd('C:\Users\moissene\Documents\Professionnel\publications\articles\1- en cours\Freslier - Gait events');
addpath(genpath('.\btk'));

% Load 3D trajectories of markers and related frame rate and frame number
c3Dfile = 'test.c3d';
btkData = btkReadAcquisition(c3Dfile);
Rmarkers = btkGetMarkers(btkData);
f = btkGetPointFrequency(btkData);
n = btkGetPointFrameNumber(btkData);

% Filter 3D trajectories of markers
% Butterworth filter, zero-phase filter, 2nd order, low-pass, cut-off 10 Hz
rightFootMarkersName = {'R_FM5','R_FCC'};
Fmarkers = [];
for i = 1:length(rightFootMarkersName)
    [B,A] = butter(2,10/(f/2),'low');
    Fmarkers.(rightFootMarkersName{i}) = filtfilt(B, A, Rmarkers.(rightFootMarkersName{i}));
end

% Calculate the velocities of the markers in sagittal plane (along X axis)
for t = 1:n-1
    for i = 1:length(rightFootMarkersName)
        Vmarkers.(rightFootMarkersName{i})(t) = (Fmarkers.(rightFootMarkersName{i})(t+1,1)- ...
                                        Fmarkers.(rightFootMarkersName{i})(t,1))/ ...
                                        (1/f);
    end
end

% Determine the average and standard deviation of sagittal velocities for
% each marker during contact period by visual inspection (from 2 subjects,
% and for each walking speed)
startContact_S1_V1.(rightFootMarkersName{1}) = 24; % FCC from test.c3d
stopContact_S1_V1.(rightFootMarkersName{1}) = 153; % FCC from test.c3d
averageVelocity_S1_V1.(rightFootMarkersName{1}) = ...
    mean(Vmarkers.R_FCC(startContact_S1_V1.(rightFootMarkersName{1}):...
                        stopContact_S1_V1.(rightFootMarkersName{1})));
stdVelocity_S1_V1.(rightFootMarkersName{1}) = ...
    std(Vmarkers.R_FCC(startContact_S1_V1.(rightFootMarkersName{1}):...
                        stopContact_S1_V1.(rightFootMarkersName{1})));
% --
startContact_S1_V1.(rightFootMarkersName{2}) = 11; % FM5 from test.c3d
stopContact_S1_V1.(rightFootMarkersName{2}) = 200; % FM5 from test.c3d
averageVelocity_S1_V1.(rightFootMarkersName{2}) = ...
    mean(Vmarkers.R_FCC(startContact_S1_V1.(rightFootMarkersName{2}):...
                        stopContact_S1_V1.(rightFootMarkersName{2})));
stdVelocity_S1_V1.(rightFootMarkersName{2}) = ...
    std(Vmarkers.R_FCC(startContact_S1_V1.(rightFootMarkersName{2}):...
                        stopContact_S1_V1.(rightFootMarkersName{2})));
                          
% Velocity thresholds (empirically set)
% 50 mm/s in the original article for barefoot gait
potentialThresholds = [];
for i = 1:length(rightFootMarkersName)
    potentialThresholds = [potentialThresholds ...
        averageVelocity_S1_V1.(rightFootMarkersName{i}) - ...
        2*stdVelocity_S1_V1.(rightFootMarkersName{i})];
end
velocityThreshold = max(potentialThresholds);

% Detect events using the velocity threshold
for i = 1:length(rightFootMarkersName)
    Tevents(i).FS = [];
    Tevents(i).FO = [];
    for t = 1:n-1
        if Vmarkers.(rightFootMarkersName{i})(t) >= velocityThreshold && ...
           isempty(Tevents(i).FS)
            Tevents(i).FS(1) = t;
        end
        if ~isempty(Tevents(i).FS)
            if t > Tevents(i).FS(1)+n*0.1 && ... % move 10% of gait cycle from 1st FS
               isempty(Tevents(i).FO) && ...
               Vmarkers.(rightFootMarkersName{i})(t) < velocityThreshold
                Tevents(i).FO(1) = t;
            end
        end
        if ~isempty(Tevents(i).FO) && length(Tevents(i).FS) <= 1
            if t > Tevents(i).FO(1)+n*0.1 && ... % move 10% of gait cycle from 1st FO
                Vmarkers.(rightFootMarkersName{i})(t) >= velocityThreshold
                Tevents(i).FS(2) = t;
            end
        end
    end
end
tFS1 = [];
tFS2 = [];
tFO1 = [];
for i = 1:length(rightFootMarkersName)
    tFS1 = [tFS1 Tevents(i).FS(1)];    
    tFS2 = [tFS2 Tevents(i).FS(2)];
    tFO1 = [tFO1 Tevents(i).FO(1)];
end
Revents.FS(1) = min(tFS1);
Revents.FS(2) = min(tFS2);
Revents.FO(1) = max(tFO1);

% Plot results
figure;
hold on;
ylim([-3000 500]);
for i = 1:length(rightFootMarkersName)
    plot(Vmarkers.(rightFootMarkersName{i}));
end
line([Revents.FS(1) Revents.FS(1)],[-1e6 1e6],'Color','red','LineStyle','-');
line([Revents.FS(2) Revents.FS(2)],[-1e6 1e6],'Color','red','LineStyle','-');
line([Revents.FO(1) Revents.FO(1)],[-1e6 1e6],'Color','red','LineStyle','--');

% Compare velocity threshold-based results with forceplate data
temp = btkGetGroundReactionWrenches(btkData);
for t = 1:length(temp(1).F(:,3))
    if temp(1).F(t,3) < 10 % 10 N threshold on Z axis (vertical)
        temp(1).F(t,:) = 0;
    end
    if temp(2).F(t,3) < 10 % 10 N threshold on Z axis (vertical)
        temp(2).F(t,:) = 0;
    end
end
Forces = interpft([temp(1).F(:,3) temp(2).F(:,3)],n); % Z axis (vertical)
temp11 = find(Forces(:,1)>1e-4);
temp21 = find(Forces(:,2)>1e-4);
temp12 = find(Forces(temp11(1):end,1)<1e-4)+temp11(1);
temp22 = find(Forces(temp21(1):end,2)<1e-4)+temp21(1);
if abs(temp11(1)-Revents.FS(1)) < abs(temp21(1)-Revents.FS(1))
    eFS = temp11(1)-Revents.FS(1);
    eFO = temp12(1)-Revents.FO(1);
    line([temp11(1) temp11(1)],[-1e6 1e6],'Color','green','LineStyle','-');
    line([temp12(1) temp12(1)],[-1e6 1e6],'Color','green','LineStyle','-');
else
    eFS = temp21(1)-Revents.FS(1);
    eFO = temp22(1)-Revents.FO(1);
    line([temp21(1) temp21(1)],[-1e6 1e6],'Color','green','LineStyle','-');
    line([temp22(1) temp22(1)],[-1e6 1e6],'Color','green','LineStyle','-');
end