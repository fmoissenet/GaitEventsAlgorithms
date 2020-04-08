% define gait events (foot strike and foot off) from Ghoussayni method
% threshold on 500 mm/s (Bruening 2014)
% input:
% - marker: 3D coordinates of a foot marker
% - gaitAxis: number of the gait axis: 1=x, 2=y
% - verticalAxis: number of the vertical axis (3=z)
% - n: number of frames
% - f: frequency of the data (markers)
%
% output:
% - FS = array of foot strikes
% - FO = array of foot offs
function [FS,FO] = Ghoussayni(marker,gaitAxis,verticalAxis,n,f)
% -------------------------------------------------------------------------
% Initialisation
% -------------------------------------------------------------------------
FS = [];
FO = [];

% -------------------------------------------------------------------------
% Calculate the 2D velocity of the marker in the plane containing
% gait (V1) and vertical (V2) axes
% -------------------------------------------------------------------------
% marker velocity
for t = 1:n-1
    velocity(t) = sqrt((marker(t+1,gaitAxis)- ...
        marker(t,gaitAxis))^2+ ...
        (marker(t+1,verticalAxis)- ...
        marker(t,verticalAxis))^2)/ ...
        (1/f);
end
% plot(velocity)
% line([1 length(velocity)],[500 500])

% -------------------------------------------------------------------------
% Velocity threshold (empirically set)
% 50 mm/s in the original article for barefoot gait
% 500 mm/s in Bruening et al., 2014
% -------------------------------------------------------------------------
vThreshold = 500;

% -------------------------------------------------------------------------
% Detect events using the velocity threshold
%   The event is defined when the marker has a velocity under
%   threshold for FS, threshold for FO
% -------------------------------------------------------------------------
twindow = fix(30/150*f); % two consecutive events must be at least distant of 30 frame (at 150 Hz)
for t = 1:n-1
    % Foot strike defined using given marker
    if isempty(FS) && isempty(FO)
        temp = [];
            if velocity(t) <= vThreshold
                temp = t;
            end
        if ~isempty(temp)
            FS = [FS temp];
        end
    elseif ~isempty(FS) && isempty(FO)
        % Do nothing: wait for a first FO (assume that we detect first a FS)
    elseif ~isempty(FS) && ~isempty(FO) && ...
            length(FS) > length(FO)
        % Do nothing: wait for the next FO (assume that we detect first a FS)
    elseif ~isempty(FS) && ~isempty(FO) && ...
            length(FS) == length(FO)
        temp = [];
            if velocity(t) <= vThreshold && ...
               t >= FO(end)+twindow
                temp = t;
            end
        if ~isempty(temp)
            FS = [FS temp];
        end
    end
    % Foot off defined using given marker
    if isempty(FS) && isempty(FO)
        % Do nothing: wait for a first FS (assume that we detect first a FS)
    elseif ~isempty(FS) && isempty(FO)
        temp = [];
            if velocity(t) >= vThreshold && ...
               t >= FS(end)+twindow
                temp = t;
            end
        if ~isempty(temp)
            FO = [FO temp];
        end
    elseif ~isempty(FS) && ~isempty(FO) && ...
            length(FO) < length(FS)
        temp = [];
            if velocity(t) >= vThreshold && ...
               t >= FS(end)+twindow
                temp = t;
            end
        if ~isempty(temp)
            FO = [FO temp];
        end
    elseif ~isempty(FS) && ~isempty(FO) && ...
            length(FO) == length(FS)
        % Do nothing: wait for the next FS (assume that we detect first a FS)
    end
end