function [FS,FO] = algorithm_Ghoussayni(Fmarkers,rightFootMarkersName,gaitAxis,verticalAxis,n,f)

% Calculate the 2D velocity of the markers in the plane containing 
% gait (V1) and vertical (V2) axes
for t = 1:n-1
    for i = 1:length(rightFootMarkersName)
        Vmarkers.(rightFootMarkersName{i})(t) = sqrt((Fmarkers.(rightFootMarkersName{i})(t+1,gaitAxis)- ...
                                                      Fmarkers.(rightFootMarkersName{i})(t,gaitAxis))^2+ ...
                                                     (Fmarkers.(rightFootMarkersName{i})(t+1,verticalAxis)- ...
                                                      Fmarkers.(rightFootMarkersName{i})(t,verticalAxis))^2)/ ...
                                                (1/f);
    end
end
                         
% Velocity threshold (empirically set)
% 50 mm/s in the original article for barefoot gait
% 500 mm/s in Bruening et al., 2014
vThreshold = 500;

% Detect events using the velocity threshold
FS = [];
FO = [];
twindow = fix(30/150*f); % two consecutive events must be at least distant of 30 frame (at 150 Hz)
for t = 1:n-1
    % Foot strike defined using heel marker
    if isempty(FS) && isempty(FO)
        if Vmarkers.RHEE(t) <= vThreshold
                FS = [FS t];
        end
    elseif ~isempty(FS) && isempty(FO)
    	% Do nothing: wait for a first FO (assume that we detect first a FS)
    elseif ~isempty(FS) && ~isempty(FO) && ...
           length(FS) > length(FO)
    	% Do nothing: wait for the next FO (assume that we detect first a FS)
    elseif ~isempty(FS) && ~isempty(FO) && ...
           length(FS) == length(FO)
        if Vmarkers.RHEE(t) <= vThreshold && ...
           t >= FO(end)+twindow
                FS = [FS t];
        end
    end
    % Foot off defined using forefoot marker
    if isempty(FS) && isempty(FO)
        % Do nothing: wait for a first FS (assume that we detect first a FS)
    elseif ~isempty(FS) && isempty(FO)
        if Vmarkers.RDMT5(t) >= vThreshold
                FO = [FO t];
        end
    elseif ~isempty(FS) && ~isempty(FO) && ...
           length(FO) < length(FS)
        if Vmarkers.RDMT5(t) >= vThreshold && ...
           t >= FS(end)+twindow
                FO = [FO t];
        end
    elseif ~isempty(FS) && ~isempty(FO) && ...
           length(FO) == length(FS)
        % Do nothing: wait for the nest FS (assume that we detect first a FS)
    end
end