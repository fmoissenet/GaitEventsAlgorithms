function [mFS,mFO,eFS,eFO] = setMeasuredEstimatedGaitEvents(btkData,FS,FO,n,f)

% Initialisation
mFS = [];
mFO = [];
eFS = [];
eFO = [];

% Compare velocity threshold-based results with forceplate data (up to 4 PFs)
temp = btkGetGroundReactionWrenches(btkData);
for i = 1:length(temp)
    for t = 1:length(temp(1).F(:,3))
        if temp(i).F(t,3) < 15 % 15 N threshold on Z axis (vertical)
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
    twindow = fix(20/150*f); % two consecutive events must be at least distant of 20 frame (at 150 Hz)
    for j = 1:length(FO)
        % Foot strike detection
        if abs(temp1(i)-FS(j)) <= twindow
            mFS = temp1(i); % measured events
            eFS = FS(j); % estimated events
        end
        % Foot off detection
        if abs(temp2(i)-FO(j)) <= twindow
            mFO = temp2(i); % measured events
            eFO = FO(j); % estimated events
        end
    end
end