% determine the estimated events from a given c3d
% input:
% - btkData: btk handle of the c3d file
% - FS: list of the estimated foot strikes
% - FO: list of the estimated foot offs
% - f: frequency of the data (markers)
% - side: leg side of the subject where the events are determined
%
% output: (=NaN if not found)
% - eFS = estimated foot strike (first estimated event from the given list
% and since a general event set in the c3d)
% - eFO = estimated foot off (first estimated event from the given list
% and since a general event set in the c3d)
function [eFS,eFO] = getEstimatedGaitEvents(btkData,FS,FO,f,side)

% Initialisation
eFS = [];
eFO = [];

% get the general event: FO and FS are the first one found from this event
events = btkGetEvents(btkData);
ff = btkGetFirstFrame(btkData);
switch side
    case 'L'
        genEvent = events.Left_Event; % in seconds
    case 'R'
        genEvent = events.Right_Event;
end
genEvent_frames = (genEvent*f + 1) - ff + 1;

% get the first FO and FS from the general event
i=1;
while isempty(eFS) && i<=length(FS)
    if FS(i) > genEvent_frames
        eFS = FS(i);
        j = 1;
        while isempty(eFO) && j<=length(FO)
            if FO(j) > eFS
                eFO = FO(j);
            else
                j = j+1;
            end
        end
    else
        i = i+1;
    end
end
if isempty(eFS)
    eFS = NaN;
end
if isempty(eFO)
    eFO = NaN;
end
