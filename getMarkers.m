% function to extract from the c3d data the listed markers, to rotate them
% so that the x-axis (1st axis) is the gait direction, to filtered them
% Input:
% - btkData = btk handle of the c3d file
% - f: frequency of the data (markers)
% - side: leg side ('L' or 'R')
% - MarkersName = list of the markers, without the side information
%
% Output:
% - nonfiltered_markers = rotated markers non filtered
% - filtered_markers = rotated markers filtered
% - pelvicMk = all pelvic markers, filtered
function [nonfiltered_markers,filtered_markers,pelvicMk] = ...
            getMarkers(btkData,f,side,MarkersName)
    
    Markers = btkGetMarkers(btkData);
    
    %% rotate the markers to the gait axis
    markers_corrected = f_rotCoordinateSystem(Markers);

    %% Filter 3D trajectories of markers
    % Butterworth filter, zero-phase filter, 2nd order, low-pass, cut-off 10 Hz
    [B,A] = butter(2,10/(f/2),'low');
    nonfiltered_markers = [];
    filtered_markers = [];
    for i = 1:length(MarkersName)
        nonfiltered_markers(:,:,i) = markers_corrected.([side MarkersName{1,i}]);
        filtered_markers(:,:,i) = filtfilt(B, A, markers_corrected.([side MarkersName{1,i}]));
    end
    
    pelvicMk.filtLASI = filtfilt(B, A, markers_corrected.LASI);
    pelvicMk.filtRASI = filtfilt(B, A, markers_corrected.RASI);
    pelvicMk.filtLPSI = filtfilt(B, A, markers_corrected.LPSI);
    pelvicMk.filtRPSI = filtfilt(B, A, markers_corrected.RPSI);
    pelvicMk.filtSACR = filtfilt(B, A, markers_corrected.SACR);