% define gait events (foot strike and foot off) from Zeni method
% input:
% - marker: 3D coordinates of a foot marker
% - pelvicMk: a struct of the 5 pelvic markers (LASI; RASI; LPSI; RPSI;
% SACR)
% - gaitAxis: number of the gait axis: 1=x, 2=y
% - verticalAxis: number of the vertical axis (3=z)
% - n: number of frames
%
% output:
% - FS = array of foot strikes
% - FO = array of foot offs
function [FS,FO] = Zeni(marker,pelvicMk,gaitAxis)%,verticalAxis,n)
% -------------------------------------------------------------------------
% Initialisation
% -------------------------------------------------------------------------
FS = [];
FO = [];

% -------------------------------------------------------------------------
% Detect events using the relative displacement of the marker to the pelvis
% -------------------------------------------------------------------------
% % relative to the pelvis origin, on the ant/posterior axis of the pelvis
% % (Adapted from SOFAMEAhack2020)
%     for frame=1:n
%         for t=1:3
%             OrPelvis(frame,t) = ((pelvicMk.filtLPSI(frame,t)+pelvicMk.filtRPSI(frame,t))/2+(pelvicMk.filtLASI(frame,t)+pelvicMk.filtRASI(frame,t))/2)/2;
%             yPelvis(frame,t) = (pelvicMk.filtLASI(frame,t)+pelvicMk.filtRASI(frame,t))/2-(pelvicMk.filtLPSI(frame,t)+pelvicMk.filtRPSI(frame,t))/2;
%         end
%     end
%     yPelvis(:,verticalAxis) = zeros(n,1);
%     rel2pelvicOrig = dot(yPelvis,marker-OrPelvis,2);
%     
%     [~,FS] = findpeaks(rel2pelvicOrig);
%     [~,FO] = findpeaks(-rel2pelvicOrig);
                
% relative to the sacrum marker
% (from Zeni paper)
    rel2sacr = marker(:,gaitAxis)-pelvicMk.filtSACR(:,gaitAxis);
    
    [~,FS] = findpeaks(rel2sacr);
    [~,FO] = findpeaks(-rel2sacr);

