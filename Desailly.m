% define gait events (foot strike and foot off) from Desailly method
% input:
% - heelMk: non filtered 3D coordinates of heel marker (for gait frequency)
% - footMk: non filtered 3D coordinates of a foot marker
% - gaitAxis: number of the gait axis: 1=x, 2=y
% - verticalAxis: number of the vertical axis (3=z)
% - f: frequence
%
% output:
% - FS = array of foot strikes
% - FO = array of foot offs
function [FS,FO] = Desailly(heelMk,footMk,gaitAxis,verticalAxis,f)
% -------------------------------------------------------------------------
    % 4th order butterworth lowpass filter (7 Hz)
% -------------------------------------------------------------------------
    [B,A] = butter(4,(7/(f/2)));
    filtered_marker = filtfilt(B, A, footMk);
    filt_heelMk = filtfilt(B, A, heelMk);

% -------------------------------------------------------------------------
    % determine the gait frequency from the vertical component of the heel marker
    % do not work with the other markers: displacements have too many peaks
% -------------------------------------------------------------------------
    [~,frame_pks]=findpeaks(filt_heelMk(:,verticalAxis),'MinPeakHeight',100); %25 is too low for the norm data
    if length(frame_pks)<2 % not enough peaks: try with the minimum peaks
        [~,frame_pks]=findpeaks(-filt_heelMk(:,verticalAxis),'MinPeakHeight',-100);
    end
    stridetime=(frame_pks(2)-frame_pks(1))/f;
    gaitFreq = 1/stridetime;
    
% -------------------------------------------------------------------------
    % highpass filter (0.5*gait frequency) of the horizontal component of the marker
% -------------------------------------------------------------------------
    [z,p,k] = butter(4,((0.5*gaitFreq)/(f/2)),'high');
    [sos,g] = zp2sos(z,p,k);
    highPass_mk_FS  = filtfilt(sos,g,filtered_marker(:,gaitAxis));
    [~,FS] = findpeaks(highPass_mk_FS);
% plot(highPass_mk_FS)
% % -------------------------------------------------------------------------
%     % original: highpass filter (1.1*gait frequency) of the horizontal component of the marker
% % -------------------------------------------------------------------------
%     [z2,p2,k2] = butter(4,((1.1*gaitFreq)/(f/2)),'high');
%     [sos2,g2] = zp2sos(z2,p2,k2);
%     highPass_mk_FO  = filtfilt(sos2,g2,filtered_marker(:,gaitAxis));
%     [~,FO] = findpeaks(-highPass_mk_FO);
    
% -------------------------------------------------------------------------
    % Bruening & Goncalves proposed 0.5*gait frequency for the FO:
% -------------------------------------------------------------------------
    [~,FO] = findpeaks(-highPass_mk_FS);

