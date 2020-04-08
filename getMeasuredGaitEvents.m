% determine the measured events from a given c3d
% input:
% - btkData: btk handle of the c3d file
% - n: number of frames
% - FP_number: forceplate number where the measured events are determined
%
% output: (=NaN if not found)
% - mFS = measured foot strike from forceplate
% - mFO = measured foot off from forceplate
function [mFS,mFO] = getMeasuredGaitEvents(btkData,n,FP_number)

% frequencies
freqKinematics = btkGetPointFrequency(btkData);
freqAnalog = btkGetAnalogFrequency(btkData);

% Compare determined results (from  some algorithm) with forceplate data (up to 4 FPs)
GRF = btkGetGroundReactionWrenches(btkData);
Force_z = interpft(GRF(FP_number).F(:,3),n); % Z axis (vertical)
clear GRF

% FS from forceplate
mFS = NaN;
% % if some artefact in the forceplate before the foot strike, it doesn't
% % work...
% if find(Force_z >= 15,1) % first index with Fz>=15 N threshold on Z axis (vertical)
%     mFS = find(Force_z >= 15,1);
% end
% check for Force_z>threshold during 5 kinematics frames
indexAnalog = 1;
while (indexAnalog < length(Force_z)) && isnan(mFS)
    if Force_z(indexAnalog)>=15
        % test if the force is >= threshold during 5 kinematics frames
        t = 1;
        FS_OK = 1;
        while FS_OK && t<(5*freqAnalog/freqKinematics)
            if Force_z(indexAnalog+t)>=15
                FS_OK = 1;
            else
                FS_OK = 0;
            end
            t = t + 1;
        end
        if FS_OK
            mFS = indexAnalog;
        end
    end
    indexAnalog = indexAnalog + 1;
end


% FO from forceplate
mFO = NaN;
if ~isnan(mFS)
    mFO = find(Force_z(mFS:end)<15,1)+mFS-1; % first index with Fz<15 N threshold on Z axis (vertical)
end
