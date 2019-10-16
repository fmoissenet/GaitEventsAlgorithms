function [FS,TO] = getEventsFromFP(Fz,threshold,firstframe,freqAnalog,freqKinematics)
% function to determine the footstrike and the footoff from a given
% forceplate and a given force threshold. It is hypothesized that the
% forceplate is only one time hit
% 
% Fz = data of the vertical force of the forceplate in N (positive value)
% threshold = given threshold (N) for determination of the events
% firstframe = first frame (s) in the trial
% freqAnalog = frequency of the forceplate data
% 
% FS = determined foot strike (s)
% TO = determined foot off (s)

    indexAnalog = 1;
    FS = 0;
    TO = 0;
    while (indexAnalog < length(Fz)) && FS == 0
        if Fz(indexAnalog)>threshold
            % test if the force is > threshold during 5 kinematics frames
            t = 1;
            FS_OK = 1;
            while FS_OK && t<(5*freqAnalog/freqKinematics)
                if Fz(indexAnalog+t)>threshold
                    FS_OK = 1;
                else
                    FS_OK = 0;
                end
                t = t + 1;
            end
            if FS_OK
                FS = (indexAnalog - 1 )/freqAnalog + (firstframe-1)/freqKinematics;
            end
        end
        indexAnalog = indexAnalog + 1;
    end
    indexAnalog = indexAnalog + 5*freqAnalog/freqKinematics; %(5 frames from kinematics)
    while (indexAnalog < length(Fz)) && TO == 0
        if Fz(indexAnalog)<threshold
            TO = (indexAnalog - 1 )/freqAnalog + (firstframe-1)/freqKinematics;
        end
        indexAnalog = indexAnalog + 1;
    end

end