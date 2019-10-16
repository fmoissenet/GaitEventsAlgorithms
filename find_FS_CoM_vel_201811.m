function FS_sec = find_FS_CoM_vel_201811(file_handle, ref_data, side )
%find_FS_CoM_vel_201811 function to calculate Foot Strikes from local first max values of
%ref_data
%   file_handle = handle from btkReadAcquisition
%   ref_data = data where the events will be searched (CoM of several
%   Markers) (first column: gait direction, second column: vertical)
%   side = 'Left' or 'Right'

    % sometimes the patient ID is in the events name, therefore it has to be corrected
    for t = 1:btkGetEventNumber(file_handle)
        btkSetEventSubject(file_handle,t,'');
    end     
    events = btkGetEvents(file_handle); % all events in seconds   
    
    FS = [];
    FS_sec = [];
%     peaks_values = [];

    freq = btkGetPointFrequency(file_handle);
    ff = btkGetFirstFrame(file_handle);
    
    if ~isfield(events,[side '_Event'])
        disp('!!!!!!!!!');
        disp('No General Events found');
        disp('!!!!!!!!!');
        return
    else
        % events is in seconds, transform in frames (0s ~ frame '1' =>
        % +1)
        begin_frame = round(events.([side '_Event'])(1) * freq) - (ff-1) + 1;
    end
    
    % peaks in the difference CoM to SACR from given begin in Trial
    [~,ind_peaks] = findpeaks(ref_data(begin_frame:length(ref_data),1));
    % Vertical Velocity of Ref_data (CoM)
    vel_ref_data = diff(ref_data(:,2));
    
    if ~isempty(ind_peaks)
        % find the nearest smallest peak of velocity of vertical
        % displacement of CoM to peak of horizontal displacement (in a
        % window +/- kinematics frequency/5 frames around the peak (=60 frames
        % for 300Hz, 30 frames for 150Hz)
        windowRange = freq/5;
        % if the peaks found is near the begin of the frame
        lowRange = max(1,ind_peaks(1)+ (begin_frame - 1) -windowRange);
        
        [peaks_vel,ind_peaks_vel_1] = findpeaks(-vel_ref_data(lowRange:ind_peaks(1)+ (begin_frame - 1) +windowRange));
        if ~isempty(ind_peaks_vel_1)
            % if several peaks, take the nearest to ind_peaks(1) one
            [~,peak_nb_1] = min(abs(ind_peaks(1)+(begin_frame-1) - (ind_peaks_vel_1+lowRange-1)));
            % in frames from c3d Trial
            FS(1) = ind_peaks_vel_1(peak_nb_1) + lowRange -1;
            % time in second
            FS_sec(1) = ((FS(1)+(ff-1)) - 1) / freq;
%             events_new = btkAppendEvent(file_handle, 'Foot Strike', FS_sec(1), side);
%             peaks_values(1) = peaks_vel(peak_nb_1);
%             if length(ind_peaks)>1
%                 [peaks_vel,ind_peaks_vel_2] = findpeaks(-vel_ref_data(ind_peaks(2)+ (begin_frame - 1) -30:ind_peaks(2)+ (begin_frame - 1) +30));
%                 if ~isempty(ind_peaks_vel_2)
%                     % if several peaks, take the nearest to ind_peaks(1) one
%                     [~,peak_nb_2] = min(abs(31-ind_peaks_vel_2));
%                     % in frames from c3d Trial
%                     FS(2) = ind_peaks_vel_2(peak_nb_2) + (ind_peaks(2)-30) -1 + (begin_frame - 1);
%                     % time in second
%                     FS_sec(2) = ((FS(2)+(ff-1)) - 1) / freq;
%                     events_new = btkAppendEvent(file_handle, 'Foot Strike', FS_sec(2), side);
%                     peaks_values(2) = peaks_vel(peak_nb_2);
%                 else
%                     disp('!!!!!!!!!');
%                     disp(['Second ' side ' Foot Strike not found.']);
%                     disp('!!!!!!!!!');
%                 end
%             else
%                 disp('!!!!!!!!!');
%                 disp(['Second ' side ' Foot Strike not found.']);
%                 disp('!!!!!!!!!');
%             end
        else
            disp('!!!!!!!!!');
            disp(['First ' side ' Foot Strike not found.']);
            disp('!!!!!!!!!');
        end
    else
        disp('!!!!!!!!!');
        disp(['First ' side ' Foot Strike not found.']);
        disp('!!!!!!!!!');
    end

    
end

