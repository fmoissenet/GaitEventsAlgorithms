function FO_sec = find_FO_CoM_201811( file_handle, ref_data, side )
%FIND_FO_COM function to calculate Foot Offs from local first max values of
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
    
    FO = [];
    FO_sec = [];
    peaks_values = [];

    freq = btkGetPointFrequency(file_handle);
    ff = btkGetFirstFrame(file_handle);

    if ~isfield(events,[side '_Event'])
        disp('!!!!!!!!!');
        disp('No General Events found');
        disp('!!!!!!!!!');
        return
    else
        % events is in seconds, transform in frames (0s ~ first frame =>
        % +1)
        begin_frame = round(events.([side '_Event'])(1) * freq) -(ff-1) + 1;
    end
    
    [peaks,ind_peaks] = findpeaks(-ref_data(begin_frame:length(ref_data),1));
    if ~isempty(ind_peaks)
        FO(1) = ind_peaks(1) + (begin_frame - 1);
        FO_sec(1) = ((FO(1)+(ff-1)) - 1) / freq; % time in second
%         events_new = btkAppendEvent(file_handle, 'Foot Off', FO_sec(1), side);
        peaks_values(1) = peaks(1);
%         if length(ind_peaks)>1
%             FO(2) = ind_peaks(2) + (begin_frame - 1);
%             FO_sec(2) = ((FO(2)+(ff-1)) - 1) / freq; % time in second
%             events_new = btkAppendEvent(file_handle, 'Foot Off', FO_sec(2), side);
%             peaks_values(2) = peaks(2);
%         else
%             disp('!!!!!!!!!');
%             disp(['Second ' side ' Foot Off not found.']);
%             disp('!!!!!!!!!');
%         end
    else
        disp('!!!!!!!!!');
        disp(['First ' side ' Foot Off not found.']);
        disp('!!!!!!!!!');
    end
        
end

