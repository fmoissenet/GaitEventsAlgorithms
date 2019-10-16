function CoMAll_toSACR = CenterOfMassAllPoints(markers, side, dataStr)
% function to calculate the movement of a center of mass of points to the
% sacrum marker
% markers = struct of markers from the c3d (btk)
% side = 'Left' or 'Right'
% dataStr = string of the subject name
    if isfield(markers,'SACR')
        prefix = '';
    else
        prefix = [dataStr '_'];
    end
    switch side
        case 'Left'
            firstLetter = 'L';
        case 'Right'
            firstLetter = 'R';
    end
    
    SACR = markers.([prefix 'SACR']);
    % Hindfoot
    HEE = markers.([prefix firstLetter 'HEE']);
    ANK = markers.([prefix firstLetter 'ANK']);
    TPR = markers.([prefix firstLetter 'TPR']);
    SITA = markers.([prefix firstLetter 'SITA']);
    % Forefoot
    PMT1 = markers.([prefix firstLetter 'PMT1']);
    DMT1 = markers.([prefix firstLetter 'DMT1']);
    PMT5 = markers.([prefix firstLetter 'PMT5']);
    DMT5 = markers.([prefix firstLetter 'DMT5']);
    CUN = markers.([prefix firstLetter 'CUN']);
    TOE = markers.([prefix firstLetter 'TOE']);
    HLX = markers.([prefix firstLetter 'HLX']);
    
    %% find walking direction
    SACR_proper = SACR;
    % delete zeros at the beginning or end of an trial
    SACR_proper(sum(SACR==0,2)>0,:) = [];
    dir_i = abs(SACR_proper(end, 1) - SACR_proper(1, 1)); 
    dir_j = abs(SACR_proper(end, 2) - SACR_proper(1, 2)); 
    walkdir = 1;  % x is walkdir
    if (dir_i < dir_j)  
        walkdir = 2;  % y is walkdir
    end;
    % pos. or neg. direktion on axis
    sgn = sign(SACR(end, walkdir) - SACR(1, walkdir));
    walkdir = walkdir * sgn;
    clear SACR_proper dir_i dir_j sgn

    %% define group of markers
    % CoG of all markers (1st column: gait direction, 2nd column:
    % vertical deplacement)
    CoMAll_toSACR = zeros(length(SACR),2);
    for i=1:length(SACR)
        % CoG of all markers
        CoMAll_toSACR(i,1) = sign(walkdir) * ( mean([HEE(i,abs(walkdir)); ANK(i,abs(walkdir)); TOE(i,abs(walkdir));...
                                                TPR(i,abs(walkdir)); SITA(i,abs(walkdir)); PMT1(i,abs(walkdir));...
                                                DMT1(i,abs(walkdir)); PMT5(i,abs(walkdir)); DMT5(i,abs(walkdir));...
                                                CUN(i,abs(walkdir)); HLX(i,abs(walkdir))]) - SACR(i,abs(walkdir)) );
        CoMAll_toSACR(i,2) = mean([HEE(i,3); ANK(i,3); TOE(i,3);...
                                TPR(i,3); SITA(i,3); PMT1(i,3);...
                                DMT1(i,3); PMT5(i,3); DMT5(i,3);...
                                CUN(i,3); HLX(i,3)]);
    end
end