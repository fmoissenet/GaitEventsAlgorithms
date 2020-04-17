% function to calculate the error between the measured events and the
% estimated events for the first gait cycle after the general event defined
% in the c3d
% Input:
% - btkData = btk handle of the c3d file
% - mFS,mFO = measured foot strikes and foot offs (in frames)
% - eFS,eFO = estimated foot strikes and foot offs (vector, in frames)
% - f: frequency of the data (markers)
% - side: leg side ('L' or 'R')
%
% Output:
% values in frames by 150Hz
% - errorFS = error for the foot strike between estimated and measured
% values
% - errorFO = error for the foot off between estimated and measured values

function [errorFS,errorFO] = getError(btkData,mFS,mFO,eFS,eFO,f,side)
    % Get the estimated gait events related to the forceplate
    [FS,FO] = getEstimatedGaitEvents(btkData,eFS,eFO,f,side);
    % !!! the FS and FO are given in frames! some measurements are done with 300Hz and not 150Hz!
    % give all the results in 150Hz : divided by 2 = 150/300
    if f==300
        errorFS = (FS-mFS)/2;
        errorFO = (FO-mFO)/2;
    else
        errorFS = FS-mFS;
        errorFO = FO-mFO;
    end