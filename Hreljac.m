% define gait events (foot strike and foot off) from Hreljac method
% input:
% - footMk: 3D coordinates of a foot marker
% - pelvicMk: a struct of the 5 pelvic markers (LASI; RASI; LPSI; RPSI;
% SACR)
% - gaitAxis: number of the gait axis: 1=x, 2=y
% - verticalAxis: number of the vertical axis (3=z)
% - f: frequence
%
% output:
% - FS = array of foot strikes
% - FO = array of foot offs
function [FS,FO] = Hreljac(footMk,pelvicMk,gaitAxis,verticalAxis,f)
FS = [];
FO = [];
% -------------------------------------------------------------------------
    % Foot strike
    % jerk of the marker in the vertical component = 0
% -------------------------------------------------------------------------
    velZ_footMk = diff(footMk(:,verticalAxis))/(1/f);
    accZ_footMk = diff(velZ_footMk)/(1/f);
    jerkZ_footMk = diff(accZ_footMk)/(1/f);
    
    % define the events from Zeni
    [FS_zeni,FO_zeni] = Zeni(footMk,pelvicMk,gaitAxis,f);
    
    % find the first max peak of acceleration of the vertical component after
    % the peaks from Zeni (= first zero of the jerk, positiv to negativ)
    % (visually determined from norm data)
    for i=1:length(FS_zeni)
        indJerk = FS_zeni(i);
        while ( indJerk<length(jerkZ_footMk) && jerkZ_footMk(indJerk)<=0 )
            indJerk = indJerk + 1;
        end % indJerk=length(jerkZ_footMk) or jerkZ_footMk(indJerk)>0
        while ( indJerk<length(jerkZ_footMk) && jerkZ_footMk(indJerk+1)>=0 )
            indJerk = indJerk + 1;
        end % indJerk=length(jerkZ_footMk) or jerkZ_footMk(indJerk+1)<0
        if indJerk < length(jerkZ_footMk) % i.e. jerkZ_footMk(indJerk+1)<0
            t2 = indJerk + 1;
            if jerkZ_footMk(indJerk)>0
                t1 = indJerk;
            else % jerkZ_footMk(indJerk)=0
                while jerkZ_footMk(indJerk)==0
                    indJerk = indJerk - 1;
                end % jerkZ_footMk(indJerk)>0
                t1 = indJerk;
            end
        else % indJerk >= length(jerkZ_footMk)
            t1 = NaN;
            t2 = NaN;
        end
        if ~isnan(t1) && ~isnan(t2)
            FS(i) = t1 + jerkZ_footMk(t1) / (jerkZ_footMk(t1)-jerkZ_footMk(t2));
        else
            FS(i) = NaN;
        end
    end
    
% plot(footMk(1:end-10,verticalAxis)*1000)
% hold on
% plot(accZ_footMk(1:end-10)*20)
% plot(jerkZ_footMk(1:end-10))
% for i=1:length(FS_zeni)
%     line([FS_zeni(i) FS_zeni(i)],[-1000000 1000000],'color','r','LineWidth',1)
% end

% -------------------------------------------------------------------------
    % Foot Off
    % jerk of the marker in the gait direction
% -------------------------------------------------------------------------
    velX_footMk = diff(footMk(:,gaitAxis))/(1/f);
    accX_footMk = diff(velX_footMk)/(1/f);
    jerkX_footMk = diff(accX_footMk)/(1/f);
    
    % define the max of acceleration in a window about the FO_zeni events
    % ([-f/5 f/5])
    for i=1:length(FO_zeni)
        begin_ = max(1,FO_zeni(i)-f/5);
        end_ = min(length(accX_footMk),FO_zeni(i)+f/5);
        [~,ind_max_accX] = max(accX_footMk(begin_:end_));
        % find the zero of jerk corresponding to this max (search from this
        % max - f/15 frames)
        indJerk = max(1,begin_ + ind_max_accX - f/15);
        while ( indJerk<length(jerkX_footMk) && jerkX_footMk(indJerk)<=0 )
            indJerk = indJerk + 1;
        end % indJerk=length(jerkX_footMk) or jerkX_footMk(indJerk)>0
        while ( indJerk<length(jerkX_footMk) && jerkX_footMk(indJerk+1)>=0 )
            indJerk = indJerk + 1;
        end % indJerk=length(jerkX_footMk) or jerkX_footMk(indJerk+1)<0
        if indJerk < length(jerkX_footMk) % i.e. jerkX_footMk(indJerk+1)<0
            t2 = indJerk + 1;
            if jerkX_footMk(indJerk)>0
                t1 = indJerk;
            else % jerkX_footMk(indJerk)=0
                while jerkX_footMk(indJerk)==0
                    indJerk = indJerk - 1;
                end % jerkX_footMk(indJerk)>0
                t1 = indJerk;
            end
        else % indJerk >= length(jerkX_footMk)
            t1 = NaN;
            t2 = NaN;
        end
        if ~isnan(t1) && ~isnan(t2)
            FO(i) = t1 + jerkX_footMk(t1) / (jerkX_footMk(t1)-jerkX_footMk(t2));
        else
            FO(i) = NaN;
        end
    end
    
    
% i=2;
% plot(footMk(FO_zeni(i)-f/5:FO_zeni(i)+f/5,gaitAxis)*1000)
% hold on
% line([f/5 f/5],[-200000 200000])
% plot(accX_footMk(FO_zeni(i)-f/5:FO_zeni(i)+f/5)*20)
% plot(jerkX_footMk(FO_zeni(i)-f/5:FO_zeni(i)+f/5))

% plot(footMk(1:end-10,gaitAxis)*1000)
% hold on
% plot(accX_footMk(1:end-10)*20)
% plot(jerkX_footMk(1:end-10))
% for i=1:length(FO_zeni)
%     line([FO_zeni(i) FO_zeni(i)],[-1000000 1000000],'color','r','LineWidth',1)
% end
