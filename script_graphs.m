load('C:\Users\Ganglabor\Documents\automatisiertenEvents\Matlab\output\results_norm.mat')
load('C:\Users\Ganglabor\Documents\automatisiertenEvents\Matlab\output\results_patients.mat')

% -------------------------------------------------------------------------
%% Sole angle
% -------------------------------------------------------------------------
%% Sole angle midstance, norm
flatfoot_solAng_mean = results_norm.soleAng_midStance_mean; % =~0
flatfoot_solAng_SD = results_norm.soleAng_midStance_std; % =~1
% flatfoot_CI95 = std(results_norm.soleAngle_midStance,0.2) ...
%                 ./sqrt(length(results_norm.soleAngle_midStance));
[Q_solAng] = quantile(results_norm.soleAngle_midStance,[0.025 0.50 0.975]);

dsf = figure;
annotation(dsf,'textbox',[0.4 0.95 0.2 0.05],'String','Sole angle', ...
            'FontSize',15,'FontWeight','bold','HorizontalAlignment','center', ...
            'LineStyle','none')
subplot(2,2,1)
edges = [-4 -3:0.5:2 4];
histogram(results_norm.soleAngle_midStance,edges)
title('Norm, mean over mid-Stance')
hold on
% mean
l1=line([flatfoot_solAng_mean flatfoot_solAng_mean],[0 18],'Color','r','Linewidth',1);
% +1 SD
line([flatfoot_solAng_mean+flatfoot_solAng_SD flatfoot_solAng_mean+flatfoot_solAng_SD],...
        [0 18],'Color','r','Linewidth',1,'LineStyle','--')
% +2 SD
line([flatfoot_solAng_mean+2*flatfoot_solAng_SD flatfoot_solAng_mean+2*flatfoot_solAng_SD],...
    [0 18],'Color','r','Linewidth',1,'LineStyle','--')
% -1 SD
line([flatfoot_solAng_mean-flatfoot_solAng_SD flatfoot_solAng_mean-flatfoot_solAng_SD],...
        [0 18],'Color','r','Linewidth',1,'LineStyle','--')
% -2 SD
line([flatfoot_solAng_mean-2*flatfoot_solAng_SD flatfoot_solAng_mean-2*flatfoot_solAng_SD],...
        [0 18],'Color','r','Linewidth',1,'LineStyle','--')
% median
l2=line([Q_solAng(2) Q_solAng(2)],[0 18],'Color','g','Linewidth',1);
% 2.5% quantile
line([Q_solAng(1) Q_solAng(1)],[0 18],'Color','g','Linewidth',1,'LineStyle','--')
% 97.5% quantile
line([Q_solAng(3) Q_solAng(3)],[0 18],'Color','g','Linewidth',1,'LineStyle','--')

legend([l1,l2],'mean +/- 1.&2. SD','median, 2.5%&97.5% quantiles')

%% sole angle before FS, norm
subplot(2,2,2)
edges = [5 flatfoot_solAng_mean+10*flatfoot_solAng_SD:2*flatfoot_solAng_SD:flatfoot_solAng_mean+28*flatfoot_solAng_SD 30];
histogram(results_norm.soleAngle_beforeFS,edges)
legend('width 2SD')
title('Norm, mean before FS')

%% sole angle before FS, CP patients
% % overview, 1SD edge
% figure
% edges = [flatfoot_solAng_mean-40*flatfoot_solAng_SD:flatfoot_solAng_SD:flatfoot_solAng_mean+25*flatfoot_solAng_SD];
% histogram(results_patients.soleAngle_beforeFS,edges)
% overview, 2SD edge
subplot(2,2,3)
edges = [-40 flatfoot_solAng_mean-40*flatfoot_solAng_SD:2*flatfoot_solAng_SD:flatfoot_solAng_mean+26*flatfoot_solAng_SD 30];
histogram(results_patients.soleAngle_beforeFS,edges)
legend('width 2SD')
title('Patients, mean before FS')
% groups
% 1st group (dorsiflexion): > mean + 2SD =~ 1.5° (4SD = 4.4°)
% 2nd group (flat foot): mean +/- 2SD
% 3rd group (mild equinus): < mean - 2SD =~ -2°
% 4th group (high equinus): < mean - 10SD =~ -10° plantarflexion (exact
% -9.75°)
for i=1:length(results_patients.soleAngle_beforeFS)
    if results_patients.soleAngle_beforeFS(i) > flatfoot_solAng_mean+2*flatfoot_solAng_SD
        grp(i) = 0;
    elseif results_patients.soleAngle_beforeFS(i) >= flatfoot_solAng_mean-2*flatfoot_solAng_SD
        grp(i) = 1;
    elseif results_patients.soleAngle_beforeFS(i) >= flatfoot_solAng_mean-10*flatfoot_solAng_SD
        grp(i) = 2;
    else
        grp(i) = 3;
    end
end
cat = categorical(grp,[3 2 1 0], ...
    {'high equinus','mild equinus','flat foot','dorsiflexion'});
subplot(2,2,4)
h = histogram(cat,'BarWidth',0.5);
title('Patients'' groups, mean before FS')
str = {'high equinus: < mean-10SD (-9.8°)','mild equinus: < mean-2SD (-2.2°)', ...
        'flat foot: [mean+/-2SD]','dorsiflexion: > mean+2SD (1.5°)'};
text(0.1,50,str,'FontSize',8,'LineStyle','none')
% edges = [-40 -37 flatfoot_solAng_mean-10*flatfoot_solAng_SD flatfoot_solAng_mean-2*flatfoot_solAng_SD ...
%         flatfoot_solAng_mean+2*flatfoot_solAng_SD 22 30];
% histogram(results_patients.soleAngle_beforeFS,edges)

% -------------------------------------------------------------------------
%% Varus angle
% -------------------------------------------------------------------------
%% Varus angle midstance, norm
flatfoot_varAng_mean = results_norm.varusAng_midStance_mean; % =~0
flatfoot_varAng_SD = results_norm.varusAng_midStance_std; % =~1
% flatfoot_CI95 = std(results_norm.soleAngle_midStance,0.2) ...
%                 ./sqrt(length(results_norm.soleAngle_midStance));
[Q_varAng] = quantile(results_norm.varusAngle_midStance,[0.025 0.50 0.975]);

varFS = figure;
annotation(varFS,'textbox',[0.3 0.95 0.4 0.05],'String','Varus angle (FS)', ...
            'FontSize',15,'FontWeight','bold','HorizontalAlignment','center', ...
            'LineStyle','none')
subplot(2,2,1)
edges = [0:1:15];
histogram(results_norm.varusAngle_midStance,edges)
axis([0 15 0 14])
title('Norm, mean over mid-Stance')
hold on
% mean
l1=line([flatfoot_varAng_mean flatfoot_varAng_mean],[0 12],'Color','r','Linewidth',1);
% +1 SD
line([flatfoot_varAng_mean+flatfoot_varAng_SD flatfoot_varAng_mean+flatfoot_varAng_SD],...
        [0 12],'Color','r','Linewidth',1,'LineStyle','--')
% +2 SD
line([flatfoot_varAng_mean+2*flatfoot_varAng_SD flatfoot_varAng_mean+2*flatfoot_varAng_SD],...
    [0 12],'Color','r','Linewidth',1,'LineStyle','--')
% -1 SD
line([flatfoot_varAng_mean-flatfoot_varAng_SD flatfoot_varAng_mean-flatfoot_varAng_SD],...
        [0 12],'Color','r','Linewidth',1,'LineStyle','--')
% -2 SD
line([flatfoot_varAng_mean-2*flatfoot_varAng_SD flatfoot_varAng_mean-2*flatfoot_varAng_SD],...
        [0 12],'Color','r','Linewidth',1,'LineStyle','--')
% median
l2=line([Q_varAng(2) Q_varAng(2)],[0 12],'Color','g','Linewidth',1);
% 2.5% quantile
line([Q_varAng(1) Q_varAng(1)],[0 12],'Color','g','Linewidth',1,'LineStyle','--')
% 97.5% quantile
line([Q_varAng(3) Q_varAng(3)],[0 12],'Color','g','Linewidth',1,'LineStyle','--')

legend([l1,l2],'mean +/- 1.&2. SD','median, 2.5%&97.5% quantiles')

%% varus angle before FS, norm
subplot(2,2,2)
edges = [0 flatfoot_varAng_mean-2*flatfoot_varAng_SD:2*flatfoot_varAng_SD:flatfoot_varAng_mean+12*flatfoot_varAng_SD 45];
histogram(results_norm.varusAngle_beforeFS,edges)
legend('width 2SD')
title('Norm, mean before FS')

%% varus angle before FS, CP patients
% % overview, 1SD edge
% figure
% edges = [flatfoot_varAng_mean-9*flatfoot_varAng_SD:flatfoot_varAng_SD:flatfoot_varAng_mean+12*flatfoot_varAng_SD];
% histogram(results_patients.varusAngle_beforeFS,edges)
% overview, 2SD edge
subplot(2,2,3)
edges = [-25 flatfoot_varAng_mean-10*flatfoot_varAng_SD:2*flatfoot_varAng_SD:flatfoot_varAng_mean+10*flatfoot_varAng_SD 40];
histogram(results_patients.varusAngle_beforeFS,edges)
axis([-25 40 0 40])
legend('width 2SD')
title('Patients, mean before FS')

% 1st group (high varus): > mean + 4SD =~ 18°
% 2nd group (mild varus): > mean + 2SD =~ 12°
% 3rd group (flat foot): mean +/- 2SD
% 4th group (mild valgus): < mean - 2SD =~ 0.6°
% 5th group (high valgus): < mean - 4SD =~ -5°
for i=1:length(results_patients.varusAngle_beforeFS)
    if results_patients.varusAngle_beforeFS(i) > flatfoot_varAng_mean+4*flatfoot_varAng_SD
        grp(i) = 0;
    elseif results_patients.varusAngle_beforeFS(i) > flatfoot_varAng_mean+2*flatfoot_varAng_SD
        grp(i) = 1;
    elseif results_patients.varusAngle_beforeFS(i) >= flatfoot_varAng_mean-2*flatfoot_varAng_SD
        grp(i) = 2;
    elseif results_patients.varusAngle_beforeFS(i) >= flatfoot_varAng_mean-4*flatfoot_varAng_SD
        grp(i) = 3;
    else
        grp(i) = 4;
    end
end
cat = categorical(grp,[4 3 2 1 0], ...
    {'high valgus','mild valgus','flat foot','mild varus','high varus'});
subplot(2,2,4)
h = histogram(cat,'BarWidth',0.5);
title('Patients'' groups, mean before FS')
str = {'high valgus: < mean-4SD (-5.2°)','mild valgus: < mean-2SD (0.6°)', ...
        'flat foot: [mean+/-2SD]','mild varus: > mean+2SD (12.2°)', ...
        'high varus: > mean+4SD (18.0°)'};
text(0.1,68,str,'FontSize',8,'LineStyle','none')

% edges = [min(results_patients.varusAngle_beforeFS) ...
%         flatfoot_varAng_mean-4*flatfoot_varAng_SD ...
%         flatfoot_varAng_mean-2*flatfoot_varAng_SD ...
%         flatfoot_varAng_mean+2*flatfoot_varAng_SD ...
%         flatfoot_varAng_mean+4*flatfoot_varAng_SD ...
%         max(results_patients.varusAngle_beforeFS)];
% histogram(results_patients.varusAngle_beforeFS,edges)

%% varus angle before FO, norm
varFO = figure;
annotation(varFO,'textbox',[0.3 0.95 0.4 0.05],'String','Varus angle (FO)', ...
            'FontSize',15,'FontWeight','bold','HorizontalAlignment','center', ...
            'LineStyle','none')
subplot(2,2,1)
edges = [0:1:15];
histogram(results_norm.varusAngle_midStance,edges)
axis([0 15 0 14])
title('Norm, mean over mid-Stance')
hold on
% mean
l1=line([flatfoot_varAng_mean flatfoot_varAng_mean],[0 12],'Color','r','Linewidth',1);
% +1 SD
line([flatfoot_varAng_mean+flatfoot_varAng_SD flatfoot_varAng_mean+flatfoot_varAng_SD],...
        [0 12],'Color','r','Linewidth',1,'LineStyle','--')
% +2 SD
line([flatfoot_varAng_mean+2*flatfoot_varAng_SD flatfoot_varAng_mean+2*flatfoot_varAng_SD],...
    [0 12],'Color','r','Linewidth',1,'LineStyle','--')
% -1 SD
line([flatfoot_varAng_mean-flatfoot_varAng_SD flatfoot_varAng_mean-flatfoot_varAng_SD],...
        [0 12],'Color','r','Linewidth',1,'LineStyle','--')
% -2 SD
line([flatfoot_varAng_mean-2*flatfoot_varAng_SD flatfoot_varAng_mean-2*flatfoot_varAng_SD],...
        [0 12],'Color','r','Linewidth',1,'LineStyle','--')
% median
l2=line([Q_varAng(2) Q_varAng(2)],[0 12],'Color','g','Linewidth',1);
% 2.5% quantile
line([Q_varAng(1) Q_varAng(1)],[0 12],'Color','g','Linewidth',1,'LineStyle','--')
% 97.5% quantile
line([Q_varAng(3) Q_varAng(3)],[0 12],'Color','g','Linewidth',1,'LineStyle','--')

legend([l1,l2],'mean +/- 1.&2. SD','median, 2.5%&97.5% quantiles')

%% varus angle before FS, norm
subplot(2,2,2)
edges = [-25 flatfoot_varAng_mean-10*flatfoot_varAng_SD:2*flatfoot_varAng_SD:flatfoot_varAng_mean+2*flatfoot_varAng_SD 15];
histogram(results_norm.varusAngle_beforeFO,edges)
legend('width 2SD')
title('Norm, mean before FO')

%% varus angle before FO, CP patients
% % overview, 1SD edge
% figure
% edges = [flatfoot_varAng_mean-21*flatfoot_varAng_SD:flatfoot_varAng_SD:flatfoot_varAng_mean+12*flatfoot_varAng_SD];
% histogram(results_patients.varusAngle_beforeFO,edges)
% overview, 2SD edge
subplot(2,2,3)
edges = [-60 flatfoot_varAng_mean-22*flatfoot_varAng_SD:2*flatfoot_varAng_SD:flatfoot_varAng_mean+12*flatfoot_varAng_SD 45];
histogram(results_patients.varusAngle_beforeFO,edges)
axis([-60 50 0 40])
legend('width 2SD')
title('Patients, mean before FO')

% 1st group (high varus): > mean + 4SD =~ 18°
% 2nd group (mild varus): > mean + 2SD =~ 12°
% 3rd group (flat foot): mean +/- 2SD
% 4th group (mild valgus): < mean - 2SD =~ 0.6°
% 5th group (high valgus): < mean - 4SD =~ -5°
for i=1:length(results_patients.varusAngle_beforeFO)
    if results_patients.varusAngle_beforeFO(i) > flatfoot_varAng_mean+4*flatfoot_varAng_SD
        grp(i) = 0;
    elseif results_patients.varusAngle_beforeFO(i) > flatfoot_varAng_mean+2*flatfoot_varAng_SD
        grp(i) = 1;
    elseif results_patients.varusAngle_beforeFO(i) >= flatfoot_varAng_mean-2*flatfoot_varAng_SD
        grp(i) = 2;
    elseif results_patients.varusAngle_beforeFO(i) >= flatfoot_varAng_mean-4*flatfoot_varAng_SD
        grp(i) = 3;
    else
        grp(i) = 4;
    end
end
cat = categorical(grp,[4 3 2 1 0], ...
    {'high valgus','mild valgus','flat foot','mild varus','high varus'});
subplot(2,2,4)
h = histogram(cat,'BarWidth',0.5);
title('Patients'' groups, mean before FO')
str = {'high valgus: < mean-4SD (-5.2°)','mild valgus: < mean-2SD (0.6°)', ...
        'flat foot: [mean+/-2SD]','mild varus: > mean+2SD (12.2°)', ...
        'high varus: > mean+4SD (18.0°)'};
text(0.1,52,str,'FontSize',8,'LineStyle','none')

% edges = [min(results_patients.varusAngle_beforeFS) ...
%         flatfoot_varAng_mean-4*flatfoot_varAng_SD ...
%         flatfoot_varAng_mean-2*flatfoot_varAng_SD ...
%         flatfoot_varAng_mean+2*flatfoot_varAng_SD ...
%         flatfoot_varAng_mean+4*flatfoot_varAng_SD ...
%         max(results_patients.varusAngle_beforeFS)];
% histogram(results_patients.varusAngle_beforeFO,edges)