% -------------------------------------------------------------------------
% Initialisation
% -------------------------------------------------------------------------
clearvars;
close all;
clc;

% Set folders
toolboxFolder = 'C:\Users\moissene\Documents\Professionnel\routines\github\GaitEventsAlgorithms';
btkFolder     = 'C:\Users\moissene\Documents\Professionnel\routines\btk';
patientFolder = 'C:\Users\moissene\Documents\Professionnel\publications\communications\2020\SOFAMEA2020\Patients';
subjectFolder = 'C:\Users\moissene\Documents\Professionnel\publications\communications\2020\SOFAMEA2020\Healthy subjects';
exportFolder  = 'C:\Users\moissene\Documents\Professionnel\publications\communications\2020\SOFAMEA2020\Outputs';
addpath(toolboxFolder);
addpath(btkFolder);

% -------------------------------------------------------------------------
% Load files per condition / case
% -------------------------------------------------------------------------
conf  = {'R_conf1' 'R_conf2' 'R_conf3' 'R_conf4' ...
         'R_confCGM10' 'R_confCGM24a' 'R_confCGM24b' 'R_confCGM24c' ...
         'R_confLEARa' 'R_confLEARb' 'R_confLEARc'};
acase = {'case1' 'case2' 'case3'};
cd(exportFolder);
for i = 1:length(conf)
    for j = 1:length(acase)
%         conditionFiles.(conf{i}).(acase{j}) = dir(['*',conf{i},'*',acase{j},'*walkNorm*','.c3d']); % Healthy subjects
        conditionFiles.(conf{i}).(acase{j}) = dir(['*',conf{i},'*',acase{j},'*EVT*','.c3d']);      % Patients
        for k = 1:length(conditionFiles.(conf{i}).(acase{j}))
            btkData.(conf{i}).(acase{j})(k) = btkReadAcquisition(conditionFiles.(conf{i}).(acase{j})(k).name);
        end
    end
end

% -------------------------------------------------------------------------
% Extract results from C3D metadata
% -------------------------------------------------------------------------
for i = 1:length(conf)
    tconfFS1 = [];
    tconfFS2 = [];
    tconfFS3 = [];
    tconfFO1 = [];
    tconfFO2 = [];
    tconfFO3 = [];
    for j = 1:length(acase)
        for k = 1:length(conditionFiles.(conf{i}).(acase{j}))
            temp = btkGetMetaData(btkData.(conf{i}).(acase{j})(k));
            FS.(conf{i}).(acase{j})(k).data = temp.children.EventDetectionErrors.children.FS.info.values;
            FO.(conf{i}).(acase{j})(k).data = temp.children.EventDetectionErrors.children.FO.info.values;
            if j == 1
                tconfFS1 = [tconfFS1 FS.(conf{i}).(acase{j})(k).data];
                tconfFO1 = [tconfFO1 FO.(conf{i}).(acase{j})(k).data];
            elseif j == 2
                tconfFS2 = [tconfFS2 FS.(conf{i}).(acase{j})(k).data];
                tconfFO2 = [tconfFO2 FO.(conf{i}).(acase{j})(k).data];
            elseif j == 3
                tconfFS3 = [tconfFS3 FS.(conf{i}).(acase{j})(k).data];
                tconfFO3 = [tconfFO3 FO.(conf{i}).(acase{j})(k).data];
            end
            clear temp;
        end
    end
    FS.(conf{i}).case1Mean = mean(tconfFS1);
    FS.(conf{i}).case1Std = std(tconfFS1);
    FS.(conf{i}).case2Mean = mean(tconfFS2);
    FS.(conf{i}).case2Std = std(tconfFS2);
    FS.(conf{i}).case3Mean = mean(tconfFS3);
    FS.(conf{i}).case3Std = std(tconfFS3);
    FO.(conf{i}).case1Mean = mean(tconfFO1);
    FO.(conf{i}).case1Std = std(tconfFO1);
    FO.(conf{i}).case2Mean = mean(tconfFO2);
    FO.(conf{i}).case2Std = std(tconfFO2);
    FO.(conf{i}).case3Mean = mean(tconfFO3);
    FO.(conf{i}).case3Std = std(tconfFO3);
end

% -------------------------------------------------------------------------
% Plot results
% -------------------------------------------------------------------------
figure;
hold on;
for i = 1:4
    subplot(1,4,i);
    title(['Config ',num2str(i)]);
    ylim([-30 30]);
    hold on;
    bar(1:3,[FS.(conf{i}).case1Mean FO.(conf{i}).case1Mean; ...
             FS.(conf{i}).case2Mean FO.(conf{i}).case2Mean; ...
             FS.(conf{i}).case3Mean FO.(conf{i}).case3Mean]);
    errorbar([FS.(conf{i}).case1Mean FO.(conf{i}).case1Mean; ...
              FS.(conf{i}).case2Mean FO.(conf{i}).case2Mean; ...
              FS.(conf{i}).case3Mean FO.(conf{i}).case3Mean], ...
             [FS.(conf{i}).case1Std FO.(conf{i}).case1Std; ...
              FS.(conf{i}).case2Std FO.(conf{i}).case2Std; ...
              FS.(conf{i}).case3Std FO.(conf{i}).case3Std],'.');
    XTickLabel = {'Case 1' ; 'Case 2'; 'Case 3'};
    XTick = 1:3;
    set(gca,'XTick',XTick);
    set(gca,'XTickLabel',XTickLabel);     
end