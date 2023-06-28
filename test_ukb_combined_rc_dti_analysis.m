%% This scripts tests whethter the current results in their essence equal the output of the original scripts

% Run current version of the analysis
ukb_combined_rc_dti_analysis

%% Load tested data
load(['/slow/projects/01_UKB/dti/rc_results_', datestr(now, 'yyyymmdd'), '_TEST.mat']);
subjectIDs = [1000010, 1000028];
tivTested = nan(length(subjectIDs), 1);
meanFaTested = nan(length(subjectIDs), 1);
for i_subjectID = 1:length(subjectIDs)
    subjectID = subjectIDs[i_subjectID];
    load(['/slow/projects/01_UKB/dti/richclub_adj_', num2str(subjectID), '_TEST.mat']);

    pathToBase = '/slow/projects/01_UKB/surface/00_batch1';
    pathToTiv = fullfile(pathToBase, num2str(subjectID), 'tiv.txt');
    tivTested(i_subjectID) = fgetl(fopen(fullfile(pathToBase, num2str(subjectID), 'tiv.txt')));
    meanFaTested(i_subjectID) = fgetl(fopen(fullfile(pathToBase, num2str(subjectID), 'meanfa.txt')));
end

%% Compare results
for i_subjectID = 1:length(subjectIDs)
    subjectID = subjectIDs[i_subjectID];
    pathToTiv = fullfile(pathToBase, num2str(subjectID), 'tiv.txt');
    isequal(tivTested, tiv)
end

