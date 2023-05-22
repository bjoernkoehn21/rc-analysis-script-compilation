%% TODO
% change to real ID-file: subjectIDs.txt
% change number of random networks and parallel workers 2 Blocks below
% adapt line of total intrakranial volume: TIV_LINE_NUMBER
% adapt file name of file with subject IDs

% add paths to relevant packages
addpath /fast/home/marketts/tools/batches  % path to fdr_bh.m
addpath /home/marketts/imaging/2017_01_15_BCT % path to BCT toolbox
% Set up freesurfer environment
setenv('FREESURFER_HOME', '/opt/freesurfer');
system(fullfile(getenv('FREESURFER_HOME'), 'SetUpFreeSurfer.sh'));

% define destination of results and other relevant paths and file names
SUBJECTS_ID_FILE = 'ukb_subjIDs.txt';
PATH_TO_SUBJECT_DIRS = '/slow/projects/01_UKB/dti/00_batch1';
RESULTS_DESTINATION = '/slow/projects/01_UKB/dti/rc_results_220622_TEST.mat';
DENSITY_DESTINATION = '/slow/projects/01_UKB/dti/rcdense_TEST.mat';

% define constants
N_PARALLEL_WORKERS = 4; %100;
N_RANDOM_NETWORKS = 10; %2500;
TIV_LINE_NUMBER = 34; % refers to line number in SUBJECT/T1w/freesurfer/stats/aseg.stats containing total intrakranial volume
EDGE_WEIGHTS={'fa' 'svd'}; % fractional anisotropy, streamline volume density
FILES={'*_connectivity_csd_dti_aparc.mat'
'*_connectivity_csd_dti_lausanne120.mat'
'*_connectivity_csd_dti_lausanne250.mat'
'*_connectivity_gqi_dti_aparc.mat'
'*_connectivity_gqi_dti_lausanne120.mat'
'*_connectivity_gqi_dti_lausanne250.mat'};

% read in subjects IDs from subjIDs.txt
fid = fopen(SUBJECTS_ID_FILE, 'r');
subjIDs = textscan(fid, '%s');
subjIDs = subjIDs{1};
fclose(fid);

% preallocate space for assembled-rc-results matrics
rcResults = cell(length(subjIDs), 1);
rcResults(:) = {struct( ...
    'id', zeros(1, 1), ...
    'max_phi', zeros(length(FILES), length(EDGE_WEIGHTS)), ...
    'max_k', zeros(length(FILES), length(EDGE_WEIGHTS)), ...
    'range', zeros(2, length(FILES), length(EDGE_WEIGHTS)), ...
    'integral', struct( ...
        'norm', zeros(length(FILES), length(EDGE_WEIGHTS)), ...
        'emp', zeros(length(FILES), length(EDGE_WEIGHTS)), ...
        'rand', zeros(length(FILES), length(EDGE_WEIGHTS)), ...
        'above', zeros(length(FILES), length(EDGE_WEIGHTS))), ...
    'odd', zeros(4, length(FILES), length(EDGE_WEIGHTS)))};

% preallocate space for density-of-empirical-networks matrix
rcDensity = cell(length(subjIDs), 1);
rcDensity(:) = {struct( ...
    'fa', zeros(length(FILES), 1), ...
    'svd', zeros(length(FILES), 1))};


%% Start parallel pool
UKB_Cluster = parcluster('local');
UKB_Cluster.NumWorkers = N_PARALLEL_WORKERS;
saveProfile(UKB_Cluster);
parpool('local', N_PARALLEL_WORKERS);

% loop over subject IDs and execute ukb_compute_rc_dti function in parallel
parfor iSubj = 1:length(subjIDs)
%for iSubj = 1:length(subjIDs)
    hcp_compute_correct_rc_dti(iSubj, subjIDs{iSubj}, N_RANDOM_NETWORKS, EDGE_WEIGHTS, FILES, rcResults, rcDensity, TIV_LINE_NUMBER);
    
    savefile=RESULTS_DESTINATION
    save(savefile, 'rcResults');
    
    savefile=DENSITY_DESTINATION
    save(savefile, 'rcDensity');
end

%delete(gcp('nocreate')); % close interactive sessions

%% compute rc with preallocation for random networks and correct the rc
function hcp_compute_correct_rc_dti(iSubj, subjID, N_RANDOM_NETWORKS, EDGE_WEIGHTS, FILES, rcResults, rcDensity, TIV_LINE_NUMBER)
    %resultsDir = fullfile('/slow/projects/01_UKB/dti/00_batch1', subjID, 'DWI_processed_v311');
    resultsDir = fullfile(PATH_TO_SUBJECT_DIRS, subjID, 'DWI_processed_v311');
    cd(resultsDir);

%% compute rc coefficient and p-values
    tic

    for iFile=1:length(FILES) % loop over files
        try
            n=dir(FILES{iFile}); load(n.name) % load file
            
            % calculate rc coefficients for empirical network
            rc{iFile}.fa.emp=rich_club_wu(connectivity(:,:,3));
            rc{iFile}.svd.emp=rich_club_wu(connectivity(:,:,13));
            
            % calculate
            rcDensity{iSubj}.fa(iFile)=density_und(connectivity(:,:,3));
            rcDensity{iSubj}.svd(iFile)=density_und(connectivity(:,:,13));
            
            % preallocate (size of rc{iFile}.fa.emp not known beforehand [before executing rich_club_wu])
            rc{iFile}.fa.rand=zeros(length(rc{iFile}.fa.emp), N_RANDOM_NETWORKS);
            rc{iFile}.svd.rand=zeros(length(rc{iFile}.svd.emp), N_RANDOM_NETWORKS);
            rc{iFile}.fa.pvals=zeros(1,1);
            rc{iFile}.svd.pvals=zeros(1,1);
        catch ME
            fprintf('\t Error in line %d in function %s: %s\n', ME.stack(1).line, ME.stack(1).name, ME.message);
            rc{iFile}=[];
        end
    end

    
    for iFile = 1:length(FILES) % loop over files
            for iRandomNetwork = 1:N_RANDOM_NETWORKS
                rc{iFile}.fa.rand(:,iRandomNetwork)=rich_club_wu(randmio_und(connectivity(:,:,3),10),length(rc{iFile}.fa.emp));
                rc{iFile}.svd.rand(:,iRandomNetwork)=rich_club_wu(randmio_und(connectivity(:,:,13),10),length(rc{iFile}.svd.emp));
            end

            rc{iFile}.fa.norm=rc{iFile}.fa.emp./mean(rc{iFile}.fa.rand,2)';
            rc{iFile}.svd.norm=rc{iFile}.svd.emp./mean(rc{iFile}.svd.rand,2)';

            rc{iFile}.fa.pvals=1-sum(rc{iFile}.fa.emp>rc{iFile}.fa.rand')/N_RANDOM_NETWORKS;
            rc{iFile}.svd.pvals=1-sum(rc{iFile}.svd.emp>rc{iFile}.svd.rand')/N_RANDOM_NETWORKS;
            % mean seems to take longer
            %rc{iFile}.fa.pvals=1 - mean(rc{1}.fa.emp > rc{1}.fa.rand', 1);
            %rc{iFile}.svd.pvals=1 - mean(rc{1}.svd.emp > rc{1}.svd.rand', 1);
            
            
%% correct rc range in cases where there is no unambiguous rc regime
        try
            for iEdgeWeight = 1:length(EDGE_WEIGHTS)
              % phi max and corresponding k  
              [rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).kmax(1,1),rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).kmax(1,2)]=max(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm) ;

              % fdr-correction of p-values (adjusted p-values)
              %[~, ~,
              %rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).adj_p]=fdr_bh(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).pvals,.05,'pdep','false');
              % original by Markett
               rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).adj_p=fdr_bh(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).pvals,.05,'pdep','false');
              %%
              % rc range
                %x=1:length(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).adj_p); xp=x(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).adj_p<.05); %helper
                xp = find(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).adj_p < 0.05); % more memory-efficient than line above
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range = [min(xp) max(xp)]; % rc range
                        idx =rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).adj_p(min(xp):max(xp)); % identify values with p>=.05 within range
                        n = length(idx);
                        idx1 = (idx > 0.05); %idx1 contains positions in RC range where p > 0.05
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Outliers = sum(idx1); rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).RangeCorrected = 0;   rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Ignored =  0;

                    if sum(idx1) > 1 % if there's more than 1 outlier    
                        %Range of RC range between minimum and maximum K in RC range > 0.05
                        x_RC = 1:n; x_RC_p = x_RC(idx1);
                        idx_bigger = idx(1,min(x_RC_p):max(x_RC_p)); 
                        idx2 = (idx_bigger < 0.05);
                        %If range of p-values above 0.05 is continous, select bigger RC
                        %range and define it as adjusted RC range
                        if sum(idx2) <= 1 
                           a = length(idx(1,1:(min(x_RC_p)-1)));
                           b = length(idx(1,(max(x_RC_p)+1):n));
                           if a > (1.5*b)
                               idx_cor = zeros(1,n);
                               idx_cor(1,1:(min(x_RC_p)-1)) = 1;
                               idx_cor = logical(idx_cor);
                               x_cor = min(xp):max(xp);
                               x_cor_p = x_cor(idx_cor);
                               rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range= [min(x_cor_p) max(x_cor_p)];
                               rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).RangeCorrected = 1; %Count of of subjects with corrections 
                           elseif b > (1.5*a)
                               idx_cor = zeros(1,n);
                               idx_cor(1,(max(x_RC_p)+1):n) = 1;
                               idx_cor = logical(idx_cor);
                               x_cor = min(xp):max(xp);
                               x_cor_p = x_cor(idx_cor);
                               rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range = [min(x_cor_p) max(x_cor_p)];
                               rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).RangeCorrected = 1; %Count of of subjects with corrections
                           %If no RC range is 1.5 times bigger than other, select RC
                           %range with higher K values
                           else
                               idx_cor = zeros(1,n);
                               idx_cor(1,(max(x_RC_p)+1):n) = 1;
                               idx_cor = logical(idx_cor);
                               x_cor = min(xp):max(xp);
                               x_cor_p = x_cor(idx_cor);
                                rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range = [min(x_cor_p) max(x_cor_p)];
                               rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).RangeCorrected = 1; 
                           end

                        elseif (length(idx2) - sum(idx2)) == 2 %If only two not
                            %continous p-values are above 0.05, ignore them
                            rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Ignored =  1;
                        else %If range of p-values above 0.05 is not continous, exclude subjects from analysis
                            rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range = [NaN NaN];
                        end
                    elseif sum(idx1) == 1
                        rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Outliers(2)=1; %Count number of subjects with single outliers
                    end
                    
                    % assemble relevant metrics
                    rcResults{iSubj}.id=subjID;
                    rcResults{iSubj}.max_phi(iFile,iEdgeWeight)=rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).kmax(1);
                    rcResults{iSubj}.max_k(iFile,iEdgeWeight)=rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).kmax(2);
                    rcResults{iSubj}.range(:,iFile,iEdgeWeight)=rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range;


                    rcResults{iSubj}.integral.norm(iFile,iEdgeWeight) = trapz(find(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm)),rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm)));
                    rcResults{iSubj}.integral.emp(iFile,iEdgeWeight) = trapz(find(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).emp)),rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).emp(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).emp)));
                    rcResults{iSubj}.integral.rand(iFile,iEdgeWeight) = trapz(find(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand)),rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand)));
                    rcResults{iSubj}.integral.above(iFile,iEdgeWeight) = rcResults{iSubj}.integral.norm(iFile,iEdgeWeight)-length(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm)));
                    
                    rcResults{iSubj}.odd(1,iFile,iEdgeWeight)=rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Outliers(1);
                    rcResults{iSubj}.odd(2,iFile,iEdgeWeight)=rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).RangeCorrected;
                    rcResults{iSubj}.odd(3,iFile,iEdgeWeight)=rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Ignored;
                    rcResults{iSubj}.odd(4,iFile,iEdgeWeight)=length(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range); 
            end
        catch ME
            fprintf('\t Error in line %d in function %s: %s\n', ME.stack(1).line, ME.stack(1).name, ME.message);
        end
    end
    
    try
%% read transkranial volume from FreeSurfer file 'stats/aseg.stats'
        %resultsDir = fullfile('C:\Users\askou\Desktop\Geschaeftliches\HU_Markett-Aktueller_Job\HCP-data_download_and_CATO\MATLAB_analysis', num2str(subjID), 'DWI_processed_v311');
        fid = fopen('../T1w/freesurfer/stats/aseg.stats');
        
        for ii = 1:(TIV_LINE_NUMBER-1) % ignore lines 1 to (TIV_LINE_NUMBER-1) (faster and without as much storing than/as with textscan)
            fgetl(fid);
        end

        tline = fgetl(fid); %read line TIV_LINE_NUMBER with relevant tiv information
        fclose(fid);

        tiv = regexp(tline, '\d+\.?\d*', 'match'); % read number (with or without decimal point) from line 35 
        tiv = strrep(tiv, ',', ''); % delete comma
        % write number of tiv into tiv.txt
        fid = fopen('tiv_TEST.txt', 'w');
        fprintf(fid, '%s', tiv{1});
        fclose(fid);
%% calculate mean FA    
    % Binarize aparc+aseg.mgz and save as wm.mask.mgz
    system('mri_binarize --i ../T1w/freesurfer/mri/aparc+aseg.mgz --wm --o wm.mask.mgz');
 
    % Convert wm.mask.mgz to wm.nii.gz using nearest-neighbor interpolation
    system(['/opt/freesurfer/bin/mri_convert -rt nearest -rl ' resultsDir '/*_fractional_anisotropy.nii.gz wm.mask.mgz wm.nii.gz']);
 
    % Compute mean FA within wm.nii.gz mask and save as meanfa.txt
    system(['fslmeants -i ' resultsDir '/*_fractional_anisotropy.nii.gz -m wm.nii.gz -o meanfa.txt']);
 
    % Clean up
    system('rm -f wm*');
    catch ME
        fprintf('\t Error in line %d in function %s: %s\n', ME.stack(1).line, ME.stack(1).name, ME.message);
    end

    savefile='richclub_adj_TEST.mat'; % adjusted richclub
    save(savefile,'rc')
    toc
end
