%% Last step Changes
% remove "TEST" from names of result matrices
% change to real ID-file: subjectIDs.txt
% change number of random networks and parallel workers 2 blocks below to higher number

%% Preamble
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
TIV_LINE_NUMBER = 35; % line number in freesurfer/stats/aseg.stats; TIV=total intracranial volume
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

% preallocate space for assembled-rc-results matrics before loop for time reasons
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

% preallocate space for density-of-empirical-networks matrix before loop for time reasons
rcDensity = cell(length(subjIDs), 1);
rcDensity(:) = {struct( ...
    'fa', zeros(length(FILES), 1), ...
    'svd', zeros(length(FILES), 1))};

%% Main process 

% Start parallel pool
UKB_Cluster = parcluster('local');
UKB_Cluster.NumWorkers = N_PARALLEL_WORKERS;
saveProfile(UKB_Cluster);
parpool('local', N_PARALLEL_WORKERS);

% loop over subject IDs and execute compute_correct_assemble_rc_dti function in parallel
parfor iSubj = 1:length(subjIDs)
    tic
    compute_correct_assemble_rc_dti(iSubj, subjIDs{iSubj}, N_RANDOM_NETWORKS, EDGE_WEIGHTS, ...
        FILES, rcResults, rcDensity, TIV_LINE_NUMBER);
    toc
end

%% subfunction
function compute_correct_assemble_rc_dti(iSubj, subjID, N_RANDOM_NETWORKS, EDGE_WEIGHTS, ...
    FILES, rcResults, rcDensity, TIV_LINE_NUMBER)
    %% compute rc coefficient and p-values
    resultsDir = fullfile(PATH_TO_SUBJECT_DIRS, subjID, 'DWI_processed_v311');
    cd(resultsDir);
    for iFile = 1:length(FILES)
        try
            n=dir(FILES{iFile}); load(n.name) % load file
            
            % calculate rc coefficients for empirical network
            rc{iFile}.fa.emp=rich_club_wu(connectivity(:,:,3));
            rc{iFile}.svd.emp=rich_club_wu(connectivity(:,:,13));
            
            % calculate densities
            rcDensity{iSubj}.fa(iFile)=density_und(connectivity(:,:,3));
            rcDensity{iSubj}.svd(iFile)=density_und(connectivity(:,:,13));
            
            % preallocate (size not known beforehand [before executing rich_club_wu])
            rc{iFile}.fa.rand=zeros(length(rc{iFile}.fa.emp), N_RANDOM_NETWORKS);
            rc{iFile}.svd.rand=zeros(length(rc{iFile}.svd.emp), N_RANDOM_NETWORKS);
            rc{iFile}.fa.pvals=zeros(1,1);
            rc{iFile}.svd.pvals=zeros(1,1);
        catch ME
            fprintf('\t Error in line %d in function %s: %s\n', ...
                ME.stack(1).line, ME.stack(1).name, ME.message);
            rc{iFile} = [];
        end
    end

    
    for iFile = 1:length(FILES)
        for iRandomNetwork = 1:N_RANDOM_NETWORKS
            rc{iFile}.fa.rand(:,iRandomNetwork) = rich_club_wu( ...
                randmio_und(connectivity(:,:,3),10),length(rc{iFile}.fa.emp));
            rc{iFile}.svd.rand(:,iRandomNetwork) = rich_club_wu( ...
                randmio_und(connectivity(:,:,13),10),length(rc{iFile}.svd.emp));
        end
        
        % normalize
        rc{iFile}.fa.norm = rc{iFile}.fa.emp./mean(rc{iFile}.fa.rand,2)';
        rc{iFile}.svd.norm = rc{iFile}.svd.emp./mean(rc{iFile}.svd.rand,2)';

        % calculate p-values
        rc{iFile}.fa.pvals = 1 - mean(rc{iFile}.fa.emp > rc{iFile}.fa.rand', 1);
        rc{iFile}.svd.pvals = 1 - mean(rc{iFile}.svd.emp > rc{iFile}.svd.rand', 1);
            
            
        %% correct rc range in cases where there is no unambiguous rc regime
        try
            for iEdgeWeight = 1:length(EDGE_WEIGHTS)
                % phi max and corresponding k  
                [rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).kmax(1,1), ...
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).kmax(1,2)] = max( ...
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm) ;
                
                % fdr-correction of p-values (adjusted p-values)
                rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).adj_p = fdr_bh( ...
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).pvals,.05,'pdep','false');
    
              xp = find(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).adj_p < 0.05); %%%%%%%%%%%%%%%%%%%%%besser eindeutigetr Variablenname: rcRangeIndexes 
              rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range = [min(xp) max(xp)]; % rc range
              
              % identify values with p>=0.05 within range
              idx = rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).adj_p(min(xp):max(xp)); %%%%%%%%%%%%%% besser eindeutigetr Variablenname: rcRangePvalues
              n = length(idx); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% besser eindeutigetr Variablenname: rcRangeLength
              idx1 = (idx > 0.05); % idx1 contains positions in RC range where p > 0.05 %%%%%%%%%%%%% besser eindeutigetr Variablenname: outlierInRcRangeBool
			  
              
              rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Outliers = sum(idx1);
              rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).RangeCorrected = 0;
              rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Ignored =  0;
                    
              if sum(idx1) > 1 % if there's more than 1 outlier
                  
                  % Range of RC range between minimum and maximum K in RC range > 0.05
                  x_RC = 1:n; x_RC_p = x_RC(idx1); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% einfacher und leichter verständlich: outliersIndexesInRcRange = find(idx1) 
                  idx_bigger = idx(1,min(x_RC_p):max(x_RC_p)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% pvaluesFirstOutlierToLastOutlier
                  idx2 = (idx_bigger < 0.05); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% rcInOutlierRangeBool
                  
                  % If range of p-values above 0.05 is continous, select bigger RC-range and define 
                  % it as adjusted RC range
                  if sum(idx2) <= 1
                      a = length(idx(1,1:(min(x_RC_p)-1))); %% The first index is nnecessary in a 1D Matrix %%%%%%%%%%%%%% nRcEntriesBeforeFirstOutlier
                      b = length(idx(1,(max(x_RC_p)+1):n)); %%%%%%%%%%%%%% nRcEntriesAfterLastOutlier
                      if a > (1.5*b) % most rc entries before outliers
                          idx_cor = zeros(1,n);	%%%%% This line and others are executed in all alternative routes --> move this out of the if-structure
                          idx_cor(1,1:(min(x_RC_p)-1)) = 1;
                          idx_cor = logical(idx_cor);
                          x_cor = min(xp):max(xp);
                          x_cor_p = x_cor(idx_cor);
                          rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range= [min(x_cor_p) max(x_cor_p)];
                          rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).RangeCorrected = 1;
                      elseif b > (1.5*a) % most rc entries after outliers %% The routine here is the same as in the else-case --> melt into one case: the else case
                          idx_cor = zeros(1,n);
                          idx_cor(1,(max(x_RC_p)+1):n) = 1;
                          idx_cor = logical(idx_cor);
                          x_cor = min(xp):max(xp);
                          x_cor_p = x_cor(idx_cor);
                          rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range = [min(x_cor_p) max(x_cor_p)];
                          rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).RangeCorrected = 1;
                          
                      % If no RC range is 1.5 times bigger than the other, select RC-range with 
                      % higher K values
                      else
                          idx_cor = zeros(1,n);
                          idx_cor(1,(max(x_RC_p)+1):n) = 1;
                          idx_cor = logical(idx_cor);
                          x_cor = min(xp):max(xp);
                          x_cor_p = x_cor(idx_cor);
                          rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range = [min(x_cor_p) max(x_cor_p)];
                          rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).RangeCorrected = 1; 
                      end
                  % If there is only two outliers all continous p-values above 0.05 are ignored %%%% Ist das einfach eine random Entscheidung, weil man eben irgendeinen cutoff braucht?
                  elseif (length(idx2) - sum(idx2)) == 2 %% Könnte man hier nicht einfacher sum(idx1) == 2 abfragen? analog  zu erster Abfrage mit "if sum(idx1) > 1 % if there's more than 1 outlier"
                      rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Ignored =  1;
                  % If range of p-values above 0.05 is not continous, exclude subjects from analysis
                  else
                      rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range = [NaN NaN];
                  end
              elseif sum(idx1) == 1 %%%%%% Dieser Fall ist doch eigentlich schon abgehandelt mit "rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Outliers = sum(idx1);", oder? Wieso wird hier noch ein zweiter Eintrag im Outliers-Feld gemacht (Outliers(2))?
                  rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Outliers(2) = 1; 
              end
              
              % assemble relevant metrics
              rcResults{iSubj}.id = subjID;
              rcResults{iSubj}.max_phi(iFile,iEdgeWeight) = rc{iFile}.( ...
                  EDGE_WEIGHTS{iEdgeWeight}).kmax(1);
              rcResults{iSubj}.max_k(iFile,iEdgeWeight) = rc{iFile}.( ...
                  EDGE_WEIGHTS{iEdgeWeight}).kmax(2);
              rcResults{iSubj}.range(:,iFile,iEdgeWeight) = rc{iFile}.( ...
                  EDGE_WEIGHTS{iEdgeWeight}).range;
              
              rcResults{iSubj}.integral.norm(iFile,iEdgeWeight) = trapz( ...
                  ind(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm)), ...
                  rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm( ...
                  ~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm)));
              rcResults{iSubj}.integral.emp(iFile,iEdgeWeight) = trapz( ...
                  find(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).emp)), ...
                  rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).emp( ...
                  ~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).emp)));
              rcResults{iSubj}.integral.rand(iFile,iEdgeWeight) = trapz( ...
                  find(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand)), ...
                  rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand( ...
                  ~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand)));
              rcResults{iSubj}.integral.above(iFile,iEdgeWeight) = rcResults{ ...
                  iSubj}.integral.norm(iFile,iEdgeWeight) - ...
                  length(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm( ...
                  ~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm)));
              
              rcResults{iSubj}.odd(1,iFile,iEdgeWeight) = rc{iFile}.( ...
                  EDGE_WEIGHTS{iEdgeWeight}).Outliers(1);
              rcResults{iSubj}.odd(2,iFile,iEdgeWeight) = rc{iFile}.( ...
                  EDGE_WEIGHTS{iEdgeWeight}).RangeCorrected;
              rcResults{iSubj}.odd(3,iFile,iEdgeWeight) = rc{iFile}.( ...
                  EDGE_WEIGHTS{iEdgeWeight}).Ignored;
              rcResults{iSubj}.odd(4,iFile,iEdgeWeight) = length(rc{iFile}.( ...
                  EDGE_WEIGHTS{iEdgeWeight}).range); 
            end
        catch ME
            fprintf('\t Error in line %d in function %s: %s\n', ...
                ME.stack(1).line, ME.stack(1).name, ME.message);
        end
    end
    
    try
        %% read transcranial volume in line (TIV_LINE_NUMBER-1) from FreeSurfer file 'aseg.stats'
        fid = fopen('../T1w/freesurfer/stats/aseg.stats');
        for ii = 1:(TIV_LINE_NUMBER-1) % skip unrelevant lines
            fgetl(fid);
        end
        tline = fgetl(fid);
        fclose(fid);

        tiv = regexp(tline, '\d+\.?\d*', 'match'); % read number (with or without decimal point)
        tiv = strrep(tiv, ',', ''); % delete comma
        
        % write tiv into tiv.txt
        fid = fopen('tiv_TEST.txt', 'w');
        fprintf(fid, '%s', tiv{1});
        fclose(fid);
    %% calculate mean FA
    
    % Binarize aparc+aseg.mgz and save as wm.mask.mgz
    system('mri_binarize --i ../T1w/freesurfer/mri/aparc+aseg.mgz --wm --o wm.mask.mgz');
 
    % Convert wm.mask.mgz to wm.nii.gz using nearest-neighbor interpolation
    system(['/opt/freesurfer/bin/mri_convert -rt nearest -rl ' resultsDir ...
        '/*_fractional_anisotropy.nii.gz wm.mask.mgz wm.nii.gz']);
 
    % Compute mean FA within wm.nii.gz mask and save as meanfa.txt
    system(['fslmeants -i ' resultsDir ...
        '/*_fractional_anisotropy.nii.gz -m wm.nii.gz -o meanfa.txt']);
 
    % Clean up
    system('rm -f wm*');
    catch ME
        fprintf('\t Error in line %d in function %s: %s\n', ...
            ME.stack(1).line, ME.stack(1).name, ME.message);
    end
    
    try 
        savefile='richclub_adj_TEST.mat'
        save(savefile,'rc')

        savefile=RESULTS_DESTINATION
        save(savefile, 'rcResults');

        savefile=DENSITY_DESTINATION
        save(savefile, 'rcDensity');
    catch ME
        fprintf('\t Error in line %d in function %s: %s\n', ...
            ME.stack(1).line, ME.stack(1).name, ME.message);
    end
end
