%% Last step Changes
% change back to parfor and parallel pooling
% change number of random networks and parallel workers 2 blocks below to higher number

%% Preamble
% add paths to relevant packages
addpath /fast/home/marketts/tools/batches  % path to fdr_bh.m
addpath /home/marketts/imaging/2017_01_15_BCT % path to BCT toolbox
% Set up freesurfer environment
setenv('FREESURFER_HOME', '/opt/freesurfer');
system(fullfile(getenv('FREESURFER_HOME'), 'SetUpFreeSurfer.sh'));

% define destination of results and other relevant paths and file names
SUBJECTS_ID_FILE = 'subjIDs.txt';
% PATH_TO_SUBJECT_DIRS = '/slow/projects/01_UKB/dti/00_batch1';
% RESULTS_DESTINATION = ['/slow/projects/01_UKB/dti/rc_results_', ...
% 	datestr(now, 'yyyymmdd'), '.mat'];
PATH_TO_SUBJECT_DIRS = '/slow/projects/HCP_1200/01_complete_batch';
RESULTS_DESTINATION = ['/slow/projects/HCP_1200/rc_results_', ...
	datestr(now, 'yyyymmdd'), '.mat'];
% Define PATH_TO_ZIPPED_FREESURFER_FILES if files "aseg.stats" and "aparc+aseg.mgz" exist in
% zipped form or outside the PATH_TO_SUBJECT_DIRS path
PATH_TO_FREESURFER_DIR = '';  

% define constants
N_PARALLEL_WORKERS = 100;
N_RANDOM_NETWORKS = 500;%2500;
EDGE_WEIGHTS = {'fa' 'svd'}; % fractional anisotropy, streamline volume density
% FILES = {'*_connectivity_csd_dti_aparc.mat'
% '*_connectivity_csd_dti_lausanne120.mat'
% '*_connectivity_csd_dti_lausanne250.mat'
% '*_connectivity_gqi_dti_aparc.mat'
% '*_connectivity_gqi_dti_lausanne120.mat'
% '*_connectivity_gqi_dti_lausanne250.mat'};
FILES = {'*_connectivity_gqi_dti_lausanne250.mat'};
RC_DOMINION_RATIO = 1.5;
SIGNIFICANCE_LEVEL = 0.5;

% read in subjects IDs from subjIDs.txt
fid = fopen(SUBJECTS_ID_FILE, 'r');
subjIDs = textscan(fid, '%s');
subjIDs = subjIDs{1};
fclose(fid);

% preallocate space for assembled-rc-results matrics before loop for time reasons
rcResults = cell(length(subjIDs), 1);
rcResults(:) = {struct( ...
    'id', nan(1, 1), ...
    'max_phi', nan(length(FILES), length(EDGE_WEIGHTS)), ...
    'max_k', nan(length(FILES), length(EDGE_WEIGHTS)), ...
    'range', nan(2, length(FILES), length(EDGE_WEIGHTS)), ...
    'integral', struct( ...
        'norm', nan(length(FILES), length(EDGE_WEIGHTS)), ...
        'emp', nan(length(FILES), length(EDGE_WEIGHTS)), ...
        'rand', nan(length(FILES), length(EDGE_WEIGHTS)), ...
        'above', nan(length(FILES), length(EDGE_WEIGHTS))), ...
    'odd', nan(4, length(FILES), length(EDGE_WEIGHTS)))};

%% Main process 

% Start parallel pool
% UKB_Cluster = parcluster('local');
% UKB_Cluster.NumWorkers = N_PARALLEL_WORKERS;
% saveProfile(UKB_Cluster);
% parpool('local', N_PARALLEL_WORKERS);

% loop over subject IDs and execute compute_correct_assemble_rc_dti function in parallel
%parfor iSubj = 1:length(subjIDs)
for iSubj = 1:length(subjIDs)
    tic
    rcResults{iSubj} = compute_correct_assemble_rc_dti(PATH_TO_SUBJECT_DIRS, PATH_TO_FREESURFER_DIR, str2double(subjIDs{iSubj}), N_RANDOM_NETWORKS, EDGE_WEIGHTS, ...
        FILES, rcResults{iSubj}, SIGNIFICANCE_LEVEL, RC_DOMINION_RATIO);
    toc
end
delete(gcp('nocreate'));
savefile = RESULTS_DESTINATION;
save(savefile, 'rcResults');

%% subfunction
function [rcResults] = compute_correct_assemble_rc_dti(PATH_TO_SUBJECT_DIRS, PATH_TO_FREESURFER_DIR, subjID, N_RANDOM_NETWORKS, EDGE_WEIGHTS, ...
    FILES, rcResults, SIGNIFICANCE_LEVEL, RC_DOMINION_RATIO)
    %% compute rc coefficient and p-values
    resultsDir = fullfile(PATH_TO_SUBJECT_DIRS, num2str(subjID), 'DWI_processed_v311');
    cd(resultsDir);
    covariates = struct(...
        'density', struct( ...
            'fa', nan(length(FILES), 1), ...
            'svd', nan(length(FILES), 1), ...
        'tiv', nan(1,1), ...
        'meanFa', nan(1,1)));

    covariates = struct(...
        'density', cell2struct( ...
            repmat({nan(length(FILES), 1)}, length(EDGE_WEIGHTS), 1), EDGE_WEIGHTS, 1), ...
        'tiv', nan(1,1), ...
        'meanFa', nan(1,1));
    

    
    for iFile = 1:length(FILES)
        try
			file = load(dir(FILES{iFile}).name); % load iFile
            
            % calculate rc coefficients for empirical network
            rc{iFile}.fa.emp=rich_club_wu(file.connectivity(:,:,3));
            rc{iFile}.svd.emp=rich_club_wu(file.connectivity(:,:,13));
            
            % calculate densities
			[covariats.density.fa(iFile), ~, ~] = density_und(file.connectivity(:,:,3));
            [covariats.density.svd(iFile), ~, ~] = density_und(file.connectivity(:,:,13));
            
            % preallocate (size not known beforehand [before executing rich_club_wu])
            rc{iFile}.fa.rand=zeros(length(rc{iFile}.fa.emp), N_RANDOM_NETWORKS);
            rc{iFile}.svd.rand=zeros(length(rc{iFile}.svd.emp), N_RANDOM_NETWORKS);
			
			% produce random networks for normalization
			for iRandomNetwork = 1:N_RANDOM_NETWORKS
				if rem(iRandomNetwork,500) == 0
					fprintf('subject %d, file %d, random network number %d\n' , ...
						subjID, iFile, iRandomNetwork)
				end
				rc{iFile}.fa.rand(:,iRandomNetwork) = rich_club_wu( ...
					randmio_und(file.connectivity(:,:,3),10),length(rc{iFile}.fa.emp));
				rc{iFile}.svd.rand(:,iRandomNetwork) = rich_club_wu( ...
					randmio_und(file.connectivity(:,:,13),10),length(rc{iFile}.svd.emp));
			end
        
			% normalize
			rc{iFile}.fa.norm = rc{iFile}.fa.emp./mean(rc{iFile}.fa.rand,2)';
			rc{iFile}.svd.norm = rc{iFile}.svd.emp./mean(rc{iFile}.svd.rand,2)';

			% calculate p-values
			rc{iFile}.fa.pvals = 1 - mean(rc{iFile}.fa.emp > rc{iFile}.fa.rand', 1);
			rc{iFile}.svd.pvals = 1 - mean(rc{iFile}.svd.emp > rc{iFile}.svd.rand', 1);
			
        catch ME
            fprintf('\t Error in line %d in function %s: %s\n', ...
                ME.stack(1).line, ME.stack(1).name, ME.message);
            rc{iFile} = [];
        end   
            
        %% correct rc range in cases where there is no unambiguous rc regime
        try
            for iEdgeWeight = 1:length(EDGE_WEIGHTS)
			
                % phi max and corresponding k  
                [rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).kmax(1,1), ...
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).kmax(1,2)] = max( ...
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm) ;
                
                % fdr-correction of p-values (adjusted p-values)
                [~, ~, rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).adj_p] = fdr_bh( ...
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).pvals, SIGNIFICANCE_LEVEL, 'pdep', 'no');
    
                rcIndexes = find(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).adj_p < SIGNIFICANCE_LEVEL);
                rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range = [min(rcIndexes) max(rcIndexes)];
				if isempty(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range)
					rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range = [NaN NaN];
				end
              
                % identify values with p>=0.05 within range
                rcRangePvalues = rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).adj_p( ...
                    min(rcIndexes):max(rcIndexes));
                rcRangeLength = length(rcRangePvalues);
                isOutlierInRcRange = (rcRangePvalues > SIGNIFICANCE_LEVEL);
			  
              
                rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Outliers = sum(isOutlierInRcRange);
                rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).RangeCorrected = false;
                rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Ignored =  false;
                    
                if sum(isOutlierInRcRange) > 1 % if there's more than 1 outlier
                  
					% Range of RC range between minimum and maximum K in RC range > 0.05
					outliersIndexesInRcRange = find(isOutlierInRcRange);
					pvaluesFirstOutlierToLastOutlier = rcRangePvalues( ...
						min(outliersIndexesInRcRange):max(outliersIndexesInRcRange));
					nNonOutliersInOutlierRange = sum( ...
						pvaluesFirstOutlierToLastOutlier < SIGNIFICANCE_LEVEL); 
                  
                    % If range of p-values above 0.05 is continous, select bigger RC-range and  
                    % define it as adjusted RC range
                    if nNonOutliersInOutlierRange <= 1
                        nRcEntriesBeforeFirstOutlier = length(1:(min(outliersIndexesInRcRange)-1)); 
                        nRcEntriesAfterLastOutlier = length( ...
                            (max(outliersIndexesInRcRange)+1):rcRangeLength);
						  
					    % most rc entries before outliers
                        if nRcEntriesBeforeFirstOutlier > ( ...
							RC_DOMINION_RATIO * nRcEntriesAfterLastOutlier)
						    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range= [ ...
                                min(rcIndexes) (min(rcIndexes) + nRcEntriesBeforeFirstOutlier -1)];
							  
					    % most rc entries after outliers or no RC range is 1.5 times bigger 
                        else 
						    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range= [ ...
                                (max(rcIndexes) - nRcEntriesAfterLastOutlier + 1) max(rcIndexes)];
                        end
					    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).RangeCorrected = true;
					  
					    % If there is only two outliers they are ignored 
                    elseif sum(isOutlierInRcRange) == 2 
                        rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Ignored =  true;
					  
                    % If range of p-vals above .05 is not continous, exclude subjects from analysis
                    else
                        rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range = [NaN NaN];
                    end
				elseif sum(isOutlierInRcRange) == 1
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).Outliers(2) = true;
                end
              
                %% assemble relevant metrics
                rcResults.id = subjID;
                rcResults.max_phi(iFile,iEdgeWeight) = rc{iFile}.( ...
                    EDGE_WEIGHTS{iEdgeWeight}).kmax(1);
                rcResults.max_k(iFile,iEdgeWeight) = rc{iFile}.( ...
                    EDGE_WEIGHTS{iEdgeWeight}).kmax(2);
                rcResults.range(:,iFile,iEdgeWeight) = rc{iFile}.( ...
                    EDGE_WEIGHTS{iEdgeWeight}).range;
              
                rcResults.integral.norm(iFile,iEdgeWeight) = trapz( ...
                    find(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm)), ...
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm( ...
                    ~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm)));
                rcResults.integral.emp(iFile,iEdgeWeight) = trapz( ...
                    find(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).emp)), ...
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).emp( ...
                    ~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).emp)));
                rcResults.integral.rand(iFile,iEdgeWeight) = trapz( ...
                    find(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand)), ...
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand( ...
                    ~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand)));
                rcResults.integral.above(iFile,iEdgeWeight) = rcResults.integral.norm( ...
					iFile,iEdgeWeight) - length(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm( ...
                    ~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm)));
              
                rcResults.odd(1,iFile,iEdgeWeight) = rc{iFile}.( ...
                    EDGE_WEIGHTS{iEdgeWeight}).Outliers(1);
                rcResults.odd(2,iFile,iEdgeWeight) = rc{iFile}.( ...
                    EDGE_WEIGHTS{iEdgeWeight}).RangeCorrected;
                rcResults.odd(3,iFile,iEdgeWeight) = rc{iFile}.( ...
                    EDGE_WEIGHTS{iEdgeWeight}).Ignored;
                rcResults.odd(4,iFile,iEdgeWeight) = ...
					rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range(2) - ...
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).range(1);
            end
        catch ME
            fprintf('\t Error in line %d in function %s: %s\n', ...
                ME.stack(1).line, ME.stack(1).name, ME.message);
        end
    end


    try
	savefile = ['/slow/projects/01_UKB/dti/richclub_adj_', num2str(subjID), '.mat']; %% This matrix should be saved in DWI_processed_v311. But currently the access rights prohibit this.
	save(savefile,'rc')
    catch ME
	fprintf('\t Error in line %d in function %s: %s\n', ...
	    ME.stack(1).line, ME.stack(1).name, ME.message);
    end

	%% read total intracranial volume from FreeSurfer file 'aseg.stats'
      
    if isempty(PATH_TO_FREESURFER_DIR)
        pathToTiv = fullfile(PATH_TO_SUBJECT_DIRS, num2str(subjID), '**', 'aseg.stats');
    else
        pathToTiv = fullfile(PATH_TO_FREESURFER_DIR, num2str(subjID), '**', 'aseg.stats');
    end

    fileList = dir(pathToTiv);
    pathToTiv = fullfile(fileList(1).folder, fileList(1).name);
	isFoundOnce = numel(fileList) == 1;
	isFoundMultipleTimes = numel(fileList) > 1;
	isZipped = false;
	if isFoundOnce
		flagAsegStatsFileFound = true;
	elseif isFoundMultipleTimes
        flagAsegStatsFileFound = false;
        fprintf('Error: Multiple files "aseg.stats" found at %s', pathToTiv);
	else
        try
            system(['unzip -qq *zip ', PATH_TO_FREESURFER_DIR, '/stats/aseg.stats']);
            pathToTiv = fullfile(PATH_TO_FREESURFER_DIR, '/stats/aseg.stats');
			flagAsegStatsFileFound = true;
			isZipped = true;
        catch ME
            fprintf(['Error unzipping file at ', PATH_TO_FREESURFER_DIR, ...
                '/stats/aseg.stats. Please check the validity of PATH_TO_FREESURFER_DIR', ...
				': Error in line %d in function %s: %s\n'], ...
                ME.stack(1).line, ME.stack(1).name, ME.message);
			fprintf('File "aseg.stats" not found in subject dirctory %s', pathToTiv, ...
				'\n Nor was it found at ', PATH_TO_FREESURFER_DIR, '/stats/aseg.stats', ...
				'\n in zipped form. Please correct variable PATH_TO_FREESURFER_DIR', ...
				' or set FLAG_CALCULATE_COVARIATES to false.');
        end
    end


    if flagAsegStatsFileFound
        fid = fopen(pathToTiv);
        tline = '';
        while ischar(tline)
            tline = fgetl(fid);
            if contains(tline, 'EstimatedTotalIntraCranialVol')
                tiv = regexp(tline, '\d+\.\d+', 'match');
                covariates.tiv = str2num(tiv{1});
                break;
            end
        end
        fclose(fid);
		if isZipped
			rmdir('FreeSurfer', 's'); % Remove the FreeSurfer directory and files
		end
    end
  
		
    %% calculate mean FA
    if isempty(PATH_TO_FREESURFER_DIR)
        pathToFa = fullfile(PATH_TO_SUBJECT_DIRS, num2str(subjID), '**', 'aparc+aseg.mgz');
    else
        pathToFa = fullfile(PATH_TO_FREESURFER_DIR, num2str(subjID), '**', 'aparc+aseg.mgz');
    end

    fileList = dir(pathToFa);
    pathToFa = fullfile(fileList(1).folder, fileList(1).name);
	isFoundOnce = numel(fileList) == 1;
	isFoundMultipleTimes = numel(fileList) > 1;
	if isFoundOnce
		flagAparcAsegMgzFileFound = true;
    elseif isFoundMultipleTimes
        flagAparcAsegMgzFileFound = false;
        fprintf('Error: Multiple files "aparc+aseg.mgz" found in subject directory %s', pathToFa);
    else  
        try
            system(['unzip -qq *zip ', PATH_TO_FREESURFER_DIR, '/mri/aparc+aseg.mgz']);
            pathToFa = fullfile(PATH_TO_FREESURFER_DIR, '/mri/aparc+aseg.mgz');
			flagAparcAsegMgzFileFound = true;
			isZipped = true;
        catch ME
            fprintf(['Error unzipping file at ', PATH_TO_FREESURFER_DIR, ...
                '/stats/aseg.stats. Please check the validity of PATH_TO_FREESURFER_DIR', ': Error in line %d in function %s: %s\n'], ...
                ME.stack(1).line, ME.stack(1).name, ME.message);
			fprintf('File "aseg.stats" not found in subject dirctory ', pathToFa, ...
				'\n Nor was it found at ', PATH_TO_FREESURFER_DIR, '/mri/aparc+aseg.mgz', ...
				'\n in zipped form. Please correct variable PATH_TO_FREESURFER_DIR', ...
				' or set FLAG_CALCULATE_COVARIATES to false.');
        end
        
    end


    if flagAparcAsegMgzFileFound
        try
            system(['mri_binarize --i ', pathToFa, ' --wm --o wm.mask.mgz']);

            % Convert wm.mask.mgz to wm.nii.gz using nearest-neighbor interpolation
            system(['/opt/freesurfer/bin/mri_convert -rt nearest -rl ' resultsDir...
                '/*_fractional_anisotropy.nii.gz wm.mask.mgz wm.nii.gz']);

            % Compute mean FA within wm.nii.gz mask and save as meanfa.txt
            [~, meanFa] = system(['fslmeants -i ' resultsDir ...
                '/*_fractional_anisotropy.nii.gz -m wm.nii.gz']); % [status, output]
			covariates.meanFa = str2num(meanFa);
			if isZipped
				rmdir('FreeSurfer', 's'); % Remove the FreeSurfer directory and files
			end
			delete('wm*');
        catch ME
            fprintf('\t Error in line %d in function %s: %s\n', ...
                ME.stack(1).line, ME.stack(1).name, ME.message);
        end
    end
    savefile = fullfile(PATH_TO_SUBJECT_DIRS, num2str(subjID), 'covariates.mat')
	save(savefile, 'covariates');
end
