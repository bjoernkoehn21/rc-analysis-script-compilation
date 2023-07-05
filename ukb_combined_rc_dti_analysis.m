%% Last step Changes
% change back to parfor and parallel pooling
% change number of random networks and parallel workers 2 blocks below to higher number

%% Preamble
% add paths to relevant packages
addpath /fast/home/marketts/tools/batches  % path to fdr_bh.m
addpath /home/marketts/imaging/2017_01_15_BCT % path to BCT toolbox
addpath /slow/projects/HCP_1200/00_scripts % path to rich_club_wu_extended
setenv('FREESURFER_HOME', '/opt/freesurfer');
system(fullfile(getenv('FREESURFER_HOME'), 'SetUpFreeSurfer.sh'));

% define destination of results and other relevant paths and file names
SUBJECTS_ID_FILE = 'ukb_subjIDs.txt';
PATH_TO_SUBJECT_DIRS = '/slow/projects/01_UKB/dti/00_batch1';
RESULTS_DESTINATION = ['/slow/projects/01_UKB/dti/rc_results_', ...
	datestr(now, 'yyyymmdd'), '.mat'];
% SUBJECTS_ID_FILE = 'hcp_subjIDs.txt';
% PATH_TO_SUBJECT_DIRS = '/slow/projects/HCP_1200/01_complete_batch';
% RESULTS_DESTINATION = ['/slow/projects/HCP_1200/rc_results_', ...
% 	datestr(now, 'yyyymmdd'), '.mat'];
% Define PATH_TO_ZIPPED_FREESURFER_FILES if files "aseg.stats" and "aparc+aseg.mgz" exist in
% zipped form or outside the PATH_TO_SUBJECT_DIRS path
PATH_TO_FREESURFER_DIR = '/slow/projects/01_UKB/surface/00_batch1';
%PATH_TO_FREESURFER_DIR = ''; 
flagCalculateCovariates = false;

% define constants
N_PARALLEL_WORKERS = 100;
N_RANDOM_NETWORKS = 2500;
EDGE_WEIGHTS = {'fa' 'svd'}; % fractional anisotropy, streamline volume density
FILES = {'*_connectivity_csd_dti_aparc.mat'
'*_connectivity_csd_dti_lausanne120.mat'
'*_connectivity_csd_dti_lausanne250.mat'
'*_connectivity_gqi_dti_aparc.mat'
'*_connectivity_gqi_dti_lausanne120.mat'
'*_connectivity_gqi_dti_lausanne250.mat'};
%FILES = {'*_connectivity_gqi_dti_lausanne250.mat'};
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
	'max_phi_numerator', nan(length(FILES), length(EDGE_WEIGHTS)), ...
	'max_phi_denominator', nan(length(FILES), length(EDGE_WEIGHTS)), ...
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
        FILES, rcResults{iSubj}, SIGNIFICANCE_LEVEL, RC_DOMINION_RATIO, flagCalculateCovariates);
    toc
end
delete(gcp('nocreate'));
savefile = RESULTS_DESTINATION;
save(savefile, 'rcResults');

%% subfunction
function [rcResults] = compute_correct_assemble_rc_dti(PATH_TO_SUBJECT_DIRS, PATH_TO_FREESURFER_DIR, subjID, N_RANDOM_NETWORKS, EDGE_WEIGHTS, ...
    FILES, rcResults, SIGNIFICANCE_LEVEL, RC_DOMINION_RATIO, flagCalculateCovariates)
    %% compute rc coefficient and p-values
    resultsDir = fullfile(PATH_TO_SUBJECT_DIRS, num2str(subjID), 'DWI_processed_v311');
    cd(resultsDir);
    covariates = struct(...
        'density', cell2struct( ...
            repmat({nan(length(FILES), 1)}, length(EDGE_WEIGHTS), 1), EDGE_WEIGHTS, 1), ...
        'tiv', nan(1,1), ...
        'meanFa', nan(1,1));
    

    
    for iFile = 1:length(FILES)
        try
			file = load(dir(FILES{iFile}).name); % load iFile
            
%             rc = cell(length(FILES), 1);
%             subfields = {'emp', 'svd', 'rand'};
%             subsubfields = {'phi', 'phiNum', 'phiDenom'};
%             rc(:) = {cell2struct(repmat({cell2struct(repmat({cell2struct(repmat({nan(1, 1)}, length(subsubfields), 1), ...
%         subsubfields, 1)}, length(subfields), 1), subfields, 1)}, length(EDGE_WEIGHTS), 1), EDGE_WEIGHTS, 1)};

            rc = cell(length(FILES), 1);
            subfields = {'emp', 'svd', 'rand'};
            subsubfields = {'phi', 'phiNum', 'phiDenom', 'Er'};
            rc(:) = {cell2struct(repmat({cell2struct(repmat({cell2struct(repmat({nan(1, 1)}, length(subsubfields), 1), ...
        subsubfields, 1)}, length(subfields), 1), subfields, 1)}, length(EDGE_WEIGHTS), 1), EDGE_WEIGHTS, 1)};
            
    
    
            rcTest = rc;
            [rcTest{iFile}.fa.emp.phi] = ...
                rich_club_wu(file.connectivity(:,:,3));
                
            [rcTest{iFile}.svd.emp.phi] = ...
                rich_club_wu(file.connectivity(:,:,13));
            
            % calculate rc coefficients for empirical network
%             [rc{iFile}.fa.emp.phi, rc{iFile}.fa.emp.phiNum, rc{iFile}.fa.emp.phiDenom] = ...
%                 rich_club_wu_extended(file.connectivity(:,:,3));
%                 
%             [rc{iFile}.svd.emp.phi, rc{iFile}.svd.emp.phiNum, rc{iFile}.svd.emp.phiDenom] = ...
%                 rich_club_wu_extended(file.connectivity(:,:,13));
            

			
			[rc{iFile}.fa.emp.phi, rc{iFile}.fa.emp.phiNum, rc{iFile}.fa.emp.phiDenom, rc{iFile}.fa.emp.Er] = ...
                rich_club_wu_extended(file.connectivity(:,:,3));
                
            [rc{iFile}.svd.emp.phi, rc{iFile}.svd.emp.phiNum, rc{iFile}.svd.emp.phiDenom, rc{iFile}.svd.emp.Er] = ...
                rich_club_wu_extended(file.connectivity(:,:,13));
            
				
            
            % calculate densities
			[covariats.density.fa(iFile), ~, ~] = density_und(file.connectivity(:,:,3));
            [covariats.density.svd(iFile), ~, ~] = density_und(file.connectivity(:,:,13));
            
            % preallocate (size not known beforehand [before executing rich_club_wu])
            [rc{iFile}.fa.rand.phi, rc{iFile}.fa.rand.phiNum, rc{iFile}.fa.rand.phiDenom, rc{iFile}.fa.rand.Er] = ...
                deal(nan(length(rc{iFile}.fa.emp.phi), N_RANDOM_NETWORKS));
            [rc{iFile}.svd.rand.phi, rc{iFile}.svd.rand.phiNum, rc{iFile}.svd.rand.phiDenom, rc{iFile}.svd.rand.Er] = ...
                deal(nan(length(rc{iFile}.svd.emp.phi), N_RANDOM_NETWORKS));
			
			% produce random networks for normalization
			for iRandomNetwork = 1:N_RANDOM_NETWORKS
				if rem(iRandomNetwork,500) == 0
					fprintf('subject %d, file %d, random network number %d\n' , ...
						subjID, iFile, iRandomNetwork)
                end
%                 [rc{iFile}.fa.rand.phi(:,iRandomNetwork), ...
%                     rc{iFile}.fa.rand.phiNum(:,iRandomNetwork), ...
%                     rc{iFile}.fa.rand.phiDenom(:,iRandomNetwork)] = ...
%                     rich_club_wu_extended(randmio_und(file.connectivity(:,:,3),10), ...
%                     length(rc{iFile}.fa.emp.phi));
%                 [rc{iFile}.svd.rand.phi(:,iRandomNetwork), ...
%                     rc{iFile}.svd.rand.phiNum(:,iRandomNetwork), ...
%                     rc{iFile}.svd.rand.phiDenom(:,iRandomNetwork)] = ...
%                     rich_club_wu_extended(randmio_und(file.connectivity(:,:,13),10), ...
%                     length(rc{iFile}.svd.emp.phi));
				
				[rc{iFile}.fa.rand.phi(:,iRandomNetwork), ...
                    rc{iFile}.fa.rand.phiNum(:,iRandomNetwork), ...
                    rc{iFile}.fa.rand.phiDenom(:,iRandomNetwork), rc{iFile}.fa.rand.Er(:,iRandomNetwork)] = ...
                    rich_club_wu_extended(randmio_und(file.connectivity(:,:,3), 10), ...
                    length(rc{iFile}.fa.emp.phi));
				%disp(rc{iFile}.fa.rand.Er(:,iRandomNetwork)')
                [rc{iFile}.svd.rand.phi(:,iRandomNetwork), ...
                    rc{iFile}.svd.rand.phiNum(:,iRandomNetwork), ...
                    rc{iFile}.svd.rand.phiDenom(:,iRandomNetwork), rc{iFile}.svd.rand.Er(:,iRandomNetwork)] = ...
                    rich_club_wu_extended(randmio_und(file.connectivity(:,:,13), 10), ...
                    length(rc{iFile}.svd.emp.phi));
				
% 				rc{iFile}.fa.rand(:,iRandomNetwork) = rich_club_wu( ...
% 					randmio_und(file.connectivity(:,:,3),10),length(rc{iFile}.fa.emp));
% 				rc{iFile}.svd.rand(:,iRandomNetwork) = rich_club_wu( ...
% 					randmio_und(file.connectivity(:,:,13),10),length(rc{iFile}.svd.emp));
			end
			%fprintf('all row elements are equal: %d\n', isequaln(rc{iFile}.svd.rand.phiDenom(1, :), rc{iFile}.svd.rand.phiDenom(1, 1)))
			%fprintf('row elements: %f\n', rc{iFile}.svd.rand.phiDenom(15, 1:10))
			%fprintf('number of connections: %d\n',  rc{iFile}.svd.rand.Er(15, 1:10))
			
			
			
			
% 			%disp(rcTest{iFile}.fa.emp.phi)
% 			%disp(rc{iFile}.fa.emp.phi)
%             %fprintf('emp is equal: %d\n', isequaln(rcTest{iFile}.fa.emp.phi, rc{iFile}.fa.emp.phi));
% 			%fprintf('size emp: %f\n', size(rc{iFile}.fa.emp.phiDenom))
% 			disp(size(rc{iFile}.fa.emp.phiDenom))
% 			%disp((rc{iFile}.fa.emp.phiDenom))
% 			disp(size(rc{iFile}.fa.emp.Er))
%             disp((rc{iFile}.fa.emp.Er))
% 			%fprintf('emp denom: %f\n', rc{iFile}.fa.emp.phiDenom);
% 			%fprintf('size rand: %f\n', size(rc{iFile}.fa.rand.phiDenom))
% 			disp(size(rc{iFile}.fa.rand.phiDenom))
% 			%disp((rc{iFile}.fa.rand.phiDenom(:,1)'))
% 			disp(size(rc{iFile}.fa.rand.Er(:,1)'))
%             disp((rc{iFile}.fa.rand.Er(:,1)'))
% 			disp((rc{iFile}.fa.rand.Er(:,2)'))
% 			%fprintf('denom is equal: %d\n', isequaln(rc{iFile}.fa.emp.phiDenom, rc{iFile}.fa.rand.phiDenom(:,1)'));
% 			%fprintf('rand first denom: %f\n', rc{iFile}.fa.rand.phiDenom(:,1));
        
			% normalize
			rc{iFile}.fa.norm = rc{iFile}.fa.emp.phi./mean(rc{iFile}.fa.rand.phi,2)';
			rc{iFile}.svd.norm = rc{iFile}.svd.emp.phi./mean(rc{iFile}.svd.rand.phi,2)';

			% calculate p-values
			rc{iFile}.fa.pvals = 1 - mean(rc{iFile}.fa.emp.phi > rc{iFile}.fa.rand.phi', 1);
			rc{iFile}.svd.pvals = 1 - mean(rc{iFile}.svd.emp.phi > rc{iFile}.svd.rand.phi', 1);
			
        catch ME
            fprintf('\t Error in line %d in function %s: %s\n', ...
                ME.stack(1).line, ME.stack(1).name, ME.message);
            rc{iFile} = [];
        end   
            
        %% correct rc range in cases where there is no unambiguous rc regime
        try
            for iEdgeWeight = 1:length(EDGE_WEIGHTS)
			
                % phi max and corresponding k  
                [rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).kmax(1), ...
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).kmax(2)] = max( ...
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
					
                rcResults.max_phi_numerator(iFile,iEdgeWeight) = rc{iFile}.( ...
                    EDGE_WEIGHTS{iEdgeWeight}).emp.phiNum(rcResults.max_k(iFile,iEdgeWeight));
                rcResults.max_phi_denominator(iFile,iEdgeWeight) = rc{iFile}.( ...
                    EDGE_WEIGHTS{iEdgeWeight}).emp.phiDenom(rcResults.max_k(iFile,iEdgeWeight));
				%fprintf('emp by hand: %f\n', rcResults.max_phi_numerator(iFile,iEdgeWeight)/rcResults.max_phi_denominator(iFile,iEdgeWeight))
				%fprintf('emp automatic: %f\n', mean(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).emp.phi(rcResults.max_k(iFile,iEdgeWeight))))

				
				%% calculate nominator and denominator of average of random network rc coefficients
% 				numerator = double(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand.phiNum(rcResults.max_k(iFile,iEdgeWeight), :));
% 				denominator = double(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand.phiDenom(rcResults.max_k(iFile,iEdgeWeight), :));
%                
% 				sortedUniqueDenom = sort(unique(denominator), 'descend');
% 				commonDenom = sortedUniqueDenom(1)*sortedUniqueDenom(2)
% 				disp(commonDenom)
% 				numeratorOfMean = sum(numerator .* (commonDenom ./ double(denominator)));
% 				denominatorOfMean = commonDenom * numel(denominator);
% 				
% 				fprintf('rand by hand: %f\n', numeratorOfMean/denominatorOfMean)
%  				fprintf('rand automatic: %f\n', mean(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand.phi(rcResults.max_k(iFile,iEdgeWeight), :)))
% 					
% 				numerator = double(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand.phiNum);
% 				denominator = double(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand.phiDenom);
% 
% 				% Calculate common denominator for each row
% 				SortedUniqueValuesPerRow = arrayfun(@(row) unique(denominator(row, :), 'stable'), 1:size(denominator, 1), 'UniformOutput', false);
% 				highestDenom = cellfun(@(values) values(end), SortedUniqueValuesPerRow);
% 				containsOneElementArray = any(cellfun(@(arr) numel(arr) == 1, SortedUniqueValuesPerRow));
% 				if containsOneElementArray
%                     secondHighestDenom = 1;
%                 else
%                     secondHighestDenom = cellfun(@(values) values(end-1), SortedUniqueValuesPerRow);
%                 end
% 				commonDenom = highestDenom .* secondHighestDenom;
% 
% 				% Calculate numerator and denominator of the mean for each row
% 				numeratorOfMean = sum(numerator .* (commonDenom' ./ denominator), 2);
% 				denominatorOfMean = commonDenom .* size(denominator, 2);
% 				
% 				meanRandCoeffByHand = numeratorOfMean./denominatorOfMean;
%  				meanRandCoeffAuto = mean(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand.phi, 2);
% 				size(meanRandCoeffByHand)
% 				size(meanRandCoeffAuto)
% 
%                 %% verify that meanRandCoeffByHand is equal to meanRandCoeffAuto
%                 % Define the tolerance
%                 tolerance = eps;
% 
%                 % Compare values within the tolerance
%                 isEqual = all(abs(meanRandCoeffByHand(~isnan(meanRandCoeffByHand)) - meanRandCoeffAuto(~isnan(meanRandCoeffAuto))) <= tolerance);
% 				%isequaln(meanRandCoeffByHand, meanRandCoeffAuto)
				
				
                numerator = double(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand.phiNum);
				denominator = double(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand.phiDenom);

				% Calculate common denominator for each row
				SortedUniqueValuesPerRow = arrayfun(@(row) unique(denominator(row, :), 'stable'), 1:size(denominator, 1), 'UniformOutput', false);
                containsOneElementArray = any(cellfun(@(arr) numel(arr) == 1, SortedUniqueValuesPerRow));
                highestDenom = cellfun(@(values) values(end), SortedUniqueValuesPerRow);
                if containsOneElementArray
                    secondHighestDenom = 1;
                else
                    secondHighestDenom = cellfun(@(values) values(end-1), SortedUniqueValuesPerRow);
                end
                
				commonDenom = highestDenom .* secondHighestDenom;

				% Calculate numerator and denominator of the mean for each row
				%fprintf('sizeNum: %d, ', size(numerator));
                %fprintf('sizeCommonDenum: %d, ', size(commonDenom));
				numeratorOfMean = sum(numerator .* (commonDenom' ./ denominator), 2);
				denominatorOfMean = commonDenom .* size(denominator, 2);
                %fprintf('sizeNumeratorOfMean: %d, ', size(numeratorOfMean));
                %fprintf('sizeDenominatorOfMean: %d\n, ', size(denominatorOfMean));
				
				meanRandCoeffByHand = numeratorOfMean./denominatorOfMean';
 				meanRandCoeffAuto = mean(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand.phi, 2);

                %% verify that meanRandCoeffByHand is equal to meanRandCoeffAuto
                % Define the tolerance
                tolerance = 1e-10;%eps;

                % Compare values within the tolerance
                isEqual = all(abs(meanRandCoeffByHand(~isnan(meanRandCoeffByHand)) - meanRandCoeffAuto(~isnan(meanRandCoeffAuto))) <= tolerance)
				%isequaln(meanRandCoeffByHand, meanRandCoeffAuto)
				
				
              
                rcResults.integral.norm(iFile,iEdgeWeight) = trapz( ...
                    find(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm)), ...
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm( ...
                    ~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm)));
                rcResults.integral.emp(iFile,iEdgeWeight) = trapz( ...
                    find(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).emp.phi)), ...
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).emp.phi( ...
                    ~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).emp.phi)));
                rcResults.integral.rand(iFile,iEdgeWeight) = trapz( ...
                    find(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand.phi)), ...
                    rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand.phi( ...
                    ~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).rand.phi)));
                rcResults.integral.above(iFile,iEdgeWeight) = rcResults.integral.norm( ...
					iFile,iEdgeWeight) - sum(~isnan(rc{iFile}.(EDGE_WEIGHTS{iEdgeWeight}).norm));
              
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
    if flagCalculateCovariates
	
		if isempty(PATH_TO_FREESURFER_DIR)
			pathToTiv = fullfile(PATH_TO_SUBJECT_DIRS, num2str(subjID), '**', 'aseg.stats');
		else
			pathToTiv = fullfile(PATH_TO_FREESURFER_DIR, num2str(subjID), '**', 'aseg.stats');
		end

		fileList = dir(pathToTiv);
		isFoundOnce = numel(fileList) == 1;
		isFoundMultipleTimes = numel(fileList) > 1;
		isZipped = false;
		flagAsegStatsFileFound = false;
		if isFoundOnce
			flagAsegStatsFileFound = true;
			pathToTiv = fullfile(fileList(1).folder, fileList(1).name);
		elseif isFoundMultipleTimes
			fprintf('Error: Multiple files "aseg.stats" found at %s', pathToTiv);
		else
			try
				pathToZip = [PATH_TO_FREESURFER_DIR, '/', num2str(subjID)];
				pathToTiv = fullfile(PATH_TO_FREESURFER_DIR, num2str(subjID), 'FreeSurfer/stats/aseg.stats');
				cd(pathToZip);
				system('unzip -qq *zip FreeSurfer/stats/aseg.stats')
				
				flagAsegStatsFileFound = true;
				isZipped = true;
			catch ME
				fprintf(['Error unzipping file at ', pathToZip, ...
					'/stats/aseg.stats. Please check the validity of PATH_TO_FREESURFER_DIR', ...
					': Error in line %d in function %s: %s\n'], ...
					ME.stack(1).line, ME.stack(1).name, ME.message);
				fprintf('File "aseg.stats" not found in subject dirctory %s', pathToZip, ...
					'\n Nor was it found at ', pathToZip, 'FreeSurfer/stats/aseg.stats', ...
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
					covariates.tiv = str2double(tiv{1});
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
		isFoundOnce = numel(fileList) == 1;
		isFoundMultipleTimes = numel(fileList) > 1;
		flagAparcAsegMgzFileFound = false;
		if isFoundOnce
			flagAparcAsegMgzFileFound = true;
			pathToFa = fullfile(fileList(1).folder, fileList(1).name);
		elseif isFoundMultipleTimes
			fprintf('Error: Multiple files "aparc+aseg.mgz" found in subject directory %s\n', pathToFa);
		else  
			try            
				%pathToZip = [PATH_TO_FREESURFER_DIR, '/', num2str(subjID)];
				%system(['unzip -qq ', pathToZip, '/*.zip /FreeSurfer/mri/aparc+aseg.mgz -d ', pathToZip, '/FreeSurfer/mri/aparc+aseg.mgz']);
				
				pathToZip = [PATH_TO_FREESURFER_DIR, '/', num2str(subjID)];
				pathToFa = fullfile(PATH_TO_FREESURFER_DIR, num2str(subjID), 'FreeSurfer/mri/aparc+aseg.mgz');
				cd(pathToZip);
				system('unzip -qq *zip FreeSurfer/mri/aparc+aseg.mgz')
				
				flagAparcAsegMgzFileFound = true;
				isZipped = true;
			catch ME
				fprintf(['Error unzipping file at ', pathToZip, ...
					'/stats/aseg.stats. Please check the validity of PATH_TO_FREESURFER_DIR', ': Error in line %d in function %s: %s\n'], ...
					ME.stack(1).line, ME.stack(1).name, ME.message);
				fprintf('File "aseg.stats" not found in subject dirctory ', pathToFa, ...
					'\n Nor was it found at ', pathToZip, 'FreeSurfer/mri/aparc+aseg.mgz', ...
					'\n in zipped form. Please correct variable PATH_TO_FREESURFER_DIR', ...
					' or set FLAG_CALCULATE_COVARIATES to false.');
			end
			
		end


		if flagAparcAsegMgzFileFound
			try
				system(['mri_binarize --i ', pathToFa, ' --wm --o wm.mask.mgz']);

				% Convert wm.mask.mgz to wm.nii.gz using nearest-neighbor interpolation
				system(['/opt/freesurfer/bin/mri_convert -rt nearest -rl ', resultsDir, ...
					'/*_fractional_anisotropy.nii.gz wm.mask.mgz wm.nii.gz']);

				% Compute mean FA within wm.nii.gz mask and save as meanfa.txt
				[~, meanFa] = system(['fslmeants -i ' resultsDir ...
					'/*_fractional_anisotropy.nii.gz -m wm.nii.gz']); % [status, output]
				covariates.meanFa = str2double(meanFa);
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
end
