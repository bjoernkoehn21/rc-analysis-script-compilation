%% TODO
% change to real ID-file: subjectIDs.txt
% activate parallel pools

% add paths
addpath /fast/home/marketts/tools/batches  % path to fdr_bh
addpath /home/marketts/imaging/2017_01_15_BCT % path to BCT toolbox
% Set up freesurfer environment
setenv('FREESURFER_HOME', '/opt/freesurfer');
system(fullfile(getenv('FREESURFER_HOME'), 'SetUpFreeSurfer.sh'));

% define constants
nParallelWorkers = 4; %100;
nRandomNetworks = 10; %2500;
tivLineNumber = 34; % this number refers to the line in SUBJECT/T1w/freesurfer/stats/aseg.stats that contains the total intrakranial volume
edge={'fa' 'svd'}; % edge weights: fractional anisotropy, 
files={'*_connectivity_csd_dti_aparc.mat'
'*_connectivity_csd_dti_lausanne120.mat'
'*_connectivity_csd_dti_lausanne250.mat'
'*_connectivity_gqi_dti_aparc.mat'
'*_connectivity_gqi_dti_lausanne120.mat'
'*_connectivity_gqi_dti_lausanne250.mat'};

% read in subjects IDs from subjIDs.txt
fid = fopen('subjIDs.txt', 'r');
subjIDs = textscan(fid, '%s');
subjIDs = subjIDs{1};
fclose(fid);

% assembled rc matrics
RCResults = cell(length(subjIDs), 1);
RCResults(:) = {struct( ...
    'id', zeros(1, 1), ...
    'max_phi', zeros(length(files), length(edge)), ...
    'max_k', zeros(length(files), length(edge)), ...
    'range', zeros(2, length(files), length(edge)), ...
    'integral', struct( ...
        'norm', zeros(length(files), length(edge)), ...
        'emp', zeros(length(files), length(edge)), ...
        'rand', zeros(length(files), length(edge)), ...
        'above', zeros(length(files), length(edge))), ...
    'odd', zeros(4, length(files), length(edge)))};

% density of empirical networks
RCdense = cell(length(subjIDs), 1);
RCdense(:) = {struct( ...
    'fa', zeros(length(files), 1), ...
    'svd', zeros(length(files), 1))};


%% Start parallel pool
UKB_Cluster = parcluster('local');
UKB_Cluster.NumWorkers = nParallelWorkers;
saveProfile(UKB_Cluster);
parpool('local', nParallelWorkers);

% loop over subject IDs and execute ukb_compute_rc_dti function in parallel
parfor iSubj = 1:length(subjIDs)
%for iSubj = 1:length(subjIDs)
    hcp_compute_correct_rc_dti(iSubj, subjIDs{iSubj}, nRandomNetworks, edge, files, RCResults, RCdense, tivLineNumber);
    
    savefile='/slow/projects/HCP_1200/rc_results_220622.mat';
    save(savefile,'RCResults')
    
    savefile='/slow/projects/HCP_1200/rcdense.mat'
    save(savefile,'RCdense')
end

%delete(gcp('nocreate')); % close interactive sessions

%% compute rc with preallocation for random networks and correct the rc
function hcp_compute_correct_rc_dti(iSubj, subjID, nRandomNetworks, edge, files, RCResults, RCdense, tivLineNumber)
    %resultsDir = fullfile('/slow/projects/01_UKB/dti/00_batch1', subjID, 'DWI_processed_v311');
    resultsDir = fullfile('/slow/projects/HCP_1200', subjID, 'DWI_processed_v311');
    cd(resultsDir);

%% compute rc coefficient and p-values
    tic

    for iFile=1:length(files) % loop over files
        try
            n=dir(files{iFile}); load(n.name) % load file
            
            % calculate rc coefficients for empirical network
            rc{iFile}.fa.emp=rich_club_wu(connectivity(:,:,3));
            rc{iFile}.svd.emp=rich_club_wu(connectivity(:,:,13));
            
            % calculate
            RCdense{iSubj}.fa(iFile)=density_und(connectivity(:,:,3));
            RCdense{iSubj}.svd(iFile)=density_und(connectivity(:,:,13));
            
            % preallocate (size of rc{iFile}.fa.emp not known beforehand [before executing rich_club_wu])
            rc{iFile}.fa.rand=zeros(length(rc{iFile}.fa.emp), nRandomNetworks);
            rc{iFile}.svd.rand=zeros(length(rc{iFile}.svd.emp), nRandomNetworks);
            rc{iFile}.fa.pvals=zeros(1,1);
            rc{iFile}.svd.pvals=zeros(1,1);
        catch ME
            fprintf('\t Error in line %d in function %s: %s\n', ME.stack(1).line, ME.stack(1).name, ME.message);
            rc{iFile}=[];
        end
    end

    
    for iFile = 1:length(files) % loop over files
            for iRandomNetwork = 1:nRandomNetworks
                rc{iFile}.fa.rand(:,iRandomNetwork)=rich_club_wu(randmio_und(connectivity(:,:,3),10),length(rc{iFile}.fa.emp));
                rc{iFile}.svd.rand(:,iRandomNetwork)=rich_club_wu(randmio_und(connectivity(:,:,13),10),length(rc{iFile}.svd.emp));
            end

            rc{iFile}.fa.norm=rc{iFile}.fa.emp./mean(rc{iFile}.fa.rand,2)';
            rc{iFile}.svd.norm=rc{iFile}.svd.emp./mean(rc{iFile}.svd.rand,2)';

            rc{iFile}.fa.pvals=1-sum(rc{iFile}.fa.emp>rc{iFile}.fa.rand')/nRandomNetworks;
            rc{iFile}.svd.pvals=1-sum(rc{iFile}.svd.emp>rc{iFile}.svd.rand')/nRandomNetworks;
            % mean seems to take longer
            %rc{iFile}.fa.pvals=1 - mean(rc{1}.fa.emp > rc{1}.fa.rand', 1);
            %rc{iFile}.svd.pvals=1 - mean(rc{1}.svd.emp > rc{1}.svd.rand', 1);
            
            
%% correct rc range in cases where there is no unambiguous rc regime
        try
            for iEdge = 1:length(edge)
              % phi max and corresponding k  
              [rc{iFile}.(edge{iEdge}).kmax(1,1),rc{iFile}.(edge{iEdge}).kmax(1,2)]=max(rc{iFile}.(edge{iEdge}).norm) ;

              % fdr-correction of p-values (adjusted p-values)
              %[~, ~,
              %rc{iFile}.(edge{iEdge}).adj_p]=fdr_bh(rc{iFile}.(edge{iEdge}).pvals,.05,'pdep','false');
              % original by Markett
               rc{iFile}.(edge{iEdge}).adj_p=fdr_bh(rc{iFile}.(edge{iEdge}).pvals,.05,'pdep','false');
              %%
              % rc range
                %x=1:length(rc{iFile}.(edge{iEdge}).adj_p); xp=x(rc{iFile}.(edge{iEdge}).adj_p<.05); %helper
                xp = find(rc{iFile}.(edge{iEdge}).adj_p < 0.05); % more memory-efficient than line above
                    rc{iFile}.(edge{iEdge}).range = [min(xp) max(xp)]; % rc range
                        idx =rc{iFile}.(edge{iEdge}).adj_p(min(xp):max(xp)); % identify values with p>=.05 within range
                        n = length(idx);
                        idx1 = (idx > 0.05); %idx1 contains positions in RC range where p > 0.05
                    rc{iFile}.(edge{iEdge}).Outliers = sum(idx1); rc{iFile}.(edge{iEdge}).RangeCorrected = 0;   rc{iFile}.(edge{iEdge}).Ignored =  0;

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
                               rc{iFile}.(edge{iEdge}).range= [min(x_cor_p) max(x_cor_p)];
                               rc{iFile}.(edge{iEdge}).RangeCorrected = 1; %Count of of subjects with corrections 
                           elseif b > (1.5*a)
                               idx_cor = zeros(1,n);
                               idx_cor(1,(max(x_RC_p)+1):n) = 1;
                               idx_cor = logical(idx_cor);
                               x_cor = min(xp):max(xp);
                               x_cor_p = x_cor(idx_cor);
                               rc{iFile}.(edge{iEdge}).range = [min(x_cor_p) max(x_cor_p)];
                               rc{iFile}.(edge{iEdge}).RangeCorrected = 1; %Count of of subjects with corrections
                           %If no RC range is 1.5 times bigger than other, select RC
                           %range with higher K values
                           else
                               idx_cor = zeros(1,n);
                               idx_cor(1,(max(x_RC_p)+1):n) = 1;
                               idx_cor = logical(idx_cor);
                               x_cor = min(xp):max(xp);
                               x_cor_p = x_cor(idx_cor);
                                rc{iFile}.(edge{iEdge}).range = [min(x_cor_p) max(x_cor_p)];
                               rc{iFile}.(edge{iEdge}).RangeCorrected = 1; 
                           end

                        elseif (length(idx2) - sum(idx2)) == 2 %If only two not
                            %continous p-values are above 0.05, ignore them
                            rc{iFile}.(edge{iEdge}).Ignored =  1;
                        else %If range of p-values above 0.05 is not continous, exclude subjects from analysis
                            rc{iFile}.(edge{iEdge}).range = [NaN NaN];
                        end
                    elseif sum(idx1) == 1
                        rc{iFile}.(edge{iEdge}).Outliers(2)=1; %Count number of subjects with single outliers
                    end
                    
                    % assemble relevant metrics
                    RCResults{iSubj}.id=subjID;
                    RCResults{iSubj}.max_phi(iFile,iEdge)=rc{iFile}.(edge{iEdge}).kmax(1);
                    RCResults{iSubj}.max_k(iFile,iEdge)=rc{iFile}.(edge{iEdge}).kmax(2);
                    RCResults{iSubj}.range(:,iFile,iEdge)=rc{iFile}.(edge{iEdge}).range;


                    RCResults{iSubj}.integral.norm(iFile,iEdge) = trapz(find(~isnan(rc{iFile}.(edge{iEdge}).norm)),rc{iFile}.(edge{iEdge}).norm(~isnan(rc{iFile}.(edge{iEdge}).norm)));
                    RCResults{iSubj}.integral.emp(iFile,iEdge) = trapz(find(~isnan(rc{iFile}.(edge{iEdge}).emp)),rc{iFile}.(edge{iEdge}).emp(~isnan(rc{iFile}.(edge{iEdge}).emp)));
                    RCResults{iSubj}.integral.rand(iFile,iEdge) = trapz(find(~isnan(rc{iFile}.(edge{iEdge}).rand)),rc{iFile}.(edge{iEdge}).rand(~isnan(rc{iFile}.(edge{iEdge}).rand)));
                    RCResults{iSubj}.integral.above(iFile,iEdge) = RCResults{iSubj}.integral.norm(iFile,iEdge)-length(rc{iFile}.(edge{iEdge}).norm(~isnan(rc{iFile}.(edge{iEdge}).norm)));
                    
                    RCResults{iSubj}.odd(1,iFile,iEdge)=rc{iFile}.(edge{iEdge}).Outliers(1);
                    RCResults{iSubj}.odd(2,iFile,iEdge)=rc{iFile}.(edge{iEdge}).RangeCorrected;
                    RCResults{iSubj}.odd(3,iFile,iEdge)=rc{iFile}.(edge{iEdge}).Ignored;
                    RCResults{iSubj}.odd(4,iFile,iEdge)=length(rc{iFile}.(edge{iEdge}).range); 
            end
        catch ME
            fprintf('\t Error in line %d in function %s: %s\n', ME.stack(1).line, ME.stack(1).name, ME.message);
        end
    end
    
    try
%% read transkranial volume from FreeSurfer file 'stats/aseg.stats'
        %resultsDir = fullfile('C:\Users\askou\Desktop\Geschaeftliches\HU_Markett-Aktueller_Job\HCP-data_download_and_CATO\MATLAB_analysis', num2str(subjID), 'DWI_processed_v311');
        fid = fopen('../T1w/freesurfer/stats/aseg.stats');
        
        for ii = 1:(tivLineNumber-1) % ignore lines 1 to (tivLineNumber-1) (faster and without as much storing than/as with textscan)
            fgetl(fid);
        end

        tline = fgetl(fid); %read line tivLineNumber with relevant tiv information
        fclose(fid);

        tiv = regexp(tline, '\d+\.?\d*', 'match'); % read number (with or without decimal point) from line 35 
        tiv = strrep(tiv, ',', ''); % delete comma
        % write number of tiv into tiv.txt
        fid = fopen('tiv.txt', 'w');
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

    savefile='richclub_adj.mat'; % adjusted richclub
    save(savefile,'rc')
    toc
end

