%% load participants
dtiall=load('/slow/projects/01_UKB/00_scripts/ukb_subjIDs.txt');
n=length(dtiall);

files={'*_connectivity_csd_dti_aparc.mat'
'*_connectivity_csd_dti_lausanne120.mat'
'*_connectivity_csd_dti_lausanne250.mat'
'*_connectivity_gqi_dti_aparc.mat'
'*_connectivity_gqi_dti_lausanne120.mat'
'*_connectivity_gqi_dti_lausanne250.mat'}

% add BCT to path
addpath /home/marketts/imaging/2017_01_15_BCT

%% Start parallel pool
delete(gcp('nocreate'))
UKB_Cluster = parcluster('local');
UKB_Cluster.NumWorkers = 100;
saveProfile(UKB_Cluster);
parpool('local',100);

%% pull rc data
edge={'fa' 'svd'}; %edge weights
parfor j=1:n  
    
    RCdense{j}.id=dtiall(j);
    
   % try
        cd(['/slow/projects/01_UKB/dti/00_batch1/' num2str(dtiall(j)) '/DWI_processed_v311/'])
    
            for i=1:length(files) % loop over files
                
                m=dir(files{i}); o=load(m.name); % load file

                RCdense{j}.fa(i)=density_und(o.connectivity(:,:,3));
                RCdense{j}.svd(i)=density_und(o.connectivity(:,:,13));
             
            end
    %catch
    %end
end
savefile='/slow/projects/01_UKB/dti/rcdense_1000010_1000028.mat'
RCdense_1000010_1000028 = RCdense
save(savefile,'RCdense_1000010_1000028')
