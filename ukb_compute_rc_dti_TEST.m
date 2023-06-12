function ukb_compute_rc_dti_TEST(subjID)
tic
seed = 1;
files={'*_connectivity_csd_dti_aparc.mat'
'*_connectivity_csd_dti_lausanne120.mat'
'*_connectivity_csd_dti_lausanne250.mat'
'*_connectivity_gqi_dti_aparc.mat'
'*_connectivity_gqi_dti_lausanne120.mat'
'*_connectivity_gqi_dti_lausanne250.mat'}

% add BCT to path
addpath /home/marketts/imaging/2017_01_15_BCT
cd(['/slow/projects/01_UKB/dti/00_batch1/', num2str(subjID), '/DWI_processed_v311'])
for i=1:length(files) % loop over files
    fprintf('------------------------ subject %d file %d ------------------------\n------------------------------------------------------------------------\n', subjID, i)
    try
    n=dir(files{i}); load(n.name) % load file
    
    rc{i}.fa.emp=rich_club_wu(connectivity(:,:,3));
    rc{i}.svd.emp=rich_club_wu(connectivity(:,:,13));
    
    for j=1:2500
	rng(seed, 'twister')
        rc{i}.fa.rand(:,j)=rich_club_wu(randmio_und(connectivity(:,:,3),10),length(rc{i}.fa.emp));
	rng(seed, 'twister')
        rc{i}.svd.rand(:,j)=rich_club_wu(randmio_und(connectivity(:,:,13),10),length(rc{i}.svd.emp));
	if rem(j, 1000) == 0
		fprintf('subject: %d, file: %d, random network: %d\n', subjID, i, j)
	end
    end
    
    rc{i}.fa.norm=rc{i}.fa.emp./mean(rc{i}.fa.rand,2)';
    rc{i}.svd.norm=rc{i}.svd.emp./mean(rc{i}.svd.rand,2)';
    
    rc{i}.fa.pvals=1-sum(rc{i}.fa.emp>rc{i}.fa.rand')/2500;
    rc{i}.svd.pvals=1-sum(rc{i}.svd.emp>rc{i}.svd.rand')/2500;
    catch
         rc{i}=[];
    end
end

%savefile='richclubcurves.mat'
savefile=['/slow/projects/01_UKB/dti/richclubcurves_', num2str(subjID), '.mat']
save(savefile,'rc')
fprintf('saved file %s', savefile)
toc
end
