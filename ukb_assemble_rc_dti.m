%% load participants
dtiall=load('/slow/projects/01_UKB/00_scripts/ukb_rc_ids.txt');
n=length(dtiall);

%% Start parallel pool
UKB_Cluster = parcluster('local');
UKB_Cluster.NumWorkers = 100;
saveProfile(UKB_Cluster);
parpool('local',100);

%% pull rc data
edge={'fa' 'svd'}; %edge weights
parfor i=1:n  
    
    RCResults{i}.id=dtiall(i);
    
    try
        S=load(['/slow/projects/01_UKB/dti/00_batch1/' num2str(dtiall(i)) '/DWI_processed_v311/richclub_adj.mat'])
        for j=1:6
            for k=1:2
                RCResults{i}.max_phi(j,k)=S.rc{j}.(edge{k}).kmax(1);
                RCResults{i}.max_k(j,k)=S.rc{j}.(edge{k}).kmax(2);
                RCResults{i}.range(:,j,k)=S.rc{j}.(edge{k}).range;


                RCResults{i}.integral.norm(j,k) = trapz(find(~isnan(S.rc{j}.(edge{k}).norm)),S.rc{j}.(edge{k}).norm(~isnan(S.rc{j}.(edge{k}).norm)));
                RCResults{i}.integral.emp(j,k) = trapz(find(~isnan(S.rc{j}.(edge{k}).emp)),S.rc{j}.(edge{k}).emp(~isnan(S.rc{j}.(edge{k}).emp)));
                RCResults{i}.integral.rand(j,k) = trapz(find(~isnan(S.rc{j}.(edge{k}).rand)),S.rc{j}.(edge{k}).rand(~isnan(S.rc{j}.(edge{k}).rand)));
                RCResults{i}.integral.above(j,k)=RCResults{i}.integral.norm(j,k)-length(S.rc{j}.(edge{k}).norm(~isnan(S.rc{j}.(edge{k}).norm)));
                
                RCResults{i}.odd(1,j,k)=S.rc{j}.(edge{k}).Outliers(1);
                RCResults{i}.odd(2,j,k)=S.rc{j}.(edge{k}).RangeCorrected;
                RCResults{i}.odd(3,j,k)=S.rc{j}.(edge{k}).Ignored;
                RCResults{i}.odd(4,j,k)=length(S.rc{j}.(edge{k}).range);
            end
        end
    catch
        
    end
end
    savefile='/slow/projects/01_UKB/dti/rc_results_220622.mat';
    save(savefile,'RCResults')