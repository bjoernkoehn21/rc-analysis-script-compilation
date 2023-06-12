function ukb_correct_rc_dti_TEST(subjID)
cd(['/slow/projects/01_UKB/dti/00_batch1/', num2str(subjID), '/DWI_processed_v311']);
%%load /slow/projects/01_UKB/dti/richclubcurves.mat
load(['/slow/projects/01_UKB/dti/richclubcurves_', num2str(subjID), '.mat'])
addpath /fast/home/marketts/tools/batches  
edge={'fa' 'svd'}; %edge weights
for i=1:6
    for j=1:2
  % phi max and corresponding k  
  [rc{i}.(edge{j}).kmax(1,1),rc{i}.(edge{j}).kmax(1,2)]=max(rc{i}.(edge{j}).norm) ;
 
  % fdr-correction of p-values
  [~, ~, rc{i}.(edge{j}).adj_p]=fdr_bh(rc{i}.(edge{j}).pvals,.05,'pdep','no');

  %%
  % rc range
    x=1:length(rc{i}.(edge{j}).adj_p);xp=x(rc{i}.(edge{j}).adj_p<.05); %helper
        rc{i}.(edge{j}).range= [min(xp) max(xp)]; % rc range
            idx =rc{i}.(edge{j}).adj_p(min(xp):max(xp)); % identify values with p>=.05 within range
            n = length(idx);
            idx1 = (idx > 0.05); %idx1 contains positions in RC range where p > 0.05
        rc{i}.(edge{j}).Outliers = sum(idx1); rc{i}.(edge{j}).RangeCorrected = 0;   rc{i}.(edge{j}).Ignored =  0;
        
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
                   rc{i}.(edge{j}).range= [min(x_cor_p) max(x_cor_p)];
                   rc{i}.(edge{j}).RangeCorrected = 1; %Count of of subjects with corrections 
               elseif b > (1.5*a)
                   idx_cor = zeros(1,n);
                   idx_cor(1,(max(x_RC_p)+1):n) = 1;
                   idx_cor = logical(idx_cor);
                   x_cor = min(xp):max(xp);
                   x_cor_p = x_cor(idx_cor);
                   rc{i}.(edge{j}).range = [min(x_cor_p) max(x_cor_p)];
                   rc{i}.(edge{j}).RangeCorrected = 1; %Count of of subjects with corrections
               %If no RC range is 1.5 times bigger than other, select RC
               %range with higher K values
               else
                   idx_cor = zeros(1,n);
                   idx_cor(1,(max(x_RC_p)+1):n) = 1;
                   idx_cor = logical(idx_cor);
                   x_cor = min(xp):max(xp);
                   x_cor_p = x_cor(idx_cor);
                    rc{i}.(edge{j}).range = [min(x_cor_p) max(x_cor_p)];
                   rc{i}.(edge{j}).RangeCorrected = 1; 
               end
               
            elseif (length(idx2) - sum(idx2)) == 2 %If only two not
                %continous p-values are above 0.05, ignore them
                rc{i}.(edge{j}).Ignored =  1;
            else %If range of p-values above 0.05 is not continous, exclude subjects from analysis
                
                rc{i}.(edge{j}).range = [NaN NaN];
            end
        elseif sum(idx1) == 1
            rc{i}.(edge{j}).Outliers(2)=1; %Count number of subjects with single outliers
        end
    end
end

%savefile='richclub_adj.mat'
savefile=['/slow/projects/01_UKB/dti/richclub_adj_', num2str(subjID), '.mat']
save(savefile,'rc')
end
