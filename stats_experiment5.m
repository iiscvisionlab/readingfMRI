%% 
allclear
load('L2_EMTbigram_searches.mat')
qt_sub = L2_str.subjinfo.ismalayalam==0;
qm_sub = L2_str.subjinfo.ismalayalam==1;

RTt = rmRToutlier(L2_str.RT(:,qt_sub,:),2,2);
RTm = rmRToutlier(L2_str.RT(:,qm_sub,:),2,2);

RT_mal = nanmean(RTm,3);
RT_tel = nanmean(RTt,3);

load('pair_RT_2baton.mat');
y_pred_t = Batonmodel(1./nanmean(RT_tel,2), 5, 2,3,0,[],[],1./rt_tt); rt = nancorrcoef(y_pred_t,1./nanmean(RT_tel,2));
y_pred_m = Batonmodel(1./nanmean(RT_mal,2), 5, 2,3,0,[],[],1./rt_mm); rm = nancorrcoef(y_pred_m,1./nanmean(RT_mal,2));


% Calculating the dissimilarity between the odd and even numbered subjects
ctt = nancorrcoef(1./nanmean(RT_tel(:,1:2:end),2),1./nanmean(RT_tel(:,2:2:end),2));
cmm = nancorrcoef(1./nanmean(RT_mal(:,1:2:end),2),1./nanmean(RT_mal(:,2:2:end),2));
