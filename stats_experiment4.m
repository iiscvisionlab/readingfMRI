allclear
load('L2fmri_bigram.mat')

L2_str.RT(L2_str.RT < .3) = NaN;  % removing accidental key press.

images = L2_str.images;
image_pairs = L2_str.img_pairs;

RT_mal = L2_str.RT(:,L2_str.subjinfo.ismalayalam == 1,:);
RT_tel = L2_str.RT(:,L2_str.subjinfo.ismalayalam == 0,:);

mRT_mal = nanmean(RT_mal,3);
mRT_tel = nanmean(RT_tel,3);

% Malayalam and telugu pairs
mal = 1:276; tel = 277:552;  

% Average across subjects
RT_mal_mean = nanmean(mRT_mal,2);
RT_tel_mean = nanmean(mRT_tel,2);
RTmm = RT_mal_mean(mal);
RTmt = RT_mal_mean(tel);
RTtt = RT_tel_mean(tel);
RTtm = RT_tel_mean(mal);

% splithalfcorr(mean(RT_tel,3)')
% splithalfcorr(mean(RT_mal,3)')

% Calculating the dissimilarity between the odd and even numbered subjects
ctt = nancorrcoef(1./nanmean(mRT_tel(tel,1:2:end),2),1./nanmean(mRT_tel(tel,2:2:end),2));
cmt = nancorrcoef(1./nanmean(mRT_mal(tel,1:2:end),2),1./nanmean(mRT_mal(tel,2:2:end),2));
cmm = nancorrcoef(1./nanmean(mRT_mal(mal,1:2:end),2),1./nanmean(mRT_mal(mal,2:2:end),2));
ctm = nancorrcoef(1./nanmean(mRT_tel(mal,1:2:end),2),1./nanmean(mRT_tel(mal,2:2:end),2));

%% Correlating with predicted data
load pred_2RT_fmri
figure; 
subplot(221); corrplot(1./RTtt, pred_tt,'Tel sub on Tel pairs',1);  xlabel('Experimental dissimialirites'); ylabel('Model dissimilarities')
subplot(224); corrplot(1./RTtm, pred_tm,'Tel sub on Mal pairs',1);  xlabel('Experimental dissimialirites'); ylabel('Model dissimilarities')
subplot(223); corrplot(1./RTmm, pred_mm,'Mal sub on Mal pairs',1);  xlabel('Experimental dissimialirites'); ylabel('Model dissimilarities')
subplot(222); corrplot(1./RTmt, pred_mt,'Mal sub on Tel pairs',1);  xlabel('Experimental dissimialirites'); ylabel('Model dissimilarities')


