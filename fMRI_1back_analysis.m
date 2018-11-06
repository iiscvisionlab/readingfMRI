% Analysing RTs for one-back task
% Accuracy and RTs are not correlated
allclear
ismal = [1 1 0 0 1 1 0 0 1 1 0 1 0 1 0 0 0 1 1 1 0 1 0 0 0 0 0 0 1 0 1 1 1 1 1];
f = dir('../preprocessing/*_READING_*');

Acc = zeros(35,64);
RT = NaN(35,64);
for sub = 1:24
    f1 = dir(['../preprocessing/' f(sub).name '/behaviour/*evt*']);
    for run = 1:8
        load(['../preprocessing/' f(sub).name '/behaviour/' f1(run).name]);
        repid = find(expt_str.data.isrep); stim = expt_str.data.stimid(repid);
        for R = 1:numel(repid)
            q = expt_str.data.keytime{repid(R)};
            if numel(q) > 50
                Acc(sub,stim(R)) = 1; RT(sub,stim(R)) = expt_str.data.keytime{repid(R)}(1) - expt_str.data.tstimon(repid(R));
            end
        end
    end
end
%%
for sub = 25:26
    f1 = dir(['../preprocessing/' f(sub).name '/behaviour/*evt*']);
    for run = 1:8
        load(['../preprocessing/' f(sub).name '/behaviour/' f1(run).name]);
        repid = find(expt_str.data.isrep); stim = expt_str.data.stimid(repid);
        for R = 1:numel(repid)
            q = expt_str.data.keytime{repid(R)};
            if numel(q) > 50
                qx = find(expt_str.data.responsekey{repid(R)} == 55);
                Acc(sub,stim(R)) = 1; RT(sub,stim(R)) = expt_str.data.keytime{repid(R)}(qx(1)) - expt_str.data.tstimon(repid(R));
            end
        end
    end
end
%%
  
for sub = 27:numel(f)
    f1 = dir(['../preprocessing/' f(sub).name '/behaviour/*evt*']);
    for run = 1:8
        load(['../preprocessing/' f(sub).name '/behaviour/' f1(run).name]);
        repid = find(expt_str.data.isrep); stim = expt_str.data.stimid(repid);
        for R = 1:numel(repid)
            q = expt_str.data.keytime{repid(R)};
            if numel(q) > 2
                Acc(sub,stim(R)) = 1; RT(sub,stim(R)) = expt_str.data.keytime{repid(R)}(1) - expt_str.data.tstimon(repid(R));
            end
        end
    end
end

%% Analysing RT and accuracy data
qmal = 1:34; qtel = 35:64;
RTmm = nanmean(RT(ismal == 1,qmal),2); acc_mm = nanmean(Acc(ismal == 1,qmal),2);
RTmt = nanmean(RT(ismal == 1,qtel),2); acc_mt = nanmean(Acc(ismal == 1,qtel),2);
RTtt = nanmean(RT(ismal == 0,qtel),2); acc_tt = nanmean(Acc(ismal == 0,qtel),2);
RTtm = nanmean(RT(ismal == 0,qmal),2); acc_tm = nanmean(Acc(ismal == 0,qmal),2);
figure; subplot(121); barweb([nanmean([RTmm RTmt]) nanmean([RTtt RTtm])], [nansem([RTmm RTmt]) nansem([RTtt RTtm])]);
ylabel('Mean Reaction Time (s)'); legend('MM','MT','TT','TM')
subplot(122); barweb([nanmean([acc_mm acc_mt]) nanmean([acc_tt acc_tm])],[nansem([acc_mm acc_mt]) nansem([acc_tt acc_tm])]); 
ylabel('Accuracy'); legend('MM','MT','TT','TM')

