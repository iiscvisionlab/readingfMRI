% Compare mean/median of two variables as appropriate 
% Required inputs
%    x             = first set of samples
%    y             = second set of samples
% Optional inputs:
%    pairedflag    = if 1, run a paired test 
%                        (defaults: 1 if x & y have same lengths, 0 otherwise)
%    pthresh       = p-value threshold for normality test (default = 0.05)
%    figflag       = if 1, generate bar graph with error bars
% Outputs:
%    pout          = p-value of test of means or medians, as appropriate
%    pxnormal      = p-value of lillietest for normality on x
%    pynormal      = p-value of lillietest for normality on y
%    xm            = mean or median of x (as applicable)
%    ym            = mean or median of y (as applicable)
%    xe            = std or inter-quartile range of x (as applicable)
%    ye            = std or inter-quartile range of y (as applicable)
% Method
%    statcomparemean performs a lillietest for normality on the two input variables. 
%    If the two variables are normal, it performs a t-test
%    If the two variables are not normal, it performs a non-parametric Wilcoxon's ranksum test

% Changelog
%    10 Aug 2013    (SPA) First version with only t-test
%    22 Apr 2015    (SPA) Appended ranksum/signrank non-parametric version
%    20 May 2015    (SPA) Added normality test
%    18 Jul 2015    (SPA) Returns means/medians and std/interquartile range as applicable
%    18 Jul 2015    (SPA) changed normality test from lillietest to adtest
%                         based on https://en.wikipedia.org/wiki/Normality_test
%    26 Nov 2015    (SPA) changed order of outputs to pout,pxnormal,pynormal first

function [pout,pxnormal,pynormal,xm,ym,xe,ye] = statcomparemean(x,y,pairedflag,pthresh,figureflag)
if(~exist('pairedflag')), if(length(x)==length(y)), pairedflag = 1; else pairedflag = 0; end; end; 
if(~exist('pthresh')), pthresh = 0.05; end; 
if(~exist('figureflag')), figureflag = 1; end; 
warning('off','stats:adtest:OutOfRangePLow'); 
dispflag = 1; if(nargout>0), figureflag = 0; dispflag = 0; end; 
x = x(:); y = y(:); 

xname = inputname(1); yname = inputname(2); 
if(isempty(xname)),xname='x';end; if(isempty(yname)),yname='y';end; 

if(pairedflag==1)
    [~,pttest] = ttest(x,y); pairstr = 'paired test'; 
    pranksum = signrank(x,y); 
else
    [~,pttest] = ttest2(x,y); pairstr = 'unpaired test'; 
    pranksum = ranksum(x,y); 
end

[~,pxnormal] = adtest(x); [~,pynormal] = adtest(y); 
isnormal = pxnormal>pthresh & pynormal>pthresh; % 1 => both are normal :)
if(isnormal)
    if(dispflag)
        fprintf('adtest for normality PASSED: %s (p=%2.2g), %s (p=%2.2g) \n',xname,pxnormal,yname,pynormal);
        fprintf('Means (t-test) : %s = %2.3g (n=%d), %s = %2.3g (n=%d), p = %2.3g, %s \n',xname,nanmean(x),length(x),yname,nanmean(y),length(y),pttest,pairstr);
    end
    xm = nanmean(x); ym = nanmean(y); xe = nanstd(x); ye = nanstd(y); 
    pout = pttest; 
else
    if(dispflag)
        fprintf('adtest for normality FAILED: %s (p=%2.2g), %s (p=%2.2g) \n', xname,pxnormal,yname,pynormal);
        fprintf('Medians (ranksum): %s = %2.2g (n=%d), %s = %2.2g (n=%d), p = %2.2g, %s \n',xname,nanmedian(x),length(x),yname,nanmedian(y),length(y),pranksum,pairstr);
    end
    xm = nanmedian(x); ym = nanmedian(y); xe = iqr(x); ye = iqr(y); pout = pranksum; 
end

if(figureflag)
    m = [xm ym]; e = [xe ye]/sqrt(length(x)); 
    bar(m); hold on; errorbar([1:2],m,e,'b.'); 
    text(1.4,max(m)*1.1,sprintf('p = %2.2g',pout)); 
    set(gca,'XTickLabel',{xname yname}); 
end

return