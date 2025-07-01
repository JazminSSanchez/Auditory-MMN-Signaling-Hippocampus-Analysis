
function CSI = CSI_bootstrap2( Dev, Std, N, reproducible )
%UNTITLED Summary of this function goes here
%   Dev: <nd x 2> single trial responses to f1 and f2
%   Std: <ns x 2> single trial responses to f1 and f2
%   N: (int) number of bootstrap samples
%   reproducible: (bool) if true, the sampling will be always the same
%   CSI : struct with fieds: mean, ci, dist

% make sure Dev and Std the same size, repeating values when necessary
% nDev = size(Dev,1);
% nStd = size(Std,1);
% nCommon = lcm(nDev, nStd);
% Dev = repmat(Dev,[nCommon/nDev 1]);
% Std = repmat(Std,[nCommon/nStd 1]);

% go
if nargin < 4
    reproducible = 1;
end

if isempty(Dev) || isempty (Std)
    CSI.distr = nan;
    CSI.mean_bstrp = nan;
    CSI.median_bstrp = nan;
    CSI.ci = [nan nan];
else
    csi = @(dev,std) ((sum(dev,2)-sum(std,2))./(sum(dev,2)+sum(std,2)));
    
    % CSI.mean = csi(mean(Dev, 'omitnan'), mean(Std, 'omitnan'));
    % CSI.median = csi(median(Dev, 'omitnan'), median(Std, 'omitnan'));
    if reproducible
        rng('default')
    end
    boot_dev = bootstrp(N,@mean,Dev);
    if reproducible
        rng('default')
    end
    boot_std = bootstrp(N,@mean,Std);
    
    CSI.csi = mean(csi(Dev, Std), 'omitnan');
    CSI.distr = csi(boot_dev, boot_std);
    CSI.mean_bstrp = mean(CSI.distr, 'omitnan');
    CSI.median_bstrp = median(CSI.distr, 'omitnan');
    CSI.ci = prctile(CSI.distr, [2.5 97.5]);
end




% 
% 
%     function csi = lcf_compute_csi(d,s)
%        
%         csi = normdiff(mean(mean(d,1),2),mean(mean(s,1),2));
%     end



end
