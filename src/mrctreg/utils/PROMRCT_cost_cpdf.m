function [cost vol_cost_mse vol_cost_cpdf1 vol_cost_cpdf2] = PROMRCT_cost_cpdf(vol_ct,volmu_ct,volsd_ct,vol_t2,ivec_mask,nlogcpdf1,nlogcpdf2,bins_ct,bins_t2,wt_mse,wtvol_mse,wtvol_cpdf)

dims = size(vol_ct);

if ~exist('wt_mse','var')
  wt_mse = 0.5;
end

if ~exist('wtvol_mse','var') | isempty(wtvol_mse)
  wtvol_mse = ones(dims);
end

if ~exist('wtvol_cpdf','var') | isempty(wtvol_cpdf)
  wtvol_cpdf = ones(dims);
end

v_mse = min(10,((vol_ct(ivec_mask)-volmu_ct(ivec_mask))./volsd_ct(ivec_mask)).^2);
cost_mse = sum(wtvol_mse(ivec_mask).*v_mse)/sum(wtvol_mse(ivec_mask));

v_t2 = max(bins_t2(1),min(bins_t2(end),vol_t2(ivec_mask)));
v_ct = max(bins_ct(1),min(bins_ct(end),vol_ct(ivec_mask)));
v_cost1 = interp2(bins_ct,bins_t2,nlogcpdf1,v_ct,v_t2,'linear');
v_cost2 = interp2(bins_ct,bins_t2,nlogcpdf2,v_ct,v_t2,'linear');

cost_cpdf = sum(wtvol_cpdf(ivec_mask).*(v_cost1+v_cost2))/sum(wtvol_cpdf(ivec_mask));

cost = wt_mse*cost_mse + (1-wt_mse)*cost_cpdf;

if nargout>1
  vol_cost_mse = zeros(dims); vol_cost_mse(ivec_mask) = wtvol_mse(ivec_mask).*v_mse;
  vol_cost_cpdf1 = zeros(dims); vol_cost_cpdf1(ivec_mask) = wtvol_cpdf(ivec_mask).*v_cost1;
  vol_cost_cpdf2 = zeros(dims); vol_cost_cpdf2(ivec_mask) = wtvol_cpdf(ivec_mask).*v_cost2;
%  showVol(vol_ct,vol_t2,vol_cost_mse,vol_cost_cpdf1,vol_cost_cpdf2)
end


