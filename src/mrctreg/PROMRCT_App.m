function result = PROMRCT_App(pathToCTDicoms, pathToT2Dicoms, outputDir, varargin)
%% Default Inputs -- Parse the default inputs 
  defaults = {
          'dispflag',0;
          'outputMatFlag',0;
          'outputDicomFlag',0;
          'ManualRegistrationMatrix',[];
	  'ProstateContourPath','';
          };

  parse_varargin_auto(varargin,defaults);
  
  starttime = now;
  fprintf(1,'***Start*** %s\n',datestr(starttime));
  
  if ~exist(outputDir, 'dir')
     mkdir(outputDir)
  end
  
%% Load atlas
  load('volmu_ct_scaled.mat');
  volmask_prostate = volmu_t2_seg>0.5; 
  volmask_prostate_dilate10 = distvol_prostate<10; 
  volmask_prostate_dilate20 = distvol_prostate<20; 
  volmask_prostate_dilate30 = distvol_prostate<30; 
  volmask_prostate_dilate40 = distvol_prostate<40; 
  volmask_prostate_dilate50 = distvol_prostate<50;

  volmu_ct_ctx = mgh2ctx(volmu_ct, M_ref);
  volsd_ct_ctx = mgh2ctx(volsd_ct, M_ref);
  volmu_ct_seg_ctx = mgh2ctx(volmu_ct_seg, M_ref);
  volmu_t2_ctx = mgh2ctx(volmu_t2, M_ref);
  volmu_t2_seg_ctx = mgh2ctx(volmu_t2_seg, M_ref);
  
  bbox_ref = getfield(regionprops(volmu_t2_seg_ctx.imgs>0.5,'BoundingBox'),'BoundingBox');
  
  atlasReferenceVolumeCTX = vol_ref_ctx;
  t2ReferenceVolumeCTX = vol_ref_ctx;
  
%% Load CT
  bins_ct = -200:1:200;
  fprintf('Loading CT: %s\n', pathToCTDicoms);
  vol_ct_ctx = QD_read_dicomdir(pathToCTDicoms);
  tmpCTVolumeCTX = vol_ct_ctx;
  vol_ct_ctx.imgs = max( min(bins_ct), min(max(bins_ct), tmpCTVolumeCTX.imgs) );
  
%% Load T2
  fprintf('Loading T2: %s\n', pathToT2Dicoms);
  vol_t2_ctx = QD_read_dicomdir(pathToT2Dicoms);
  
%% Generate prostate mask from T2
  fname_T2 = fullfile(outputDir, 'T2_corrected.mgz');
  QD_ctx_save_mgh( vol_t2_ctx, fname_T2 );
  
  if isempty(ProstateContourPath)
    fprintf('%s -- %s.m:    Segmenting prostate from T2 volume using CMIG software...\n',datestr(now),mfilename);
    if isdeployed
      container = 'docker';
    else
      container = 'singularity';
    end
    contour = contour_prostate_cmig(fname_T2, container);
  else
    contour = QD_ctx_load_mgh(ProstateContourPath);
  end
  vol_t2_seg_ctx = contour;
  
%% Do some intensity corrections for the T2 Volume
  fprintf('N3 corrections...\n');
  vol_t2_ctx.imgs = vol_t2_ctx.imgs*200/median(vol_t2_ctx.imgs(vol_t2_seg_ctx.imgs>0.5));
  for iter_n3 = 1:4 
    vol_t2_ctx = vol_correct_bias_field_n3_amd(vol_t2_ctx,[],[],[],[],[],[],[],50);
  end 
  vol_t2_ctx.imgs = vol_t2_ctx.imgs*200/median(vol_t2_ctx.imgs(vol_t2_seg_ctx.imgs>0.5));
  
%% Calculate Transformation for T2 atlas registration
  patientT2Centroid = vol_t2_seg_ctx.Mvxl2lph(:,[2 1 3 4]) * [regionprops3(vol_t2_seg_ctx.imgs>0.5,'Centroid').Centroid 1]';
  t2ReferenceVolumeMidpoint = t2ReferenceVolumeCTX.Mvxl2lph(:,[2 1 3 4]) * [(size(t2ReferenceVolumeCTX.imgs)+1)/2 1]';
  M_trans_t2 = Mtrans(patientT2Centroid(1)-t2ReferenceVolumeMidpoint(1), patientT2Centroid(2)-t2ReferenceVolumeMidpoint(2), patientT2Centroid(3)-t2ReferenceVolumeMidpoint(3));
  
%% Do some rescaling+shifting based on prostate mask
  vol_seg = vol_resample(vol_t2_seg_ctx, t2ReferenceVolumeCTX, M_trans_t2); 
  vol_seg = vol_seg.imgs;
  bbox_subj = getfield(regionprops(vol_seg>0.5,'BoundingBox'),'BoundingBox');
  svec = bbox_subj([5 4 6])./bbox_ref([5 4 6]);
  M_scale = Mscale(svec(1),svec(2),svec(3));
  
  vol_t2_ctr_ctx = vol_resample(vol_t2_ctx, t2ReferenceVolumeCTX, M_trans_t2);
  vol_t2_sca_ctx = vol_resample(vol_t2_ctr_ctx, volmu_t2_seg_ctx, M_scale);
  
  vol_t2_seg_ctr_ctx = vol_resample(vol_t2_seg_ctx, t2ReferenceVolumeCTX, M_trans_t2);
  vol_t2_seg_sca_ctx = vol_resample(vol_t2_seg_ctr_ctx, atlasReferenceVolumeCTX, M_scale);
  
  
%% Load in joint pdf
  fprintf('Loading joint PDF...\n');
  dilstr = 'dilate30';
  fname_cpdf = sprintf('cpdf_%s.mat',dilstr);
  load(fname_cpdf);
  
  volmask = eval(sprintf('volmask_prostate_%s',dilstr))>0.5;
  ivec_mask = find(volmask);
  
  M_scale_t2 = M_scale;
  
  
%% Scale CT volume
  c0 = vol_ct_ctx.Mvxl2lph * [1 1 1 1]';
  c1 = vol_ct_ctx.Mvxl2lph * [vol_ct_ctx.dimc vol_ct_ctx.dimr vol_ct_ctx.dimd 1]'; 
  zmin = min(c0(3),c1(3));
  M_tmp = Mtrans(0,0,zmin+200); 
  M_scale_ct = M_scale*M_tmp; % Assume that prostate is approximately some fixed distance from bottom, shift and scale accordingly
  vol_ct_sca2_ctx = vol_resample(vol_ct_ctx, atlasReferenceVolumeCTX, M_scale_ct);
  
  
%% Create surface meshes for prostate mask in different spaces
  smf = 3;
  
  vol = vol_t2_seg_ctx;
  volsm = vol; volsm.imgs = smooth3(volsm.imgs,'gaussian',1+2*ceil(smf),smf);
  S = isosurface(volsm.imgs,0.5);
  V = cat(2,S.vertices,ones(size(S.vertices,1),1))';
  S.vertices = (vol.Mvxl2lph(1:3,[2 1 3 4])*V)';
  S_native_t2 = S;

  vol = vol_t2_seg_ctr_ctx;
  volsm = vol; volsm.imgs = smooth3(volsm.imgs,'gaussian',1+2*ceil(smf),smf);
  S = isosurface(volsm.imgs,0.5);
  V = cat(2,S.vertices,ones(size(S.vertices,1),1))';
  S.vertices = (vol.Mvxl2lph(1:3,[2 1 3 4])*V)';
  S_ctr_t2 = S;

  vol = vol_t2_seg_sca_ctx;
  volsm = vol; volsm.imgs = smooth3(volsm.imgs,'gaussian',1+2*ceil(smf),smf);
  S = isosurface(volsm.imgs,0.5);
  V = cat(2,S.vertices,ones(size(S.vertices,1),1))';
  S.vertices = (vol.Mvxl2lph(1:3,[2 1 3 4])*V)';
  S_scale_t2 = S;

%% Perform initial rough rigid-body registration -- Should replace with prostate segmentation method (Use DL?)
  fprintf('Rigid body registration...\n');
  % Find body in y-dimension
  smf = 5;
  vol_tmp = vol_ct_sca2_ctx.imgs;
  v3 = squeeze(mean(mean(abs(vol_tmp),2),1));
  c3_min = min(find(v3>0)); c3_max = max(find(v3>0));
  v2 = zeros(1,size(vol_tmp,2));
  for i2 = 1:length(v2)
    im = squeeze(vol_tmp(:,i2,:));
    v2(i2) = mean(var(im(:,round(c3_min+0.5*(c3_max-c3_min)):c3_max),[],2));
  end
  v2_sm = real(smooth1d_amd(v2,50));
  [mv c2_mid] = max(v2_sm); c2_min = min(find(v2_sm>0.3*v2_sm(c2_mid))); c2_max = max(find(v2_sm>0.3*v2_sm(c2_mid))); 
  c2 = c2_mid;
  
  volmu = volmu_ct; volsd = volsd_ct; volmu_seg = volmu_t2_seg;
  volmu_ctx = mgh2ctx(volmu, M_ref);
  
  % Find location of butt crack to initialize 3D grid search -- should allow tightening of search range
  %   Make coronal mean image from ctVolumeScaled_ctx
  %   Look for local minimum 
  im_mean_cor = squeeze(mean(vol_ct_sca2_ctx.imgs(:,c2_min:c2_max,:),2));
  im_mip_cor = squeeze(max(vol_ct_sca2_ctx.imgs(:,c2_min:c2_max,:),[],2)); % Get rid of table in back
  im_mip_sag = squeeze(max(vol_ct_sca2_ctx.imgs,[],1));
  v1 = mean(im_mean_cor(:,round(0.75*max(find(max(abs(im_mean_cor))>0))):end),2); v1_sm = real(smooth1d_amd(v1,20)); % v1_sm = smooth(v1,200);
  im = double(im_mean_cor(:,round(0.75*max(find(max(abs(im_mean_cor))>0))):end));
  tmp = getfield(regionprops(true(size(im)),im-min(im(:)),'WeightedCentroid'),'WeightedCentroid'); c1 = round(tmp(2)); % Should subtract min(im_mean_cor(:)) from im_mean_cor (not clear what WeightedCentroid does  with negative intensity values)
  v3 = mean(im_mip_cor(c1+[-20:20],:),1); v3_sm = smooth(v3,10);
  c3 = max(find(v3_sm>150))-20;
  
  im_mip_axi = squeeze(max(vol_ct_sca2_ctx.imgs(:,:,c3_min:c3_max),[],3));
  
  if dispflag
    figure(1); clf;
    subplot(2,2,1); imagesc(im_mean_cor); colormap(gray); line(xlim,[c1 c1],'color',[0 1 0]); line([c3 c3],ylim,'color',[0 1 0]);
    subplot(2,2,2); imagesc(im_mip_cor); colormap(gray); line(xlim,[c1 c1],'color',[0 1 0]); line([c3 c3],ylim,'color',[0 1 0]);
    subplot(2,2,3); imagesc(im_mip_sag); colormap(gray); line(xlim,[c2_min c2_min],'color',[1 0 0]); line(xlim,[c2 c2],'color',[0 1 0]); line(xlim,[c2_max c2_max],'color',[0 0 1]); line([c3 c3],ylim,'color',[0 1 0]);
    subplot(2,2,4); imagesc(im_mip_axi); colormap(gray); line([c2_min c2_min],ylim,'color',[1 0 0]); line([c2 c2],ylim,'color',[0 1 0]); line([c2_max c2_max],ylim,'color',[0 0 1]);
    figure(2); clf;
    subplot(2,2,1); plot(v1_sm); line([c1 c1],ylim,'color',[0 1 0])
    subplot(2,2,2); plot(v2_sm); line([c2_min c2_min],ylim,'color',[1 0 0]); line([c2 c2],ylim,'color',[0 1 0]); line([c2_max c2_max],ylim,'color',[0 0 1])
    subplot(2,2,3); plot(v3_sm); line([c3 c3],ylim,'color',[0 1 0])
    drawnow
  end
  
  xhat0 = rowvec(vol_ct_sca2_ctx.Mvxl2lph(1:3,:)*[c1 c2 c3 1]');
  
  
%% Perform 3D grid search -- look for faster algorithm
  
  nlogcpdf1 = -log(max(1e-3,jpdf_sum./max(10000,sum(jpdf_sum,2)))); % Conditional pdf of ct given t2
  nlogcpdf2 = -log(max(1e-3,jpdf_sum./max(10000,sum(jpdf_sum,1)))); % Conditional pdf of t2 given ct
  
  distvol_prostate_subj = bwdist(vol_t2_seg_sca_ctx.imgs); % Use distance from OnQ Prostate mask for spatial weighting of cpdf terms
  
  options = optimset('Display','none','MaxFunEvals',100);
  
  d1vec = -30:10:30;
  d2vec = -30:10:30;
  d3vec = -60:10:60; % May be able to tighten this range
  xhat = xhat0;
  for pass = 1:2
    xhat_bak = xhat;
    if pass==1
      wt_mse = 1.0; % use mse only 
      wtvol_mse = distvol_prostate_subj<=30; 
      wtvol_cpdf = distvol_prostate_subj<=30; 
    else
      wt_mse = 0.5; % Mix of mse and cpdf
      wtvol_mse = distvol_prostate_subj<=5; 
      wtvol_cpdf = distvol_prostate_subj<=5; 
    end
    volmask = max(wtvol_mse,wtvol_cpdf);
    volmask_box = false(size(volmask)); indvec = find(volmask>0.5); [ivec2 ivec1 ivec3] = ind2sub(size(volmask),indvec); volmask_box(min(ivec1):max(ivec1),min(ivec2):max(ivec2),min(ivec3):max(ivec3)) = true;
    dims = size(atlasReferenceVolumeCTX.imgs);
    dims_cropped = [1+max(ivec1)-min(ivec1) 1+max(ivec2)-min(ivec2) 1+max(ivec3)-min(ivec3)];
    vol_cropped = zeros(dims_cropped); [vol_ref M_ref] = ctx2mgh(atlasReferenceVolumeCTX); M_cropped = M_ref; M_cropped(1:3,4) = M_ref(1:3,:)*[min(ivec1)-1 min(ivec2)-1 min(ivec3)-1 1]'; vol_cropped_ctx = mgh2ctx(vol_cropped,M_cropped);
    volmask_cropped = volmask(min(ivec1):max(ivec1),min(ivec2):max(ivec2),min(ivec3):max(ivec3)); ivec_mask_cropped = find(volmask_cropped>0.5);
    volmu_cropped = volmu(min(ivec1):max(ivec1),min(ivec2):max(ivec2),min(ivec3):max(ivec3));
    volsd_cropped = volsd(min(ivec1):max(ivec1),min(ivec2):max(ivec2),min(ivec3):max(ivec3));
    volmu_seg_cropped = volmu_seg(min(ivec1):max(ivec1),min(ivec2):max(ivec2),min(ivec3):max(ivec3));
    volmu_cropped_ctx = mgh2ctx(volmu_cropped,M_cropped);
    vol_t2_sca_cropped = vol_t2_sca_ctx.imgs(min(ivec1):max(ivec1),min(ivec2):max(ivec2),min(ivec3):max(ivec3));
    wtvol_mse_cropped = wtvol_mse(min(ivec1):max(ivec1),min(ivec2):max(ivec2),min(ivec3):max(ivec3)); ivec_mask_cropped = find(volmask_cropped>0.5);
    wtvol_cpdf_cropped = wtvol_cpdf(min(ivec1):max(ivec1),min(ivec2):max(ivec2),min(ivec3):max(ivec3)); ivec_mask_cropped = find(volmask_cropped>0.5);
    distvol_prostate_cropped = distvol_prostate(min(ivec1):max(ivec1),min(ivec2):max(ivec2),min(ivec3):max(ivec3));
    resfun = @(vol,x) vol_resample(vol,atlasReferenceVolumeCTX,Mtrans(x(1),x(2),x(3)));
    resfun_cropped = @(vol,x) vol_resample(vol,vol_cropped_ctx,Mtrans(x(1),x(2),x(3)));
    costfun = @(x) PROMRCT_cost_cpdf(getfield(resfun_cropped(vol_ct_sca2_ctx,x),'imgs'),volmu_cropped,volsd_cropped,vol_t2_sca_cropped,ivec_mask_cropped,nlogcpdf1,nlogcpdf2,bins_ct,bins_t2,wt_mse,wtvol_mse_cropped,wtvol_cpdf_cropped);
  
    if pass==1
      costmat = NaN(length(d1vec),length(d2vec),length(d3vec));
      tic
      for di3 = 1:length(d3vec)
        for di1 = 1:length(d1vec)
          for di2 = 1:length(d2vec)
            x = xhat_bak;
            x(1) = xhat_bak(1) + d1vec(di1);
            x(2) = xhat_bak(2) + d2vec(di2);
            x(3) = xhat_bak(3) + d3vec(di3);
            costmat(di1,di2,di3) = costfun(x);
  %          fprintf(1,'ci=%d d=%f: cost=%e\n',ci,dvec(di),costvec(di));
          end
        end
      end
      toc
      [cost1 mi] = min(costmat(:));
      [di1 di2 di3] = ind2sub(size(costmat),mi);
      xhat_bak = xhat_bak + [d1vec(di1) d2vec(di2) d3vec(di3)];
    end
  
    tic
    [xhat cost] = fminsearch(costfun,xhat_bak,options);
    toc
    fprintf(1,'pass=%d xhat=[%f %f %f] cost=%e\n',pass,xhat,cost);
  end
  
  M_trans_ct = Mtrans(xhat(1),xhat(2),xhat(3));
  M_ct_atlas = M_scale_ct*M_trans_ct; % Mapping from subject native CT to scaled atlas space
  
  M_reg = M_trans_t2*M_scale_t2*inv(M_ct_atlas);
  S = S_native_t2;
  V = cat(2,S.vertices,ones(size(S.vertices,1),1))';
  M = inv(M_reg);
  S.vertices = (M(1:3,:)*V)';
  S_native_ct = S; S_native_ct.color = [0 1 1]; S_native_ct.linewidth = 2;
  
  vol_ct_ctr_ctx = vol_resample(vol_ct_ctx,volmu_ct_ctx,M_ct_atlas); % map ct to atlas coordinates
  vol_t2_res_ctx = vol_resample(vol_t2_ctx,vol_ct_ctx,M_reg); % Map from T2 native coordinates to CT native coordinates
  vol_t2_seg_res_ctx = vol_resample(vol_t2_seg_ctx,vol_ct_ctx,M_reg); % Map from OnQ seg in native T2 space to native CT space
  volmu_ct_res_ctx = vol_resample(volmu_ct_ctx,vol_ct_ctx,inv(M_ct_atlas)); % map mu ct from atlas back to CT native coordinates
  volsd_ct_res_ctx = vol_resample(volsd_ct_ctx,vol_ct_ctx,inv(M_ct_atlas)); % map mu ct stdv from atlas back to CT native coordinates
  volmu_ct_seg_res_ctx = vol_resample(volmu_ct_seg_ctx,vol_ct_ctx,inv(M_ct_atlas)); % map mean ct seg from atlas back to CT native coordinates
  volmu_t2_res_ctx = vol_resample(volmu_t2_ctx,vol_ct_ctx,inv(M_ct_atlas)); % map mu t2 from atlas back to CT native coordinates
  volmu_t2_seg_res_ctx = vol_resample(volmu_t2_seg_ctx,vol_ct_ctx,inv(M_ct_atlas)); % map mean t2 seg from atlas back to CT native coordinates
  vol_ct_sca_res_ctx = resfun(vol_ct_sca2_ctx,xhat);
  [cost vol_cost_mse vol_cost_cpdf1 vol_cost_cpdf2] = PROMRCT_cost_cpdf(vol_ct_sca_res_ctx.imgs,volmu,volsd,vol_t2_sca_ctx.imgs,ivec_mask,nlogcpdf1,nlogcpdf2,bins_ct,bins_t2,wt_mse,wtvol_mse,wtvol_cpdf);


  %% Upsampling stuff to address blurry T2 images post registration
  %% Crop CT volume to resampled T2 data, and then created upsampled version
  BB = getfield(regionprops(vol_t2_res_ctx.imgs~=0),'BoundingBox');
  i1_start = floor(BB(2)); i2_start = floor(BB(1)); i3_start = floor(BB(3));
  i1_end = i1_start+BB(5); i2_end = i2_start+BB(4); i3_end = i3_start+BB(6);
  [vol M] = ctx2mgh(vol_ct_ctx);
  vol_cropped = vol(i1_start:i1_end,i2_start:i2_end,i3_start:i3_end,:);
  M_cropped = M; M_cropped(1:3,4) = M_cropped(1:3,4) + M(1:3,:)*[i1_start i2_start i3_start 1]' - M_cropped(1:3,:)*[1 1 1 1]';

  vx_dims_cropped = rowvec(ceil(abs(diag(M_cropped(1:3,1:3)))/1));
  vx_dims_orig_t2 = rowvec(abs(diag(vol_t2_ctx.Mvxl2lph(1:3,1:3))/1));
  ssvec = round(vx_dims_cropped./vx_dims_orig_t2);
  ssvec(3) = ceil(abs(vol_t2_ctx.Mvxl2lph(3,3)));
  ssvec(ssvec<1) = 1; % Don't downsample if resolution of cropped CT is higher than original T2

  M_cropped_up = M_cropped;
  M_cropped_up(1:3,1:3) = M_cropped(1:3,1:3)./rowvec(ssvec);
  vol_cropped_up = zeros(size(vol_cropped).*ssvec);
  vol_cropped_up_ctx = mgh2ctx(vol_cropped_up,M_cropped_up);

  fprintf('----- Resampling registered T2 to a higher resolution -----\n')
  fprintf('Rescaling vector: [%s %s %s]\n', num2str(ssvec(1)), num2str(ssvec(2)), num2str(ssvec(3)));


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Resample T2 to upsampled, cropped CT space
  vol_t2_res_up_ctx = vol_resample(vol_t2_ctx, vol_cropped_up_ctx, M_reg, 2); % This is the upsampled T2 volume that should be saved
  %% Resample autoseg to upsampled, cropped CT space
  vol_t2_seg_res_up_ctx = vol_resample(vol_t2_seg_ctx, vol_cropped_up_ctx, M_reg, 2);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
  %% QA stuff
  result = struct();
  coordmat = RSI_find_centers(vol_t2_seg_res_ctx.imgs); 
  result.coordmat = coordmat; result.M_reg = M_reg; result.vol_ct_ctx = vol_ct_ctx; result.vol_t2_res_ctx = vol_t2_res_ctx; result.S_native_ct = S_native_ct;
  if ~isempty(ManualRegistrationMatrix)
    vol_t2_seg_res2_ctx = vol_resample(vol_t2_seg_ctx,vol_ct_ctx,ManualRegistrationMatrix);
    vol_t2_res2_ctx = vol_resample(vol_t2_ctx,vol_ct_ctx,ManualRegistrationMatrix);
    coordmat2 = RSI_find_centers(vol_t2_seg_res2_ctx.imgs);
    disp(coordmat2(1,:)-coordmat(1,:)) % This is the displacement of the center of the prostate mas
    S = S_native_t2;
    V = cat(2,S.vertices,ones(size(S.vertices,1),1))';
    M = inv(ManualRegistrationMatrix);
    S.vertices = (M(1:3,:)*V)';
    S_native_ct2 = S; S_native_ct2.color = [0 1 0]; S_native_ct2.linewidth = 2;
    result.coordmat2 = coordmat2; result.ManualRegistrationMatrix = ManualRegistrationMatrix; result.S_native_ct2 = S_native_ct2;
  end
 
  if dispflag
    if size(coordmat,1)>0
      if ~isempty(ManualRegistrationMatrix)
        showVol(vol_ct_ctx,vol_t2_res_ctx,vol_t2_res2_ctx,S_native_ct,S_native_ct2,struct('RCS',round(coordmat(1,:)))); % Resample data using ManualRegistrationMatrix, display result, with contour
      else
        showVol(vol_ct_ctx,vol_t2_res_ctx,S_native_ct,struct('RCS',round(coordmat(1,:)))); % Resample data using ManualRegistrationMatrix, display result, with contour
      end
    end
  end
   
  
  if isequal(outputMatFlag,1)
     fprintf('Saving Output .mat Files \n')
     outputMatDir = fullfile(outputDir,'/mat/');
     if ~exist(outputMatDir, 'dir')
      mkdir(outputMatDir)
     end
    save(fullfile(outputMatDir,'output.mat'),'M_scale_t2','cost','M_ct_atlas','M_trans_t2','vol_t2_res_ctx','vol_t2_seg_res_ctx', ...
        'vol_ct_ctx','vol_ct_sca_res_ctx','volmu_t2_res_ctx', ...
                    'volmu_ct_res_ctx','volsd_ct_ctx','volmu_ct_seg_res_ctx', ...
        'volmu_t2_seg_res_ctx','vol_t2_seg_ctr_ctx', 'vol_t2_seg_sca_ctx', ...
        'S_native_ct','vol_cost_mse','vol_cost_cpdf1','vol_cost_cpdf2','-v7.3');
  
     % Write out desired native volumes and segmentations
     save(fullfile(outputMatDir,'vol_ct_ctx.mat'),'vol_ct_ctx','-v7.3');
     save(fullfile(outputMatDir,'vol_t2_ctx.mat'),'vol_t2_ctx','-v7.3');
     save(fullfile(outputMatDir,'vol_t2_seg_ctx.mat'),'vol_t2_seg_ctx','-v7.3');   
  
     % Write out desired registration output volumes and segmentations
     save(fullfile(outputMatDir,'vol_t2_res_ctx.mat'),'vol_t2_res_ctx','-v7.3');
     save(fullfile(outputMatDir,'vol_t2_seg_res_ctx.mat'),'vol_t2_seg_res_ctx','-v7.3');
     save(fullfile(outputMatDir,'volmu_t2_seg_res_ctx.mat'),'volmu_t2_seg_res_ctx','-v7.3');
     save(fullfile(outputMatDir,'volmu_ct_seg_res_ctx.mat'),'volmu_ct_seg_res_ctx','-v7.3');
  
     save(fullfile(outputMatDir,'volmu_t2_res_ctx.mat'),'volmu_t2_res_ctx','-v7.3');
     save(fullfile(outputMatDir,'volmu_ct_res_ctx.mat'),'volmu_ct_res_ctx','-v7.3');
 
     save(fullfile(outputMatDir,'vol_t2_res_up_ctx.mat'),'vol_t2_res_up_ctx','-v7.3');
 
     % Write out atlas volumes and segmentations
     save(fullfile(outputMatDir,'volmu_t2_ctx.mat'),'volmu_t2_ctx','-v7.3');
     save(fullfile(outputMatDir,'volmu_ct_ctx.mat'),'volmu_ct_ctx','-v7.3');
     save(fullfile(outputMatDir,'volmu_t2_seg_ctx.mat'),'volmu_t2_seg_ctx','-v7.3');
     save(fullfile(outputMatDir,'volmu_ct_seg_ctx.mat'),'volmu_ct_seg_ctx','-v7.3');
  
     % Write out desired midput -- transformations, cost, etc.
     save(fullfile(outputMatDir,'vol_ct_sca_res_ctx.mat'),'vol_ct_sca_res_ctx','-v7.3');
     save(fullfile(outputMatDir,'volsd_ct_ctx.mat'),'volsd_ct_ctx','-v7.3');
     save(fullfile(outputMatDir,'vol_t2_seg_sca_ctx.mat'),'vol_t2_seg_sca_ctx','-v7.3');
  
     save(fullfile(outputMatDir,'vol_t2_seg_ctr_ctx.mat'),'vol_t2_seg_ctr_ctx','-v7.3');
     save(fullfile(outputMatDir,'vol_t2_ctr_ctx.mat'),'vol_t2_ctr_ctx','-v7.3');
     save(fullfile(outputMatDir,'vol_ct_ctr_ctx.mat'),'vol_ct_ctr_ctx','-v7.3');
  
     save(fullfile(outputMatDir,'S_native_ct.mat'),'S_native_ct','-v7.3');
     save(fullfile(outputMatDir,'vol_cost_mse.mat'),'vol_cost_mse','-v7.3');
     save(fullfile(outputMatDir,'vol_cost_cpdf1.mat'),'vol_cost_cpdf1','-v7.3');
     save(fullfile(outputMatDir,'vol_cost_cpdf2.mat'),'vol_cost_cpdf2','-v7.3');
     save(fullfile(outputMatDir,'cost.mat'),'cost','-v7.3');
     
     save(fullfile(outputMatDir,'M_scale_t2.mat'),'M_scale_t2','-v7.3');
     save(fullfile(outputMatDir,'M_ct_atlas.mat'),'M_ct_atlas','-v7.3');
     save(fullfile(outputMatDir,'M_trans_t2.mat'),'M_trans_t2','-v7.3');
  end
 
  result.vol_ct_ctx = vol_ct_ctx; % Native Patient CT Image
  result.vol_t2_res_up_ctx = vol_t2_res_up_ctx; % Output Registered T2 Image
  result.vol_t2_seg_res_up_ctx = vol_t2_seg_res_up_ctx; % Output Registered T2 Prostate Segmentation 
  result.vol_cropped_up_ctx = vol_cropped_up_ctx; % Output Upsampled CT Image
  result.M_trans_t2 = M_trans_t2; % Transformation Vector for the T2 to T2 Atlas Registration
  result.M_scale_t2 = M_scale_t2; % Scaling Factor for the Native T2 Image
  result.M_ct_atlas = M_ct_atlas; % Transformation Vector for the CT to CT Atlas Registration 

  fprintf('%s -- %s.m:    Completed MR-CT Registration. Exiting PROMRCT Registration Pipeline...\n',datestr(now),mfilename);
  
  endtime = now;
  fprintf(1,'***Done*** %s (%0.2f mins)\n',datestr(endtime),(endtime-starttime)*60*24);
  
end
  
