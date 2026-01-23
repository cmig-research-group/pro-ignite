function pro_ignite(path_input, path_config_file)

params = eval_file(path_config_file);


% Directory management -----------------------------------------------------------
path_tmp = params.tmp_directory;
if ~exist(path_tmp, 'dir')
  mkdir(path_tmp);
end
fprintf('Temporary directory for incoming data: %s\n', path_tmp);
path_tmp_raw = fullfile(path_tmp, 'raw');
if ~exist(path_tmp_raw, 'dir')
  mkdir(path_tmp_raw);
end
path_tmp_proc = fullfile(path_tmp, 'proc');
if ~exist(path_tmp_proc, 'dir')
  mkdir(path_tmp_proc);
end

% Delete old directories in tmp directory
tmp_age_limit = params.lifetime_tmp_directory;
delete_old_directories(path_tmp_raw, tmp_age_limit);
delete_old_directories(path_tmp_proc, tmp_age_limit);


% Sort incoming data and collect patient ID
data_info = unpack_data(path_input, path_tmp_raw, 'PatientID');
subj_id = data_info.subject_ids{1};
scan_dates = data_info.scan_dates;


% Fetch paths to relevant series
paths = struct('RSI_raw', [], 'T2_ax', [], 'CT', [], 'RT', [], 'generic', []);
for i = 1:length(scan_dates)
  paths_date = fetch_paths_pro_ignite( fullfile(path_tmp_raw, subj_id, scan_dates{i}) );
  paths.RSI_raw = cat(2, paths.RSI_raw, paths_date.RSI_raw);
  paths.T2_ax = cat(2, paths.T2_ax, paths_date.T2_ax);
  paths.CT = cat(2, paths.CT, paths_date.CT);
  paths.RT = cat(2, paths.RT, paths_date.RT);
  paths.generic = cat(2, paths.generic, paths_date.generic);
end
paths.RSI_raw = paths.RSI_raw(~cellfun(@isempty, paths.RSI_raw));
paths.T2_ax = paths.T2_ax(~cellfun(@isempty, paths.T2_ax));
paths.CT = paths.CT(~cellfun(@isempty, paths.CT));
paths.RT = paths.RT(~cellfun(@isempty, paths.RT));
paths.generic = paths.generic(~cellfun(@isempty, paths.generic));

disp(paths);
fprintf('\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                Processing                                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path_sr_out = fullfile(path_tmp_proc, subj_id, 'structured_reports');

% Process T2 ----------------------------------------------------------------------
if ~isempty(paths.T2_ax)

  path_in_T2 = paths.T2_ax{1};
  path_out_T2 = strrep(path_in_T2, '/raw/', '/proc/');

  t_now = datetime('now');
  message_text = sprintf('%s -- T2 volume received for processing\n(T2 processing may take upwards of 60 minutes)\n\n', datestr(t_now, 'HH:MM:SS'));
  writeSRdicom(path_in_T2, path_sr_out, message_text);
  send_dicoms(path_sr_out, destinations);

  fname_contour_prostate_mgz = fullfile(path_out_T2, 'prostate_contour_T2_space.mgz');
  fname_T2_GUW_mgz = fullfile(path_out_T2, 'T2_corrected_GUW.mgz');
  fname_contour_urethra_mgz = fullfile(path_out_T2, 'urethra_contour_T2_space.mgz');

  if ~exist(fname_contour_prostate_mgz, 'file')
    fprintf('Processing T2 volume: %s\n', path_in_T2);
    params_T2 = struct;
    params_T2.ProstateSegContainer = params.container;
    process_T2(path_in_T2, path_out_T2, params_T2);
  else
    fprintf('T2 volume already processed: %s\n', fname_contour_prostate_mgz);
  end

  path_out_T2_dcms = fullfile(path_out_T2, 'DICOM_Output', 'PPro_T2W');
  files_T2 = dir(path_out_T2_dcms);
  dcminfo_T2 = dicominfo(fullfile(files_T2(3).folder, files_T2(3).name));

  paths.RT{end+1} = fullfile(path_out_T2, 'DICOM_Output', 'RT_DICOM');

  t_now = datetime('now');
  message_text = sprintf('%s%s -- T2 processing complete\n\n', message_text, datestr(t_now, 'HH:MM:SS'));
  writeSRdicom(path_in_T2, path_sr_out, message_text);
  send_dicoms(path_sr_out, destinations);

end


% Process RSI ---------------------------------------------------------------------
if ~isempty(paths.RSI_raw)

  t_now = datetime('now');
  message_text = sprintf('%s%s -- RSI data received for processing\n(RSI processing may take upwards of 60 minutes)\n\n', message_text, datestr(t_now, 'HH:MM:SS'));
  writeSRdicom(path_in_T2, path_sr_out, message_text);
  send_dicoms(path_sr_out, destinations);

  fprintf('Processing raw RSI series\n');
  
  % Pass entire patient directory for RSI processing
  % RSI pipeline code will handle the rest
  path_in_rsi = fullfile(path_tmp_raw, subj_id); 
  path_out_rsi = path_tmp_proc;
  rsi_output_dir_list = RSI_pipeline(path_in_rsi, path_out_rsi, params.rsi_config_file);

  for i = 1:length(rsi_output_dir_list)
    path_rsi_out = rsi_output_dir_list{i};
    path_rsi_out_dcms = fullfile(path_rsi_out, 'DICOM_Output');
    path_rsi_out_dcms = dir(path_rsi_out_dcms);
    for j = 3:length(path_rsi_out_dcms)
      if strcmp(path_rsi_out_dcms(j).name, 'RSI_Visual_Report')
	continue
      elseif strcmp(path_rsi_out_dcms(j).name, 'RT_DICOM')
	paths.RT{end+1} = fullfile(path_rsi_out_dcms(j).folder, path_rsi_out_dcms(j).name);
      else
	paths.generic{end+1} = fullfile(path_rsi_out_dcms(j).folder, path_rsi_out_dcms(j).name);
      end
    end
  end

  t_now = datetime('now');
  message_text = sprintf('%s%s -- RSI processing complete\n\n', message_text, datestr(t_now, 'HH:MM:SS'));
  writeSRdicom(path_in_T2, path_sr_out, message_text);
  send_dicoms(path_sr_out, destinations);

end


% MR-CT registration -------------------------------------------------------------
if ~isempty(paths.CT) && ~isempty(paths.T2_ax)

  t_now = datetime('now');
  message_text = sprintf('%s%s -- MR-CT registration underway\n(MR-CT registration may take upwards of 20 minutes)\n\n', message_text, datestr(t_now, 'HH:MM:SS'));
  writeSRdicom(path_in_T2, path_sr_out, message_text);
  send_dicoms(path_sr_out, destinations);

  path_in_CT = paths.CT{1};
  path_out_reg = strrep(path_in_CT, '/raw/', '/proc/');
  path_out_reg = nixify([path_out_reg '_reg']);
  path_out_reg_dcms = fullfile(path_out_reg, 'DICOM_Output');

  fprintf('Performing MR-CT registration\n');
  result = PROMRCT_App(path_in_CT, path_out_T2_dcms, path_out_reg, 'ProstateContourPath', fname_contour_prostate_mgz);
  M_T2_to_CT = result.M_reg;
  ctx_T2_CT = result.vol_t2_res_up_ctx; % Output Registered T2 Image
  ctx_contour_CT = result.vol_t2_seg_res_up_ctx; % Output Registered T2 Prostate Segmentation 

  % Write DICOM files for T2 outputs
  fnames_CT = dir(path_in_CT);
  dcminfo_CT = dicominfo(fullfile(fnames_CT(3).folder, fnames_CT(3).name));
  date_CT = dcminfo_CT.StudyDate;
  datetime_CT = datetime(date_CT, 'InputFormat', 'yyyyMMdd', 'Format', 'yyyy.MM.dd');

  SeriesDescription_T2 = dcminfo_T2.SeriesDescription;
  date_T2 = dcminfo_T2.StudyDate;
  datetime_T2 = datetime(date_T2, 'InputFormat', 'yyyyMMdd', 'Format', 'yyyy.MM.dd');

  dcm_hdr_struct = dcminfo_CT;
  dcm_hdr_struct.Modality = 'MR';
  dcm_hdr_struct.SeriesDescription = sprintf('PPro Fusion | %s %s <> CT %s', SeriesDescription_T2, datetime_T2, datetime_CT);
  SeriesNumber_T2 = num2str(dcminfo_T2.SeriesNumber);
  SeriesNumber_T2_new = str2double([SeriesNumber_T2 '10']);
  dcm_hdr_struct.SeriesNumber = SeriesNumber_T2_new;
  dcm_hdr_struct.ScanningSequence = dcminfo_T2.ScanningSequence;
  dcm_hdr_struct.SequenceVariant = dcminfo_T2.SequenceVariant;
  dcm_hdr_struct.MRAcquisitionType = dcminfo_T2.MRAcquisitionType;
  dcm_hdr_struct.EchoTime = dcminfo_T2.EchoTime;

  path_out_T2_reg = fullfile(path_out_reg_dcms, sprintf('%s_reg', nixify(SeriesDescription_T2)));
  dicomwrite_cmig(ctx_T2_CT, path_out_T2_reg, dcm_hdr_struct);

  % Apply registration to every series in the 'generic' list
  if ~isempty(paths.generic)
    ctx_T2_orig = QD_ctx_load_mgh(fname_T2_GUW_mgz);
  end
  for i = 1:length(paths.generic)
      path_in_generic = paths.generic{i};
      ctx_gen = QD_read_dicomdir(path_in_generic);
      ctx_gen_T2 = vol_resample(ctx_gen, ctx_T2_orig, eye(4));
      ctx_gen_CT = vol_resample(ctx_gen_T2, ctx_T2_CT, M_T2_to_CT, 2);

      fnames_generic = dir(path_in_generic);
      dcminfo_generic = dicominfo(fullfile(fnames_generic(3).folder, fnames_generic(3).name));
      SeriesNumber_generic = num2str(dcminfo_generic.SeriesNumber);
      SeriesNumber_generic_new = str2double([SeriesNumber_generic '10']);
      dcm_hdr_struct.SeriesNumber = SeriesNumber_generic_new;
      SeriesDescription_generic = dcminfo_generic.SeriesDescription;
      date_generic = dcminfo_generic.StudyDate;
      datetime_generic = datetime(date_generic, 'InputFormat', 'yyyyMMdd', 'Format', 'yyyy.MM.dd');
      dcm_hdr_struct.SeriesDescription = sprintf('PPro Fusion | %s %s <> CT %s', SeriesDescription_generic, datetime_generic, datetime_CT);

      dcm_hdr_struct.ScanningSequence = dcminfo_generic.ScanningSequence;
      dcm_hdr_struct.SequenceVariant = dcminfo_generic.SequenceVariant;
      dcm_hdr_struct.MRAcquisitionType = dcminfo_generic.MRAcquisitionType;
      dcm_hdr_struct.EchoTime = dcminfo_generic.EchoTime;

      path_out_gen = fullfile(path_out_reg_dcms, sprintf('%s_reg', nixify(SeriesDescription_generic)));
      dicomwrite_cmig(ctx_gen_CT, path_out_gen, dcm_hdr_struct);
  end

  % Apply registration to every RT contour
  if ~isempty(paths.RT)
    ctx_T2_orig = QD_ctx_load_mgh(fname_T2_GUW_mgz);
  end

  candidate_series = [paths.T2_ax paths.CT, paths.generic];
  label_list = {};
  seg_indx = 1;
  for i = 1:length(paths.RT)
    path_RT = paths.RT{i};
    rtlist = recursive_dir(path_RT);
    for j = 1:length(rtlist)
      fname_RT = rtlist{j};
      rtinfo = dicominfo(fname_RT);
      ref_uid = rtinfo.ReferencedFrameOfReferenceSequence.Item_1.RTReferencedStudySequence.Item_1.ReferencedSOPInstanceUID;
      path_ref_dcms = RT_find_ref_series(ref_uid, candidate_series);

      ctx_rt = loadSegmentation(fname_RT, path_ref_dcms);
      rt_labels = ctx_rt.labels;

      ctx_rt_T2 = vol_resample(ctx_rt, ctx_T2_orig, eye(4));
      ctx_rt_CT = vol_resample(ctx_rt_T2, ctx_T2_CT, M_T2_to_CT, 2);
      
      for s = 1:length(rt_labels)
	if any(strcmp(rt_labels{s},label_list))
	  continue
	end
	segSTRUCT(seg_indx).number = seg_indx;
	segSTRUCT(seg_indx).name = rt_labels{s};
	segSTRUCT(seg_indx).seg = ctx_rt_CT.imgs(:,:,:,s);
	seg_indx = seg_indx + 1;
	label_list = [label_list rt_labels{s}];
      end

    end
  end

  path_out_rt = fullfile(path_out_reg_dcms, 'RT_DICOM');
  RK_write_segSTRUCT(segSTRUCT, path_out_T2_reg, path_out_rt, 'Prostate_contours', 963);

  t_now = datetime('now');
  message_text = sprintf('%s%s -- MR-CT registration complete\n\n', message_text, datestr(t_now, 'HH:MM:SS'));
  writeSRdicom(path_in_T2, path_sr_out, message_text);
  send_dicoms(path_sr_out, destinations);

end


% DICOM routing -------------------------------------------------------------------
destinations = eval_file(params.routing_config_file, 'destinations');

if exist(path_out_reg_dcms, 'dir')
  send_dicoms(path_out_reg_dcms, destinations);
else
  send_dicoms(fullfile(path_tmp_proc, subj_id), destinations);
end

t_now = datetime('now');
message_text = sprintf('%s%s -- Entire Pro-IGNITE pipeline complete', message_text, datestr(t_now, 'HH:MM:SS'));
writeSRdicom(path_in_T2, path_sr_out, message_text);
send_dicoms(path_sr_out, destinations);


end
