function paths = fetch_paths_pro_ignite(exam_dir)

% Enumerate data to be collected -----------------------------
paths = struct;

paths.RSI_raw = {''};
paths.T2_ax = {''};
paths.CT = {''};
paths.RT = {''};
paths.generic = {''};


% Patterns used to ID data -----------------------------------
% Ignore
global_exclude = {'.', '..', '.DS_Store'};

% RSI general
patterns_RSI_raw = {'RSI'};
patterns_exclude_RSI_raw = {'trace', 'adc', 'color', 'fa', 'tensor', '[\s_]+T2[\s_]+', 'NOT.FOR.CLINICAL.USE', 'Apparent.Diffusion.Coefficient', 'Restricted.Signal.Map', 'RSI_anatomic_T2W', 'RSI_DWI_averages', 'RSI_C', 'RSIrs_Experimental', 'RSI_Visual_Report'};

% RSI GE 
seqID_RSI_raw_GE = {'epi2_pepolarFOCUSFLEX', 'epi2_pepolarFLEX', 'epi2_ART', 'epi2_revART', 'epi2alt', 'epi2altoff'};

% RSI Siemens
patterns_seqID_RSI_raw_Siemens = {'ep_b', 'ez_b', 'WI_b'};

% RSI Philips
patterns_seqID_RSI_raw_philips = {'DwiSE'};

% Axial T2 for RSI overlay
patterns_T2_ax_plane = {'ax', 'tra'};
patterns_T2_ax_contrast = {'T2', 'FSE'};
patterns_exclude_T2_ax = {'water', 'fat', 'flex', 'sag', 'cor', 'reformat', 'cube', 'T1', 'bifurcation', 'prop', '3d', 'space'};


% -----------------------------------------------------------
RSI_raw_path_num = 1;
T2_ax_path_num = 1;
CT_path_num = 1;
RT_path_num = 1;
generic_path_num = 1;

acqs = dir(exam_dir);
for i = 1:length(acqs)
  if ~any(strcmp(acqs(i).name, global_exclude))
    acq_path = fullfile(exam_dir, acqs(i).name);

    ims = dir(acq_path);
    try
      im = ims(3).name;
      info = dicominfo(fullfile(acq_path, im));
      info = fix_impax_dcm_tags(info);
    catch ME
      fprintf('%s\n', ME.message)
      continue
    end

    if ~isfield(info, 'ImageType')
       continue
    end
    im_type = info.ImageType;
    
    SeriesDescription = info.SeriesDescription;
    Modality = info.Modality;
    
    % Sort out vendor-specific stuff -----------------------------------------------
    manufacturer = info.Manufacturer;
    if strcmpi(manufacturer, 'ge medical systems')
      manufacturer = 'ge';
    elseif any(strcmpi(manufacturer, {'siemens', 'siemens healthineers'}))
      manufacturer = 'siemens';
    elseif any(strcmpi(manufacturer, {'philips', 'philips healthcare'}))
      manufacturer = 'philips';
    end

    seq_name = '';
    switch manufacturer
      case 'ge'
	if isfield(info, 'Private_0019_109c')
	  seq_name = info.Private_0019_109c;
	end

      case 'siemens'
	if isfield(info, 'SequenceName')
	  seq_name = info.SequenceName;
	end

      case 'philips'
	if isfield(info, 'MRSeriesScanningTechniqueDesc')
	  seq_name = info.MRSeriesScanningTechniqueDesc;
	end

      otherwise
	seq_name = '';
    end

    % Check for RSI data -----------------------------------------------------
    is_RSI_raw = 0;
    if strcmp(manufacturer, 'ge')
      match_RSI_raw_seq = any(~cellfun(@isempty, regexpi(seq_name, seqID_RSI_raw_GE)));
      match_RSI_raw_name = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_RSI_raw)));
      match_exclude = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_exclude_RSI_raw)));
      if (match_RSI_raw_seq || match_RSI_raw_name) && ~match_exclude
	paths.RSI_raw{RSI_raw_path_num} = acq_path;
	series_description_list{RSI_raw_path_num} = SeriesDescription; % Save for later filtering of reverse acquisitions
	RSI_raw_path_num = RSI_raw_path_num + 1;
	is_RSI_raw = 1;
      end
    end
    if strcmp(manufacturer, 'siemens')
      match_RSI_raw_seq = any(~cellfun(@isempty, regexpi(seq_name, patterns_seqID_RSI_raw_Siemens)));
      match_RSI_raw_name = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_RSI_raw)));
      match_exclude = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_exclude_RSI_raw)));
      if (match_RSI_raw_seq || match_RSI_raw_name) && ~match_exclude
	paths.RSI_raw{RSI_raw_path_num} = acq_path;
	series_description_list{RSI_raw_path_num} = SeriesDescription; % Save for later filtering of reverse acquisitions
	RSI_raw_path_num = RSI_raw_path_num + 1;
	is_RSI_raw = 1;
      end
    end
    if strcmp(manufacturer, 'philips')
      match_RSI_raw_seq = any(~cellfun(@isempty, regexpi(seq_name, patterns_seqID_RSI_raw_philips)));
      match_RSI_raw_name = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_RSI_raw)));
      match_exclude = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_exclude_RSI_raw)));
      if (match_RSI_raw_seq || match_RSI_raw_name) && ~match_exclude
	paths.RSI{RSI_raw_path_num} = acq_path;
	series_description_list{RSI_raw_path_num} = SeriesDescription; % Save for later filtering of reverse acquisitions
	RSI_raw_path_num = RSI_raw_path_num + 1;
	is_RSI_raw = 1;
      end
    end

    % Check for anatomical T2 DICOMs
    is_T2 = 0;
    match_T2_plane = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_T2_ax_plane)));
    match_T2_contrast = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_T2_ax_contrast)));
    match_exclude = any(~cellfun(@isempty, regexpi(SeriesDescription, patterns_exclude_T2_ax)));
    match_exclude_reformat = ~isempty(regexpi(im_type, 'REFORMAT'));
    match_app_output = ~isempty(regexpi(SeriesDescription, 'RSI_anatomic_T2W'));
    if match_app_output || (match_T2_plane && match_T2_contrast && ~match_exclude && ~match_exclude_reformat)
      paths.T2_ax{T2_ax_path_num} = acq_path;
      T2_ax_path_num = T2_ax_path_num + 1;
      is_T2 = 1;
    end

    % Check for CT DICOMs
    is_CT = 0;
    match_CT = strcmp(Modality, 'CT');
    if match_CT
      paths.CT{CT_path_num} = acq_path;
      CT_path_num = CT_path_num + 1;
      is_CT = 1;
    end

    % Check for RT DICOMs
    is_RT = 0;
    match_RT = strcmp(Modality, 'RTSTRUCT');
    if match_RT
       paths.RT{RT_path_num} = acq_path;
       RT_path_num = RT_path_num + 1;
       is_RT = 1;
    end

    % If none of the above, put in list of "generic" scans
    if ~is_RSI_raw && ~is_T2 && ~is_CT && ~is_RT
      paths.generic{generic_path_num} = acq_path;
      generic_path_num = generic_path_num + 1;
    end

  end
end


end
