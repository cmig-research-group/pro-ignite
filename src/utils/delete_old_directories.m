function delete_old_directories(path_in, age_limit)

% age_limit is specified as an integer number of hours
% e.g., 12, 24, 48 

dir_in = dir(path_in);

t_now = datetime('now');

for i = 1:length(dir_in)
  if strcmp(dir_in(i).name, '.') | strcmp(dir_in(i).name, '..')
    continue
  end

  t_mod = datetime(dir_in(i).date, 'InputFormat', 'dd-MMM-yyy HH:mm:ss');
  t_diff = t_now - t_mod;

  if hours(t_diff) > age_limit
    path2delete = fullfile(dir_in(i).folder, dir_in(i).name);
    fprintf('Deleting expired directory: %s\n', path2delete);
    rmdir(path2delete, 's');
  end

end

end
