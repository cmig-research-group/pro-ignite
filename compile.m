entry_point = '/home/ccconlin/apps/pro-ignite/pro_ignite.m';
exe_name = 'pro_ignite';
destination_folder = '/home/ccconlin/apps/pro-ignite/compiled';

string_add = '';
paths_add = {
'/home/ccconlin/code/rsi_pelvis/autoseg'
'/home/ccconlin/code/cmig_core_utils/dicom/dicts'
'/home/ccconlin/code/cmig_core_utils/gradient_nonlinearity'
'/home/ccconlin/code/rsi_pelvis/protocols/artpro_protocol_reference.mat'
'/home/ccconlin/apps/pro-ignite/share'
'/home/ccconlin/apps/pro-ignite/src'
};
for i = 1:length(paths_add)
  string_add = [string_add '-a ' paths_add{i} ' '];
end

cmd = sprintf('mcc -v -o %s -W %s -T link:exe -d %s %s %s' , exe_name, ['main:' exe_name], destination_folder, entry_point, string_add);
fprintf('%s\n', cmd);
eval(cmd);
