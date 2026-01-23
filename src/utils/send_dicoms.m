function send_dicoms(dicom_dir, destinations)

for d = 1:length(destinations)
  aet = destinations(d).ae_title;
  aec = destinations(d).called_ae_title;
  ip = destinations(d).ip_address;
  port = destinations(d).port;
  cmd = sprintf('/bin/bash -c ''storescu -v -nh -aet %s -aec %s +sd +r %s %s %s''', aet, aec, ip, port, dicom_dir);
  fprintf('%s\n', cmd);
  system(cmd);
end

end
