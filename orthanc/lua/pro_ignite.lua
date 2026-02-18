-- Change this path! It needs to point to the directory where config_app.m is stored.
PATH_CONFIG_DIR = '/home/ccconlin/work/pro-ignite'

-- Only change these paths if you aren't happy with the default configuration
PATH_WRITE_TARGET = '/home/orthanc/pro_ignite'
PATH_LOGS = '/home/orthanc/pro_ignite_logs'
PATH_PROIGNITE_TMP = '/home/tmp_pro_ignite'


function ToAscii(s)
   -- http://www.lua.org/manual/5.1/manual.html#pdf-string.gsub
   -- https://groups.google.com/d/msg/orthanc-users/qMLgkEmwwPI/6jRpCrlgBwAJ
   return s:gsub('[^a-zA-Z0-9-/-: ]', '_')
end


function WritePatient(patientId, tags, metadata)
   print('This patient is now stable, writing its instances on the disk: ' .. patientId)

   local studies = ParseJson(RestApiGet('/patients/' .. patientId)) ['Studies']

   for i, study in pairs(studies) do
      local series = ParseJson(RestApiGet('/studies/' .. study)) ['Series']

      for j, serie in pairs(series) do
	 local instances = ParseJson(RestApiGet('/series/' .. serie)) ['Instances']

	 for k, instance in pairs(instances) do
	    local dicom = RestApiGet('/instances/' .. instance .. '/file')
	    local path = ToAscii(PATH_WRITE_TARGET .. '/' .. patientId)
	    os.execute('mkdir -p "' .. path .. '"')
	    local target = assert(io.open(path .. '/' .. instance .. '.dcm', 'wb'))
	    target:write(dicom)
	    target:close()
	 end

      end
      
   end
   
end


function OnStablePatient(patientId, tags, metadata)

   WritePatient(patientId, tags, metadata)

   -- Delete so that OnStablePatient is triggered again when re-sent from modality
   RestApiDelete('/patients/' .. patientId) 

   os.execute(string.format('mkdir -p %s', PATH_PROIGNITE_TMP))
   os.execute(string.format('mkdir -p %s', PATH_LOGS))
   
   local path_in = ToAscii(PATH_WRITE_TARGET .. '/' .. patientId)

   local cmd = string.format('sudo docker run -v /var/run/docker.sock:/var/run/docker.sock -v %s:%s -v %s:%s -v %s:%s ghcr.io/cmig-research-group/pro_ignite %s %s/config_app.m >%s/%s.log 2>&1 &', PATH_WRITE_TARGET, PATH_WRITE_TARGET, PATH_PROIGNITE_TMP, PATH_PROIGNITE_TMP, PATH_CONFIG_DIR, PATH_CONFIG_DIR, path_in, PATH_CONFIG_DIR, PATH_LOGS, patientId)
   print(cmd)
   os.execute(cmd)
   
end
