function writeSRdicom(dicomdir, outputdir, message_text, fixed_uids)
% dicomdir - template dicom
% outputdir - directory to save created dicoms
    
load('srTemplate.mat');

if exist('fixed_uids', 'var')
  uid.sop = fixed_uids.SOPInstanceUID;
  uid.series = fixed_uids.SeriesInstanceUID;
else
  uid.sop = dicomuid;
  uid.series = dicomuid;
end

mkdir(outputdir);
if isfolder(dicomdir)
  % If dicomdir is a directory, retrieve metadata from the first file in the directory
  filelist = dir(dicomdir);
  filename = fullfile(filelist(3).folder, filelist(3).name);
else
  % If dicomdir is a file, use that file to retrieve metadata
  filename = dicomdir;
end
pt_metadata = dicominfo(filename);

seriesDesc = sprintf('PPro Status Report');

% Change relevant tags of interest
srTemplate.SeriesDescription = seriesDesc;
srTemplate.StudyInstanceUID = pt_metadata.StudyInstanceUID;
srTemplate.StudyID = pt_metadata.StudyID;
srTemplate.StudyDate = datestr(now, 'yyyymmdd');
srTemplate.StudyTime = datestr(now, 'HHMMSS');
srTemplate.SeriesDate = datestr(now, 'yyyymmdd');

srTemplate.MediaStorageSOPInstanceUID = uid.sop;
srTemplate.SOPInstanceUID = uid.sop;
srTemplate.SeriesInstanceUID = uid.series;

srTemplate.PatientID = pt_metadata.PatientID;
srTemplate.PatientName = pt_metadata.PatientName;

% Alter relevant text entries in Report
currentDateTime = datetime('now');
formattedDateTime = char(currentDateTime);

srTemplate.ContentSequence.Item_6.UID = uid.series;

srTemplate.ContentSequence.Item_12.ConceptNameCodeSequence.Item_1.CodeMeaning = 'Precision Pro Status Report';
srTemplate.ContentSequence.Item_12.ContentSequence.Item_1.TextValue = sprintf('%s', message_text);
srTemplate.ContentSequence.Item_12.ContentSequence.Item_2.TextValue = sprintf('Content last updated: %s \n ', formattedDateTime);

% write to DICOM
dicomwrite([], sprintf('%s/SR.dcm', outputdir), srTemplate, 'CreateMode', 'Copy'); 

end
    
