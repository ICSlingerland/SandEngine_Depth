% TempGCPs2DB.m
% script to upload temporary GCPs to data base for cam9 to cam12

% Note: names of gcp include the camera number and an id-number that relates 
% to a field sketch on which locations of temporary GCPs are indicated 
% Example: gcp0901 = gcp for camera 9, id-number 1 on field sketch

% first run argusInit so you can connect to DB

load 'gcpnames.mat'
load 'xyzArg_tempGCP.mat'
clear GCP % contains raw surveydata (in RD)
 
startseq=5100;
idstart=100;
 
for i=1:length(GCP_xyzArg)
    new_gcp=DBCreateEmptyStruct('gcp')
    
    id = ['ZMXX0',num2str(idstart+i)]
    name = gcpnames(i,:)
    
    new_gcp.seq = startseq+i
    new_gcp.id = id
    new_gcp.name = name
    new_gcp.siteID = 'ZMXXXXX'
    new_gcp.timeIN = matlab2Epoch(datenum(2016,6,17,01,0,0))
    new_gcp.timeOUT = matlab2Epoch(datenum(2016,6,17,18,00,0))
    new_gcp.x = GCP_xyzArg(i,1)
    new_gcp.y = GCP_xyzArg(i,2)
    new_gcp.z = GCP_xyzArg(i,3)
    new_gcp.xDim = 0
    new_gcp.yDim = 0
    new_gcp.intensity = 0
    new_gcp.eccentricity = 0
    new_gcp.PPPable = 0
    new_gcp.timestamp = datestr(now)
    
    DBInsert('gcp',new_gcp);
    clear new_gcp
end    
