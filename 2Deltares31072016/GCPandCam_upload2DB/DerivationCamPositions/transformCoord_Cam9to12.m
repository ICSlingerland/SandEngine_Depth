% transformCoord_Cam9to12.m

% Estimate position of cam8 to cam12 in their UP position (=in top of the
% mast) by using average xyz offsets (DOWN vs UP position)for cam 1 to cam8
%
% NOTE we correct vertical location from measured z-coordinate on top of
% camera housing to a camera elevation, by using average elevation 
% difference of GPS surveyed DOWN elevation on top of camera housing and 
% camera elevation as in existing data base file for cam 1 to cam 8
% (so no need to correct elevation measurement in DOWN position to actual 
% camera elevation in UP position)

xyzCams_DOWNandUP = ...
   [-1.7884  -1.2135   15.0010  -1.8946     -1.0807     42.735; ...
   -0.6124   -0.7854   15.0010  -0.71339    -0.68998    42.749;...
    0.0058   -0.4690   15.0010  -0.11346    -0.37195    42.74;...
    0.5418   -0.1276   15.0090  0.41514     -0.049806   42.742;...
    0.3981    1.1669   15.0000  0.30185     1.2271      42.742;...
    0.0021    1.4443   14.9840  -0.088923   1.4906      42.735;...
   -1.2516    2.1232   14.9730  -1.3141     2.1961      42.734;...
   -1.8781    2.1184   14.9470  -1.9372     2.1889      42.739;...
   -2.2803    1.4310   14.9590      NaN         NaN       NaN ;...    
   -2.2948    0.7466   14.9430      NaN         NaN       NaN ;... 
   -2.3136    0.2074   15.0070      NaN         NaN       NaN ;... 
   -2.2697   -0.4672   14.9630      NaN         NaN       NaN ] 
      
% cam1 to cam8 : x y z frame down, new survey; x y z frame from DB
A= [-1.7884   -1.2135   15.0010  -1.8946     -1.0807     42.735; ...
   -0.6124   -0.7854   15.0010  -0.71339    -0.68998    42.749 ;...
    0.0058   -0.4690   15.0010  -0.11346    -0.37195    42.74 ;...
    0.5418   -0.1276   15.0090  0.41514     -0.049806   42.742;...
    0.3981    1.1669   15.0000  0.30185     1.2271      42.742;...
    0.0021    1.4443   14.9840  -0.088923   1.4906      42.735;...
   -1.2516    2.1232   14.9730  -1.3141     2.1961      42.734;...
   -1.8781    2.1184   14.9470  -1.9372     2.1889      42.739  ] 

% estimate UP position for cam8 to cam12 by using the average differences 
% for x,y,z for cam1 to cam8 (DOWN-UP)
% their is also a small offset in x,y direction (so not only for z)

xyzCams_DOWNandUPcalc = xyzCams_DOWNandUP;
xyzCams_DOWNandUPcalc(9:12,4) = xyzCams_DOWNandUPcalc(9:12,1)-(mean(A(:,1)-A(:,4)))
xyzCams_DOWNandUPcalc(9:12,5) = xyzCams_DOWNandUPcalc(9:12,2)-(mean(A(:,2)-A(:,5)))
xyzCams_DOWNandUPcalc(9:12,6) = xyzCams_DOWNandUPcalc(9:12,3)-(mean(A(:,3)-A(:,6)))

%add column with estimated elevation for ALL cameras (for plotting)
xyzCams_DOWNandUPcalc(:,7)= xyzCams_DOWNandUPcalc(:,3)-(mean(A(:,3)-A(:,6))) 
% for plotting purposes add cam 1 at bottom of list    
camsDOWN=xyzCams_DOWNandUP(1:12,1:3); camsDOWN(13,:)= camsDOWN(1,:)
camsUP=xyzCams_DOWNandUPcalc(1:12,4:6); camsUP(13,:)= camsUP(1,:)

% In addition to a translation, there also appears to be some tilting
% happening as the frame is raised to its UP position.
% Correct cam9-cam12 elevation for tilting of frame
% use zcam8-zcam1/ycam8-ycam1 as approximation (y=y_Argus)
z_corr_c1=camsUP(1,3)-camsDOWN(1,3)-27.75;
z_corr_c8=camsUP(8,3)-camsDOWN(8,3)-27.75;
y_corr_c1=camsDOWN(1,2)-(mean(A(:,2)-A(:,5)));
y_corr_c8=camsDOWN(8,2)-(mean(A(:,2)-A(:,5)));

S=(z_corr_c8-z_corr_c1)/(y_corr_c8-y_corr_c1); 

y_corr=camsDOWN(1:12,2)-(mean(A(:,2)-A(:,5)));

z_corr_c9 = xyzCams_DOWNandUPcalc(9,6)+S*(y_corr(9)-y_corr(1))+z_corr_c1;
z_corr_c10 = xyzCams_DOWNandUPcalc(10,6)+S*(y_corr(10)-y_corr(1))+z_corr_c1;
z_corr_c11 = xyzCams_DOWNandUPcalc(11,6)+S*(y_corr(11)-y_corr(1))+z_corr_c1;
z_corr_c12 = xyzCams_DOWNandUPcalc(12,6)+S*(y_corr(12)-y_corr(1))+z_corr_c1;
newcamsZ=[z_corr_c9 z_corr_c10 z_corr_c11 z_corr_c12]

newCamsUP=xyzCams_DOWNandUPcalc(:,4:6);
newCamsUP(9:12,3) = newcamsZ';
% for plotting purposes add cam1 at end of the file
newCamsUP(13,:)=newCamsUP(1,:);

figure(103),clf
hold on
    plot3(newCamsUP(:,2),newCamsUP(:,1),newCamsUP(:,3),'bo-')  
    plot3(newCamsUP(9:12,2),newCamsUP(9:12,1),newCamsUP(9:12,3),'b*-')  
    axis([-2 3 -3.5 1 42.5 43]) 
    set(gca, 'xdir','reverse')
    grid on
    title('camera positions')
    xlabel('y Argus (m)'), ylabel('x Argus (m)')
   
xyz_cam9tocam12=newCamsUP(9:12,:);