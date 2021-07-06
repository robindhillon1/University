% lab5_initial.m

% load a demonstration file (already in the matlab code base) that
% contains a 180x360 matrix with elevations above sea level (in meters)
load topo
% Make lat/long vectors with values corresponding to the rows/columns of 
% topo which range from -89.5 to 89.5, and 0.5 to 359.5, respectively.
lat = -89.5:89.5;
long= 0.5:359.5;

% Display the topography, with axes labelled in lat/longs
clf;
figure(1); 
imagesc(long,lat,topo)
set(gca,'ydir','normal');
xlabel('Longitude');
ylabel('Latitude');
                   
% Grab a single point by clicking with the mouse and find 
% the row and column indexes corresponding to the position 
% in topo nearest to the point clicked.
[X,Y]=ginput(1);
[~,ix]=min(abs(long-X)); %the ~ is a dummy variable indicating we ...
[~,iy]=min(abs(lat-Y));  %  ...don't need this output variable

% Draw a small black square marker at that point, and then add a label 
% beside the marker with its height, latitude, and longitude.
line(long(ix),lat(iy),'color','k','marker','s','linewidth',3);
text(long(ix),lat(iy),{['height (m)='    num2str(topo(iy,ix))],...
                       [' at lat '       num2str(lat(iy))  ],...
                       [' and long '     num2str(long(ix)) ]});%%%%%


dHeight = topo(iy,ix+1)-topo(iy,ix-1);
dLong = long(ix+1)-long(ix-1);
slope = atand((dHeight)/(dLong*111000*cosd(lat(iy))));




if slope >= -0.1 && slope <= 0.1 
    disp('The slope of the land is: Flat');
elseif slope < -0.1 && slope > -1
    disp(['The slope of the land is: East Facing ' num2str(slope)]);
elseif slope > 0.1 && slope < 1
    disp(['The slope of the land is: West Facing ' num2str(slope)]);
else %%ix >= 359.5 || ix <= 1
    disp('The slope of the land is out of the DEM ');
end


                   
                   
                   