% lab5_initial.m

% load a demonstration file (already in the matlab code base) that
% contains a 180x360 matrix with elevations above sea level (in meters)
load topo
% Make lat/long vectors with values corresponding to the rows/columns of 
% topo which range from -89.5 to 89.5, and 0.5 to 359.5, respectively.
lat =[-89.5:89.5];
long=[0.5:359.5];

% Display the topography, with axes labelled in lat/longs
figure(1); clf;
imagesc(long,lat,topo)
set(gca,'ydir','normal');
xlabel('Longitude');
ylabel('Latitude');

% Grab a single point by clicking with the mouse and find 
% the row and column indexes corresponding to the position 
% in topo nearest to the point clicked.
[X,Y]  = ginput(1);
[~,ix] = min(abs(long-X)); %the ~ is a dummy variable indicating we ...
[~,iy] = min(abs(lat-Y));  %  ...don't need this output variable
if ix == 360
    ixright = 1;
else
    ixright = ix+1;
end
if ix == 1
    ixleft = 360;
else
    ixleft = ix-1;
end
slope = atand((topo(iy,ixright)-topo(iy,ixleft))/((long(ixright)-long(ixleft))*111000*cos(long(ix))));

% Draw a small black square marker at that point, and then add a label 
% beside the marker with its height, latitude, and longitude.
line(long(ix),lat(iy),'color','k','marker','s','linewidth',3);
text(long(ix),lat(iy),{['height (m)='      num2str(topo(iy,ix))],...
                       [' at lat '         num2str(lat(iy))  ],...
                       [' and long '       num2str(long(ix)) ],...
                       [' east/west slope' num2str(slope)    ] });
                   
disp(strcat('lat=',num2str(lat(iy))))
disp(strcat('long=',num2str(long(ix))))
disp(strcat('height=',num2str(topo(iy,ix))))
disp(strcat('slope=',num2str(slope)))
if slope <= -0.1
    disp('East-facing')
else if slope >= 0.1
        disp('West-facing')
    else
        disp('Flat')
    end
end

partner.name = 'Kyle Zell & Connor Lyons';
Time_spent   = 3;