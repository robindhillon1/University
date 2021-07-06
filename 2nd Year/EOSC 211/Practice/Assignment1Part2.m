clear
clf
% load CTD
% load cstlne



for k=1:932; %specific point selected for testing that is applicable to all matrices
dayofyear=(CTD.mtime(k) - datenum(CTD.gtime(1,k),1,0)); %locates a specific date and time from k value
xtk=datenum(0,1:2:12,1); %later used to convert x values into monthly names
[dist1,i1]=min(abs(CTD.pr(:,k)-2));
[dist2,i2]=min(abs(CTD.pr(:,k)-20));
[dist3,i3]=min(abs(CTD.pr(:,k)-200));
[dist4,i4]=min(abs(CTD.pr(:,k)-350));

%Calculations
CTD.location(k)=0;%Decision to plot based on latitude and longitude
plotFlag=0;
if(CTD.latitude(k)<=48.8) %Creates a specific area surrounding Juan de Fuca
        if (CTD.longitude(k)<=-123.6)
                colour = 'b';
                ptType = 'x';
                plotFlag = 'JDF';
                CTD.location(k) = 2;
                %colourPoint='bx';
        else
                colour = 'g';
                ptType = '.';
                
                %colourPoint='g.';  
        end
        
elseif (CTD.latitude(k)>=48.8)  %Creates a specific area surrounding Strait of Georgia
        if CTD.longitude(k)
                colour = 'r';
                ptType = 'x';
                plotFlag = 'SOG';
                CTD.location(k) = 1;
                %colourPoint='rx';
         
        else
                colour = 'g';    %if not within ranges, all outide points become green dots
                ptType = '.';               
                %colourPoint='g.';  
        end
else
        colour = 'g';   %All values not specified become a green dot
        ptType = '.';
        %colourPoint='g.';
end


%outputs
%figure


subplot(2,2,1) %first subplot
plot(cstlne(:,1),cstlne(:,2)) %plots the coastline figure
hold on
plot(CTD.longitude(k),CTD.latitude(k),strcat(colour,ptType)) %places dot at K value and specifies colour and shape
hold on

subplot(2,2,2) %second subplot location
plot(dayofyear,CTD.temp(k)) %plots day of year vs temperature
if (plotFlag == 'SOG') %if within Strait of Georgia range, values will be plotted, otherwise graph stays empty
        if dist1<=3,
        plot(dayofyear,CTD.temp(i1,k),'ro');
        hold on
        else NaN;
        end
        if dist2<=3,
        plot(dayofyear,CTD.temp(i2,k),'bx');
        hold on
        else NaN;
        end
        if dist3<=3,
        plot(dayofyear,CTD.temp(i3,k),'gd');
        else NaN;
        end
        if dist4<=3,
        plot(dayofyear,CTD.temp(i4,k),'ms');
        hold on
        else NaN;
        end      
end
hold on


subplot(2,2,3) %Third Subplot
plot(CTD.temp(:,k),CTD.pr(:,k)); %Plotting temperature vs depth
plot(CTD.temp(:,k),CTD.pr(:,k),strcat(colour, '-')) %Changes colour of line relative to location
hold on

subplot(2,2,4) %Subplot 4
plot(dayofyear,CTD.temp(k)) %Plotting day of year vs temperature
if (plotFlag == 'JDF') %If K value is within Juan de Fuca range, points are plotted. Otherwise graph remains clear
    if dist1<=3,
        plot(dayofyear,CTD.temp(i1,k),'ro');
        hold on
        else NaN;
        end
        if dist2<=3,
        plot(dayofyear,CTD.temp(i2,k),'bx');
        hold on
        else NaN;
        end
        if dist3<=3,
        plot(dayofyear,CTD.temp(i3,k),'gd');
        hold on
        else NaN;
        end
        if dist4<=3,
        plot(dayofyear,CTD.temp(350,k),'ms');
        hold on
        else NaN;
        end      
end
hold on
tempstorage(1,k) = CTD.temp(i1,k);
tempstorage(2,k) = CTD.temp(i2,k);
tempstorage(3,k) = CTD.temp(i3,k);
tempstorage(4,k) = CTD.temp(i4,k);
end

subplot(2,2,1) %first subplot
xlabel('Longitude');
ylabel('Latitude');
title('Station Locations')
hold on

subplot(2,2,2) %second subplot location
set(gca,'xtick',xtk,'xticklabel',{'Jan','Mar','May','July','Sep','Nov'}) %renames x-axis and values as months
xlim([1,max(xtk)])
ylim([4,25])
xlabel('Day of Year');
ylabel('Temperature ^oC')
title('Strait of Georgia Temperatures')
hold on

subplot(2,2,3) %Third Subplot
set(gca,'YDir','Reverse') %Reverses y-axis so that depth is grows downward
xlim([4,25])
xlabel('Temperature ^oC')
ylabel('Depth m')
title('Temperature Profiles')
hold on

subplot(2,2,4) %Subplot 4
set(gca,'xtick',xtk,'xticklabel',{'Jan','Mar','May','July','Sep','Nov'}) %renames x-axis and values as months
xlim([1,max(xtk)])
ylim([4,25])
xlabel('Day of Year');
ylabel('Temperature ^oC');
title('Juan de Fuca Strait Temperature');
hold on

for i=1:12;
      
%straight of georgia - average temp
tavg(1,i)=nanmean(tempstorage(1,CTD.gtime(2,:)==i & CTD.location==1));
tavg(2,i)=nanmean(tempstorage(2,CTD.gtime(2,:)==i & CTD.location==1));
tavg(3,i)=nanmean(tempstorage(3,CTD.gtime(2,:)==i & CTD.location==1));
tavg(4,i)=nanmean(tempstorage(4,CTD.gtime(2,:)==i & CTD.location==1));
%juan de fuca straight - average temp
tavg(5,i)=nanmean(tempstorage(1,CTD.gtime(2,:)==i & CTD.location==2));
tavg(6,i)=nanmean(tempstorage(2,CTD.gtime(2,:)==i & CTD.location==2));
tavg(7,i)=nanmean(tempstorage(3,CTD.gtime(2,:)==i & CTD.location==2));
tavg(8,i)=nanmean(tempstorage(4,CTD.gtime(2,:)==i & CTD.location==2));
 
 
end      
      
subplot(2,2,2) %plot average temp for strait of georgia (subplot c)
plot(datenum(0,1:12,15),tavg(1,:),'color','r','linewi',3);
plot(datenum(0,1:12,15),tavg(2,:),'color','b','linewi',3);
plot(datenum(0,1:12,15),tavg(3,:),'color','g','linewi',3);
plot(datenum(0,1:12,15),tavg(4,:),'color','m','linewi',3);
 
subplot(2,2,4) %plot average temp for Juan de Fuca (subplot d)
plot(datenum(0,1:12,15),tavg(5,:),'color','r','linewi',3);
plot(datenum(0,1:12,15),tavg(6,:),'color','b','linewi',3);
plot(datenum(0,1:12,15),tavg(7,:),'color','g','linewi',3);
plot(datenum(0,1:12,15),tavg(8,:),'color','m','linewi',3);