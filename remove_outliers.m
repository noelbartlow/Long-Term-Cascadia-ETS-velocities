%remove extreme outlies

%also detects and removes offsets!

stations2=stations;

%ID_G=char([stations(:).name]');

thresh=4/1000;
wind=20;

wind2=floor(wind/2);

%remove stations without secular velocity, ID_G comes from UNRsecular2019_1.mat
for i=1:length(stations)
    dur(i)=size(stations(i).data, 1);
    ind=GetIndex(ID_G, stations(i).name);
    inde(i)=ind(1);
end

stations=stations(inde>0 & dur>100);

%cycle through GPS stations
for i=1:length(stations)
   % stations(i).year=deciyear(datenum(stations(i).dates2)); %have to
   % remove this because otherwise the QC in the incoming secular code
   % breaks things
    ind=(1:wind)';
    %inspect each data point
    for j=wind+1:length(stations(i).data)-wind2
        %mean positon of last wind non-outlier days
       % meanpos=median(stations(i).data(ind(end-wind+1):ind(end), 1:3));
        
        %mean positon of last wind days and future wind/2-1 days too
        meanpos=median([stations(i).data((j-wind):j-1, 1:3); stations(i).data(j+1:j+wind2-1, 1:3)]);
        
        %position minus running mean
        pos_dist = stations(i).data(j, 1:3)-meanpos;
        
        %correct for higher noise in vertical (downweight by 2)
        %pos_dist(3)=pos_dist(3)/2;
        
        %don't use data if any component is more that thresh out of mean (*3 for vertical) of last wind points
        if(norm(pos_dist(1:2))<thresh && abs(pos_dist(3))<thresh*3.5)
            ind=[ind; j];
      %  else if(stations(i).year(j)-stations(i).year(ind(end))>.15)
       %         ind=[ind;j];
       %     end
                
        end
    end
    
    stations(i).data_nooutliers=stations(i).data(ind, :);
    %stations2(i).time=stations(i).time(ind, :);
    stations(i).year_nooutliers=stations(i).year(ind, :);
    stations(i).dates2_nooutliers=stations(i).dates2(ind, :);
    stations(i).Jdays_nooutliers=stations(i).Jdays(ind);
    
    s1(i)=j;
    s2(i)=length(stations2(i).data);
    
end
%% detect and remove offsets, from antenna changes (PANGA cleaned files) and earthquakes (USGS database over M 5.5)
%read offset times from PANGA files
%some antenna offsets are missing from PANGA files.  Therefore I added
%the UNR step database.  Still fitting the steps myself.

%Get UNR steps database
%websave('UNR_steps.txt', 'http://geodesy.unr.edu/NGLStationPages/steps.txt');
UNRsteps=importdata('/Users/bartlowno/research/UNR_data/steps.txt');
%extract station names
UNRsteps_names=char();
for k=1:length(UNRsteps)
    UNRsteps_names(k, :)=UNRsteps{k}(1:4);
    UNRsteps_dates(k, :)=UNRsteps{k}(7:13);
    UNRsteps_deciyear(k, :)=deciyear(datenum(UNRsteps_dates(k, :), 'yymmmdd'));
end

%cycle through offsets and correct
for i=1:length(stations)
    ind_UNR=find(strcmp(stations(i).name, string(UNRsteps_names)));
    for j=1:length(ind_UNR)
        ind_offset=find(stations(i).year_nooutliers>UNRsteps_deciyear(ind_UNR(j)));
        data=stations(i).data_nooutliers;

        if(~isempty(ind_offset) && ind_offset(1)>wind+1 && length(ind_offset)>wind)
            offset_mags=median(data(ind_offset(1:wind), 1:3))-median(data(ind_offset(1)-wind-1:ind_offset(1)-1, 1:3));
            stations(i).data_nooutliers(ind_offset, 1:3)=stations(i).data_nooutliers(ind_offset, 1:3)-offset_mags;
      %  else %not enough data to fix the step, so remove data before or after step as necessary.  Program later.
            
        end
    end
end

%antenna offsets from PANGA, if available
no_panga_files=[];
for i=1:length(stations)
    

    Panga_filename=strcat('../panga_cleaned/', ID_G(i, :), '.lon');
    fid=fopen(Panga_filename);
    ind=[];
    if(fid>0)
        fclose(fid);
        S = textread(Panga_filename,'%s', 500);
        TF=strcmp('OFFSET', S);
        ind=find(TF)+1;
    else
        no_panga_files=[no_panga_files; ID_G(i, :)];
    end
     
    
    %cycle over offsets, correct them
    for j=1:length(ind)
        offset_time=str2double(S{ind(j)});
        
        ind_offset=find(stations(i).year_nooutliers>offset_time);
        data=stations(i).data_nooutliers;
        
        if(~isempty(ind_offset) && ind_offset(1)>wind)
            offset_mags=median(data(ind_offset(1:wind), 1:3))-median(data(ind_offset(1)-wind:ind_offset(1)-1, 1:3));
            stations(i).data_nooutliers(ind_offset, 1:3)=stations(i).data_nooutliers(ind_offset, 1:3)-offset_mags;
        end
    end
    

    

end




%Earthquake offsets from USGS
load nearbyquakes
%cycle over earthquakes
%correct offsets at stations within r km of quake epicenters
r=10.^nearbyquakes.mag/1.5e4;
for k=1:length(nearbyquakes.lat)
    
    quake_year=deciyear(datenum(nearbyquakes.Date(k)));

    %cycle over stations, calculate distance to quake
    for i=1:length(stations)
        E=referenceEllipsoid('wgs84', 'km');
        [ARCLEN, AZ] = distance([stations(i).Lat(1), stations(i).Lon(1)], ...
            [nearbyquakes.lat(k), nearbyquakes.lon(k)], E);
        
        %determine if station is close enough to quake, correct offset
        if(ARCLEN<r(k))
            disp(quake_year)
            disp(stations(i).name)
            ind_offset=find(stations(i).year_nooutliers>quake_year-0.5/365);
            data=stations(i).data_nooutliers;
        
            if(~isempty(ind_offset) && length(ind_offset)>wind && ind_offset(1)>wind+1)
                offset_mags=median(data(ind_offset(2:wind), 1:3))-median(data(ind_offset(1)-wind-1:ind_offset(1)-1, 1:3));
                stations(i).data_nooutliers(ind_offset, 1:3)=stations(i).data_nooutliers(ind_offset, 1:3)-offset_mags;
            else if(~isempty(ind_offset) && ind_offset(1)>wind+21)
                    offset_mags=median(data(ind_offset(2:end), 1:3))-median(data(ind_offset(1)-wind-1:ind_offset(1)-1, 1:3));
                    stations(i).data_nooutliers(ind_offset, 1:3)=stations(i).data_nooutliers(ind_offset, 1:3)-offset_mags;
          
                 end
            end

        end

    end
end







        
        
%%

c=1;
%c=84;


%c=58;
c=59;
c=43;

c=GetIndex(ID_G, 'CBRV');
c=GetIndex(char([stations(:).name]'), 'P327');

%c=1;

figure
inipot=stations(c).data(1, 1:3);
long=length(stations(c).data);
long2=length(stations(c).data_nooutliers);



 plot(stations(c).year, detrend(stations(c).data(:, 1:3))*1000 ...
     -repmat(inipot+[0 -20 -50], long, 1), '.')
 hold on
  plot(stations(c).year_nooutliers, detrend(stations(c).data_nooutliers(:, 1:3))*1000 ...
       -repmat(inipot+[10 -10 -40], long2, 1), 'k.')
   
   title(stations(c).name)
