%new secular velocity code for Cascadia

%written by Noel Bartlow, March 6, 2018
%re-written to solve in one step following Paul Segall suggestion, 5/1/2018
clear all



load Tremor_Mar2019_custom
load UNR_stations_new


% %set how far from station tremor detection should be to mask data, in km
 tremor_rad=250;
 
  %length to declare SSE
SSE_L=6;
% 
% %minimum number of tremor detections within radius
T_min=50;



save_file='secular_UNR_Nov19_8.mat';
clear durations

for k=1:length(stations) %inital length QC
    durations(k)=length(stations(k).year);
end

stations=stations(durations>4*365);

good_stations=[];
for k=1:length(stations)
    largest_gap(k)=max(stations(k).year(2:end)-stations(k).year(1:end-1)); 
    if(largest_gap(k) < 0.5)    %greater than 4 years of total data with no gaps over 6 months
        good_stations=[good_stations; k];
    else %if enough data but gaps are the problem check timeseries after last big gap, see if that qualifies
        disp('checking')
        disp(k)
        gap_ind=find((stations(k).year(2:end)-stations(k).year(1:end-1))>= 0.5);
        gap_ind=max(gap_ind);
        stations(k).data=stations(k).data(gap_ind+1:end, :);
        stations(k).year=stations(k).year(gap_ind+1:end);
       
        if(length(stations(k).year)>4*365)
            good_stations=[good_stations; k];
       
        end
   end
end

stations=stations(good_stations);




%load tremor information
%load Tremor_Nov2017
%load Tremor_Jun2018


%load data (through June 2017)

%UNR
%ID_G=char([stations(:).name]');
%PBO
ID_G=char({stations(:).name});

num_sta=length(stations);

%initialize parameters
%initial values
E0=zeros(num_sta, 1);
N0=zeros(num_sta, 1);
U0=zeros(num_sta, 1);
%velocities
Ev=zeros(num_sta, 1);
Nv=zeros(num_sta, 1);
Uv=zeros(num_sta, 1);
%annual sine
Ea1=zeros(num_sta, 1);
Na1=zeros(num_sta, 1);
Ua1=zeros(num_sta, 1);
%annual cosine
Ea2=zeros(num_sta, 1);
Na2=zeros(num_sta, 1);
Ua2=zeros(num_sta, 1);
%semiannual sine
Es1=zeros(num_sta, 1);
Ns1=zeros(num_sta, 1);
Us1=zeros(num_sta, 1);
%semiannual cosine
Es2=zeros(num_sta, 1);
Ns2=zeros(num_sta, 1);
Us2=zeros(num_sta, 1);


%% first loop to identiy data masks

%remove outliers
remove_outliers

for i=1:num_sta
    
    %index that gives tremor from tremor database which is within R_T km of
    %station
    ind_tremor=find(distance(Tremor_lats, Tremor_lons, repmat(stations(i).Lat(1), length(Tremor_lats), 1),...
        repmat(stations(i).Lon(1), length(Tremor_lats), 1))*110<tremor_rad);
    
    %insert minimum tremor detections here
    %group tremor times by day  
    
    for coz=1:length(datetime(stations(i).dates2_nooutliers))
        mo=datetime(stations(i).dates2_nooutliers(coz, :));
        
        Tremor_coz=find(Tremor_times(ind_tremor).Day==mo.Day & Tremor_times(ind_tremor).Month==mo.Month & Tremor_times(ind_tremor).Year==mo.Year);
        
        if length(Tremor_coz)<T_min
            ind_tremor(Tremor_coz)=NaN;
        end
        
        ind_tremor=ind_tremor(~isnan(ind_tremor));
        
    end
    

    
    ind_mask=[];
    for j=1:length(ind_tremor)
        ind_mask=[ind_mask; find(abs(datetime(stations(i).dates2_nooutliers)-Tremor_times(ind_tremor(j)))<hours(12))];
    end
    
    ind_mask=unique(ind_mask);
    
    stations(i).ind_mask=ind_mask;
    stations(i).Masked_data=stations(i).data_nooutliers;
    stations(i).Masked_data(ind_mask, 1)=NaN;
    stations(i).Mask_times=stations(i).year_nooutliers(ind_mask);
    stations(i).Masked_year=stations(i).year_nooutliers(~isnan(stations(i).Masked_data(:, 1)));
    stations(i).Masked_data=stations(i).Masked_data(~isnan(stations(i).Masked_data(:, 1)), :);
    
    if(mod(i, 10)==0)
        disp('masking data, station number ')
        disp(i)
    end
    
end


save masked_stations_UNR_new1
%% 

%remove questionable data after May 27, 2018 when UNR has reference frame
%change which somehow was not in steps database but now is argh
for i=1:length(stations)
    ind1=find(stations(i).year_nooutliers<=2018.4026);
    stations(i).year_nooutliers=stations(i).year_nooutliers(ind1);
    stations(i).data_nooutliers=stations(i).data_nooutliers(ind1, :);
    
    ind2=find(stations(i).Masked_year<=2018.4026);
    stations(i).Masked_year=stations(i).Masked_year(ind2);
    stations(i).Masked_data=stations(i).Masked_data(ind2, :);
end

%% solve for everything at once
num_sta=length(stations);
clear Uv Ev Nv

%number of bootstrap tries to determine population of inter-ETS velocities
N=20;

for i=1:num_sta
    %detect data gaps. Gap_inds is index of just before data gap.  Only
    %detect gaps larger than SSE_L days
    
    if(mod(i,10)==0)
        disp('fitting station  ')
        disp(i)
    end
        
            %choose begin time based on tremor catalog
    if (stations(i).Lat(1)<49.5 && stations(i).Lat(1)>40.5)
        ind_times=find(stations(i).Masked_year>2009.6 & stations(i).Masked_year<2019.211); %to match end of tremor catalog
    else
        ind_times=find(stations(i).Masked_year>2012.15 & stations(i).Masked_year<2019.211); %to match end of tremor catalog
    end
    
    %if not enough (less than 1 year) after tremor catalog start, then at least fit seasonals to
    %existing data
    bad_vel_inds=[];
    
    if(length(ind_times)<365)
        ind_times=1:length(stations(i).Masked_year);
        bad_vel_inds=[bad_vel_inds; i];
    end
    
    t=stations(i).Masked_year(ind_times);
    
    gap_inds=find(abs(t(2:end, :)-t(1:end-1, :))>SSE_L/365);
    for j=length(gap_inds):-1:1 %check if gap is due to tremor, not outliers or missing data
        check=find((stations(i).Mask_times-t(gap_inds(j)))<2/365);
        if(isempty(check))
            gap_inds=gap_inds([1:j-1, j+1:end]);
        end
    end
    
    
    
    %best soluition, use all data
    %build G matrix

    G=ones(length(t), length(gap_inds)+6);
    %velocity term
    G(:, 2)=t-2000;
    %seasonal terms
    G(:, 3)=sin(2*pi*t);
    G(:, 4)=cos(2*pi*t);
    G(:, 5)=sin(4*pi*t);
    G(:, 6)=cos(4*pi*t);
    %ETS offset terms
    for j=1:length(gap_inds)
        G(:, j+6)=zeros(length(t), 1);
        G(gap_inds(j)+1:end, j+6)=ones(length(t)-gap_inds(j), 1);
    end
    
    dE=stations(i).Masked_data(ind_times, 1);
    dN=stations(i).Masked_data(ind_times, 2);
    dU=stations(i).Masked_data(ind_times, 3);
    
    mE=G\dE;
    mN=G\dN;
    mU=G\dU;
    
    stations(i).secular_mE=mE;
    stations(i).secular_mN=mN;
    stations(i).secular_mU=mU;
    stations(i).secular_fitE=G*mE;
    stations(i).secular_fitN=G*mN;
    stations(i).secular_fitU=G*mU;
    stations(i).secular_G=G;
    
     %calculate residuals (d-Gm)
    %Large residuals indicate unmodeled deformation sources
    resE=stations(i).Masked_data(ind_times, 1) - G*mE; 
    resN=stations(i).Masked_data(ind_times, 2) - G*mN; 
    resU=stations(i).Masked_data(ind_times, 3) - G*mU; 

    stations(i).RMSE=rms(resE);
    stations(i).RMSN=rms(resN);
    stations(i).RMSU=rms(resU);
    
    %turns out it's better to use post-ETS fit residuals
    
    Ev(i)=mE(2);
    Nv(i)=mN(2);
    Uv(i)=mU(2);
    Ea1(i)=mE(3);
    Na1(i)=mN(3);
    Ua1(i)=mU(3);
    Ea2(i)=mE(4);
    Na2(i)=mN(4);
    Ua2(i)=mU(4);
    Es1(i)=mE(5);
    Ns1(i)=mN(5);
    Us1(i)=mU(5);
    Es2(i)=mE(6);
    Ns2(i)=mN(6);
    Us2(i)=mU(6);
    
    %BOOTSTRAP LOOP
    %remove 20% of data from the ends to produce population of inter-ETS
    %velocities
    
    boot_frac=0.2;
    
    stations(i).secular_Ev_bootpop=[];
    stations(i).secular_Nv_bootpop=[];
    stations(i).secular_Uv_bootpop=[];
    
    for boot=1:N
        %build G matrix
        G=ones(length(t), length(gap_inds)+6);
        %velocity term
        G(:, 2)=t-2000;
        %seasonal terms
        G(:, 3)=sin(2*pi*t);
        G(:, 4)=cos(2*pi*t);
        G(:, 5)=sin(4*pi*t);
        G(:, 6)=cos(4*pi*t);
        %ETS offset terms
        for j=1:length(gap_inds)
            G(:, j+6)=zeros(length(t), 1);
            G(gap_inds(j)+1:end, j+6)=ones(length(t)-gap_inds(j), 1);
        end
        
        dE=stations(i).Masked_data(ind_times, 1);
        dN=stations(i).Masked_data(ind_times, 2);
        dU=stations(i).Masked_data(ind_times, 3);
    
        %remove 20% randomly from end on both data and G
        begin_frac=rand(1)*boot_frac;
        end_frac=boot_frac-begin_frac;
        begin_ind=round(length(dE)*begin_frac)+1;
        end_ind=round(length(dE)-length(dE)*end_frac);
        
        dE=dE(begin_ind:end_ind);
        dN=dN(begin_ind:end_ind);
        dU=dU(begin_ind:end_ind);
        G=G(begin_ind:end_ind, :);
        
        ind_cols=1; %always keep the first col of all 1s
        %remove rows of all 1s or all 0s from G
        for mo=1:size(G, 2)
            if(sum(G(:, mo))~=0 && sum(G(:, mo))~=size(G, 1) )
                ind_cols=[ind_cols; mo];
            end
        end
        
        G=G(:, ind_cols);
    
        mEboot=G\dE;
        mNboot=G\dN;
        mUboot=G\dU;

        stations(i).secular_Ev_bootpop=[stations(i).secular_Ev_bootpop; mEboot(2)];
        stations(i).secular_Nv_bootpop=[stations(i).secular_Nv_bootpop; mNboot(2)];
        stations(i).secular_Uv_bootpop=[stations(i).secular_Uv_bootpop; mUboot(2)];

    end
    

end

%% If station doesn't have enough (or any) data after start of tremor catalog, 
%then use inter-ETS velocity from nearby stations instead
clear Lats Lons
for i=1:num_sta
    Lats(i)=stations(i).Lat(1);
    Lons(i)=stations(i).Lon(1);

end

for i=bad_vel_inds

    dists=distance(Lats(i), Lons(i), Lats([1:i-1, i+1:end]), Lons([1:i-1, i+1:end]));
    A=find(dists<.2);
    if(isempty(A))
        A=find(dists==min(dists));
    end
    
        A(A>i)=A(A>i)+1;
    
    
    Ev(i)=mean(Ev(A));
    Nv(i)=mean(Nv(A));
    Uv(i)=mean(Uv(A));
    sig_Ev(i)=mean(sig_Ev(A));
    sig_Nv(i)=mean(sig_Nv(A));
    sig_Uv(i)=mean(sig_Uv(A));
    
    stations(i).secular_mE(2)=Ev(i);
    stations(i).secular_mN(2)=Nv(i);
    stations(i).secular_mN(2)=Uv(i);
    
    stations(i).secular_fitE=stations(i).secular_G*stations(i).secular_mE;
    stations(i).secular_fitN=stations(i).secular_G*stations(i).secular_mN;
    stations(i).secular_fitU=stations(i).secular_G*stations(i).secular_mU;
    
   
    
end




%% 
% 
figure; worldmap([38 52], [-130 -114])
title('Inter-SSE velocities');
load coast

plotm(lat, long, 'k')

hold on
sc=50;
scatterm(Lats, Lons, sc, Uv', 'filled')



quiverm([Lats, 40], [Lons, 360-127], sc*[Nv, 0], sc*[Ev, 0.01], 'k', 0)

%plot error ellipses on map
% ind_Nv=find(sig_Nv>sig_Ev); %semimajor axis N
% ind_Ev=find(sig_Ev>=sig_Nv); %semimajor axis E
% 
% [elat, elon]=ellipse1(Lats(ind_Nv)'+sc*Nv(ind_Nv), Lons(ind_Nv)'-360+sc*Ev(ind_Nv), [sc*sig_Nv(ind_Nv)', sqrt(1-sig_Ev(ind_Nv).^2./sig_Nv(ind_Nv).^2)'], 90);
% plotm(elat, elon, 'r')
% [elat, elon]=ellipse1(Lats(ind_Ev)'+sc*Nv(ind_Ev), Lons(ind_Ev)'-360+sc*Ev(ind_Ev), [sc*sig_Ev(ind_Ev)', sqrt(1-sig_Nv(ind_Ev).^2./sig_Ev(ind_Ev).^2)']);
% plotm(elat, elon, 'r')




caxis([-0.005 0.005])
%% Solve for and plot total SSE velocities (displacements/time, time varies by station)
clear Ev_SSE Nv_SSE Uv_SSE Ev_SSE Nv_SSE Uv_SSE sig_Ev_SSE sig_Nv_SSE sig_Uv_SSE


bad_data=[];

for i=1:length(stations)
     %build G matrix for secular & seasonal removal
    t=stations(i).year_nooutliers;
    G=ones(length(stations(i).data_nooutliers), 6);
    %velocity term
    G(:, 2)=t-2000;
    %seasonal terms
    G(:, 3)=sin(2*pi*t);
    G(:, 4)=cos(2*pi*t);
    G(:, 5)=sin(4*pi*t);
    G(:, 6)=cos(4*pi*t);
    
    %build G matrix for SSE vel fit
    t=stations(i).year_nooutliers;
    G2=ones(length(stations(i).data_nooutliers), 2);
    %velocity term
    G2(:, 2)=t-2000;
    
    
    mE2=G2\(stations(i).data_nooutliers(:, 1) - G*stations(i).secular_mE(1:6));
    mN2=G2\(stations(i).data_nooutliers(:, 2) - G*stations(i).secular_mN(1:6));
    mU2=G2\(stations(i).data_nooutliers(:, 3) - G*stations(i).secular_mU(1:6));
    
    Ev_SSE(i)=mE2(2);
    Nv_SSE(i)=mN2(2);
    Uv_SSE(i)=mU2(2);
    
%     if(length(stations(i).year_nooutliers)<365*2)
%         Ev_SSE(i)=NaN;
%         Nv_SSE(i)=NaN;
%         Uv_SSE(i)=NaN; 
%         bad_data=[bad_data; i];
%         
%     end
    
         %calculate residuals (d-Gm)
    %Large residuals indicate unmodeled deformation sources
    resE=(stations(i).data_nooutliers(:, 1) - G*stations(i).secular_mE(1:6)) - G2*mE2; 
    resN=(stations(i).data_nooutliers(:, 2) - G*stations(i).secular_mN(1:6)) - G2*mN2; 
    resU=(stations(i).data_nooutliers(:, 3) - G*stations(i).secular_mU(1:6)) - G2*mU2; 

    stations(i).ETS_RMSE=rms(resE);
    stations(i).ETS_RMSN=rms(resN);
    stations(i).ETS_RMSU=rms(resU);
    
    if(sqrt(rms(resE).^2+rms(resN).^2)>3.5e-3 | rms(resU)>10.5e-3)
        bad_data=[bad_data; i];
    end
    
end

        
        
%% Remove bad stations



for c=1:length(bad_data)
        k=bad_data(end-c+1);
        stations=stations([1:k-1, k+1:end]);
        Ev=Ev([1:k-1, k+1:end]);
        Ev_SSE=Ev_SSE([1:k-1, k+1:end]);
        Ea1=Ea1([1:k-1, k+1:end]);
        Ea2=Ea2([1:k-1, k+1:end]);
        Es1=Es1([1:k-1, k+1:end]);
        Es2=Es2([1:k-1, k+1:end]);
        Nv=Nv([1:k-1, k+1:end]);
        Nv_SSE=Nv_SSE([1:k-1, k+1:end]);
        Na1=Na1([1:k-1, k+1:end]);
        Na2=Na2([1:k-1, k+1:end]);
        Ns1=Ns1([1:k-1, k+1:end]);
        Ns2=Ns2([1:k-1, k+1:end]);        
        Uv=Uv([1:k-1, k+1:end]);
        Uv_SSE=Uv_SSE([1:k-1, k+1:end]);
        Ua1=Ua1([1:k-1, k+1:end]);
        Ua2=Ua2([1:k-1, k+1:end]);
        Us1=Us1([1:k-1, k+1:end]);
        Us2=Us2([1:k-1, k+1:end]);
end

num_sta=length(stations);
Lats1=Lats;
Lons1=Lons;
ID_G1=ID_G;
clear Lats Lons ID_G
for i=1:num_sta
    Lats(i)=stations(i).Lat(1);
    Lons(i)=stations(i).Lon(1);

end
%ID_G=char([stations(:).name]');
ID_G=char({stations(:).name});
        %% PLOT SSE VELOCITIES

figure; worldmap([38 51], [-129 -118])
title('SSE velocities')
load coast
plotm(lat, long, 'k')
for i=1:num_sta
    Lats(i)=stations(i).Lat(1);
    Lons(i)=stations(i).Lon(1);

end
load blue_red_map
set(gcf, 'colormap', blue_red_map)
scatterm(Lats, Lons, 20, Uv_SSE, 'filled')
caxis([-0.012 0.012])
sc=70;
quiverm([Lats, 40], [Lons, 360-127], sc*[Nv_SSE, 0], sc*[Ev_SSE, 0.01], 'k', 0)

%plot error ellipses on map
% ind_Nv=find(sig_Nv_SSE>sig_Ev_SSE & sig_Ev_SSE>0); %semimajor axis N
% ind_Ev=find(sig_Ev_SSE>=sig_Nv_SSE & sig_Nv_SSE>0); %semimajor axis E
% 
% [elat, elon]=(Lats(ind_Nv)'+sc*Nv_SSE(ind_Nv)', Lons(ind_Nv)'-360+sc*Ev_SSE(ind_Nv)', [sc*sig_Nv_SSE(ind_Nv), sqrt(1-sig_Ev_SSE(ind_Nv).^2./sig_Nv_SSE(ind_Nv).^2)]);
% plotm(elat, elon, 'k')
% [elat, elon]=ellipse1(Lats(ind_Ev)'+sc*Nv_SSE(ind_Ev)', Lons(ind_Ev)'-360+sc*Ev_SSE(ind_Ev)', [sc*sig_Ev_SSE(ind_Ev), sqrt(1-sig_Nv_SSE(ind_Ev).^2./sig_Ev_SSE(ind_Ev).^2)], 90);
% plotm(elat, elon, 'k')









%% Make some plots

%ID_G=char([stations(:).name]');
%ID_G=char({stations(:).name});

c=166;
c=GetIndex(ID_G, 'P324'); %784, 370, 380, 663, QUIN, 154, 179

     %build G matrix
    t=stations(c).year_nooutliers;
    G=ones(length(stations(c).data_nooutliers), 6);
    %velocity term
    G(:, 2)=t-2000;
    %seasonal terms
    G(:, 3)=sin(2*pi*t);
    G(:, 4)=cos(2*pi*t);
    G(:, 5)=sin(4*pi*t);
    G(:, 6)=cos(4*pi*t);
    
            %choose begin time based on tremor catalog
    if (stations(c).Lat(1)<49.5 && stations(c).Lat(1)>40.5)
        ind_times=find(stations(c).Masked_year>2009.6 & stations(c).Masked_year<2019.211); %to match end of tremor catalog
    else
        ind_times=find(stations(c).Masked_year>2012.15 & stations(c).Masked_year<2019.211); %to match end of tremor catalog
    end
    
    
    
    t=stations(c).Masked_year(ind_times);
    


figure;

subplot(3, 1, 1)
set(gca, 'fontsize', 16)
hold on
title([ID_G(c, :), ' East'])
%plot(stations(c).year, stations(c).data(:, 1), '.');
plot(stations(c).Masked_year, stations(c).Masked_data(:, 1)*1000, '.');
plot(t, stations(c).secular_fitE*1000, 'k', 'Linewidth', 2)
plot(stations(c).year_nooutliers, G*stations(c).secular_mE(1:6)*1000, '-.r', 'Linewidth', 2)
%plot(stations(c).Masked_year, stations(c).Corrected_data(:, 1), 'r.');
%plot(stations(c).year, feval(model, E0(c), Ev(c), Ea1(c), Ea2(c), Es1(c), Es2(c), stations(c).year), 'k', 'Linewidth', 2)
%xlim([2009 2017.8])

%mark mendocino EQ times
% plot([2005+166/365, 2005+166/365], [-1 1], 'k')
% plot([2010+9/365, 2010+9/365], [-1 1], 'k')
% plot([2014+69/365, 2014+69/365], [-1 1], 'k')
% plot([2016+343/366, 2016+343/366], [-1 1], 'k')

subplot(3, 1, 2)
set(gca, 'fontsize', 16)
hold on
title([ID_G(c, :), ' North'])
%plot(stations(c).year, stations(c).data(:, 2), '.');
plot(stations(c).Masked_year, stations(c).Masked_data(:, 2)*1000, '.');
plot(t, stations(c).secular_fitN*1000, 'k', 'Linewidth', 2)
plot(stations(c).year_nooutliers, G*stations(c).secular_mN(1:6)*1000, '-.r', 'Linewidth', 2)
%plot(stations(c).Masked_year, stations(c).Corrected_data(:, 2), 'r.');
%plot(stations(c).year, feval(model, N0(c), Nv(c), Na1(c), Na2(c), Ns1(c), Ns2(c), stations(c).year), 'k', 'Linewidth', 2)
%xlim([2009 2017.8])

subplot(3, 1, 3)
set(gca, 'fontsize', 16)
hold on
title([ID_G(c, :), ' Up'])
%plot(stations(c).year, stations(c).data(:, 3), '.');
plot(stations(c).Masked_year, stations(c).Masked_data(:, 3)*1000, '.');
plot(t, stations(c).secular_fitU*1000, 'k', 'Linewidth', 2)
plot(stations(c).year_nooutliers, G*stations(c).secular_mU(1:6)*1000, '-.r', 'Linewidth', 2)
%plot(stations(c).Masked_year, stations(c).Corrected_data(:, 3), 'r.');
%plot(stations(c).year, feval(model, U0(c), Uv(c), Ua1(c), Ua2(c), Us1(c), Us2(c), stations(c).year), 'k', 'Linewidth', 2)
%xlim([2009 2017.8])

    
figure;



subplot(3, 1, 1)
set(gca, 'fontsize', 16)
hold on
title([ID_G(c, :), ' East minus inter-ETS and seasonals'])
yE=stations(c).data_nooutliers(:, 1) - G*stations(c).secular_mE(1:6);
plot(stations(c).year_nooutliers, yE*1000, '.');
%plot(stations(c).year_nooutliers, smooth(yE, 50, 'rloess'), '-g', 'linewidth', 2);

%mark mendocino EQ times
% plot([2005+166/365, 2005+166/365], [min(yE)-0.003 max(yE)+0.003], 'k')
%  plot([2010+9/365, 2010+9/365], [min(yE)-0.003 max(yE)+0.003]*1000, 'k')
%  plot([2014+69/365, 2014+69/365], [min(yE)-0.003 max(yE)+0.003]*1000, 'k')
%  plot([2016+343/366, 2016+343/366], [min(yE)-0.003 max(yE)+0.003]*1000, 'k')


subplot(3, 1, 2)
set(gca, 'fontsize', 16)
hold on
title([ID_G(c, :), ' North minus inter-ETS and seasonals'])
yN=stations(c).data_nooutliers(:, 2) - G*stations(c).secular_mN(1:6);
plot(stations(c).year_nooutliers, yN*1000, '.');

%mark mendocino EQ times
% plot([2005+166/365, 2005+166/365], [min(yN)-0.003 max(yN)+0.003], 'k')
%  plot([2010+9/365, 2010+9/365], [min(yN)-0.003 max(yN)+0.003]*1000, 'k')
%  plot([2014+69/365, 2014+69/365], [min(yN)-0.003 max(yN)+0.003]*1000, 'k')
%  plot([2016+343/366, 2016+343/366], [min(yN)-0.003 max(yN)+0.003]*1000, 'k')

subplot(3, 1, 3)
set(gca, 'fontsize', 16)
hold on
title([ID_G(c, :), ' Up minus inter-ETS and seasonals'])
yU=stations(c).data_nooutliers(:, 3) - G*stations(c).secular_mU(1:6);
plot(stations(c).year_nooutliers, yU*1000, '.');

figure;
subplot(3, 1, 1)
set(gca, 'fontsize', 16)
hold on
title([ID_G(c, :), ' East minus full model (no ETS left)'])
yE=stations(c).Masked_data(ind_times, 1) - stations(c).secular_fitE;
plot(t, yE*1000, '.');

%mark mendocino EQ times
 %plot([2005+166/365, 2005+166/365], [min(yE)-0.003 max(yE)+0.003], 'k')
 %plot([2010+9/365, 2010+9/365], [min(yE)-0.003 max(yE)+0.003]*1000, 'k')
%  plot([2014+69/365, 2014+69/365], [min(yE)-0.003 max(yE)+0.003]*1000, 'k')
%  plot([2016+343/366, 2016+343/366], [min(yE)-0.003 max(yE)+0.003]*1000, 'k')


subplot(3, 1, 2)
set(gca, 'fontsize', 16)
hold on
title([ID_G(c, :), ' North minus full model (no ETS left)'])
% yN=stations(c).data_nooutliers(ind_times, 2) - stations(c).secular_G*stations(c).secular_mN;
% plot(t, yN*1000, 'o');
yN=stations(c).Masked_data(ind_times, 2) - stations(c).secular_fitN;
plot(t, yN*1000, '.');
%plot(stations(c).year_nooutliers, smooth(yN, 50, 'rloess'), '-g', 'linewidth', 2);

%mark mendocino EQ times
% plot([2005+166/365, 2005+166/365], [min(yN)-0.003 max(yN)+0.003], 'k')
 %plot([2010+9/365, 2010+9/365], [min(yN)-0.003 max(yN)+0.003]*1000, 'k')
%  plot([2014+69/365, 2014+69/365], [min(yN)-0.003 max(yN)+0.003]*1000, 'k')
%  plot([2016+343/366, 2016+343/366], [min(yN)-0.003 max(yN)+0.003]*1000, 'k')

subplot(3, 1, 3)
set(gca, 'fontsize', 16)
hold on
title([ID_G(c, :), ' Up minus full model (no ETS left)'])
yU=stations(c).Masked_data(ind_times, 3) - stations(c).secular_fitU;
plot(t, yU*1000, '.');
%plot(stations(c).year_nooutliers, smooth(yU, 100, 'rloess'), '-g', 'linewidth', 2);



%% %% 

save(save_file)
