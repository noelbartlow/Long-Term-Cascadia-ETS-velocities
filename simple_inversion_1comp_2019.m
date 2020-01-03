
matpath

 clear all




 
 %%load William_Gorda_N
   load homog_GFs_Nov2019.mat

 
  %load combined_ETS
  load SSE_secular_tot_PBO
  load Heter_GFs_May2019
  

ID_G=ID_G_master;
Lats=Lats_master;
Lons=Lons_master-360;
 
% load all_IDs

zone_w=1.82;%1.82; %transition zone width, in decimal degrees.  1.82=200 km wide, 0.91 = 100 km
Gorda_rake=90; %direction of Gorda plate motion
Gorda_cont=0; %continuously varying Gorda plate direction (N-S plate extension)


skinny=0; %1 for just triangles near tremor, 0 for all triangles
Gorda=0; %different Gorda plate motions given by parameters above

save_file='ETSslip_het_sigw_MCPole_comb_Nov19_PBO';
 
lands=shaperead('~/research/gshhg-shp-2.3.6/GSHHS_shp/i/GSHHS_i_L1.shp', 'UseGeoCoords', true);



%do you want to use all the stations?  West of 120 can save some time
ind_GPS=find(Lons<-120);


ind_Kern=[ind_GPS*3-2; ind_GPS*3-1; ind_GPS*3];
ind_Kern=ind_Kern(:);

%Are you using Charles' GFs?  If so compare ID_G variables, reorder

ID_G2=ID_G(ind_GPS, :);

ind_Charles=[];
ind_Kern=[];
tot=0;
ind_GPS2=[];
indii=[];
for mo=1:length(ID_G2)
    c=GetIndex(ID_G_Charles, ID_G2(mo, :));
    c=c(1);
    if(c>0) %station found in Charles GFs

        ind_Kern=[ind_Kern; c*3-2; c*3-1; c*3];
        disp(c)
        tot=tot+1;
        
        ind_GPS2=[ind_GPS2; ind_GPS(mo)];
        indii=[indii; mo];
        
    end
    
end

ind_GPS=ind_GPS2;
Nsites=size(G1_py, 1)/3;

G1=G1_py;
G2=G2_py;


%% build data vector and Kernel



d=[];
sig=[];
for i = 1:length(ind_GPS)
d=[d; Ev_SSE(ind_GPS(i)); Nv_SSE(ind_GPS(i)); Uv_SSE(ind_GPS(i))];
%d=[d; eV2(i)-eV1(i); nV2(i)-nV1(i); uV2(i)-uV1(i)];
sig=[sig; sig_Ev_SSE(ind_GPS(i)); sig_Nv_SSE(ind_GPS(i)); sig_Uv_SSE(ind_GPS(i))];
end

sig(sig<3e-4)=3e-4;


Pole=[-31.96, 68.30];     %Euler pole for JdF - NA motion from MORVEL 56 the 2010 paper table 3 (more precise value in MORVEL calculator)
%Pole is Lat, Lon
%possible pole with forearc from McCaffrey et al 2007 section 5.7 249.9E, 12.3N, ? =
%0.55° Myr?1 

%that pole is forearc - JdF so flip it:
%-12.3, 69.9


Pole2=[-12.3, 69.9];
Pole=Pole2;
    


    Dist_pole=zeros(length(el), 1);
    Az=zeros(length(el), 1);
    trilon=zeros(length(el), 1);
    trilat=zeros(length(el), 1);
    tridepth=zeros(length(el), 1);
    rake=zeros(length(el), 1);
    Dist_GN=zeros(length(el), 1);
    Az_GN=zeros(length(el), 1);
    Dist_GS=zeros(length(el), 1);
    Az_GS=zeros(length(el), 1);
    for i=1:length(el)
        trilon(i)=mean(nd_ll(el(i, :), 1));
        trilat(i)=mean(nd_ll(el(i, :), 2));
        tridepth(i)=mean(nd_ll(el(i, :), 3));
        
        [Dist_pole(i), Az(i)]=distance(trilat(i), trilon(i)+360, Pole(1), Pole(2));
        
        %degrees CW from north
        rake(i)=90-(360-Az(i));
       
        [A, B]=distance(trilat(i), trilon(i), Gorda_lats+zone_w/2, Gorda_lons);
        Dist_GN(i)=min(A);
        Az_GN(i)=B(A==min(A));
        
        [A, B]=distance(trilat(i), trilon(i), Gorda_S_lats, Gorda_S_lons);
        Dist_GS(i)=min(A);
        Az_GS(i)=B(A==min(A));
    end
    

    
      %solve for strike and dip of each element
    for i = 1:size(el,1)
    norm_vec = getnormal(nd(el(i,1),:),nd(el(i,2),:),nd(el(i,3),:));
    [dip(i),strike(i)] = normal2dipstrike(norm_vec);
    end

  
        %this block of code changes the rake in the south to match the data
    
            dip_dir=strike-90;
            
            
if(Gorda == 1)
    ind_transition = find(Dist_GN<zone_w & (Az_GN<90 | Az_GN>270)); 
    ind_Gorda = find(Dist_GN>zone_w & (Az_GN<90 | Az_GN>270));
else
    ind_transition = []; 
    ind_Gorda = [];
end

clear zeta zeta_frac
if(Gorda_cont==1)
    offsetN=0;
    offsetS=0;
    ind_Gorda = find(Dist_GN>offsetN & (Az_GN<90 | Az_GN>270));
    Gorda_rake=115; %angle of southern edge CW from N on map
    Gorda_x=sind(Gorda_rake);
    Gorda_y=cosd(Gorda_rake);
    JdF_x=sind(rake);
    JdF_y=cosd(rake);

    %add vectors
    zeta_frac=(Dist_GN-offsetN)./(Dist_GN-offsetN+Dist_GS-offsetS);
    zeta_frac(zeta_frac>1)=1;
    x_tot=JdF_x+Gorda_x*zeta_frac;
    y_tot=JdF_y+Gorda_y*zeta_frac;
    %normalize and calculate rake
    x_tot=x_tot./sqrt(x_tot.^2+y_tot.^2);
    y_tot=y_tot./sqrt(x_tot.^2+y_tot.^2);
    rake(ind_Gorda)=90-atan2d(y_tot(ind_Gorda), x_tot(ind_Gorda));
    rake=rake';
else
    zeta=Gorda_rake-rake(ind_transition);
    zeta_frac=(Dist_GN)/zone_w;
    rake(ind_transition)=rake(ind_transition)+zeta_frac(ind_transition).*zeta;
    rake(ind_Gorda)=Gorda_rake;
    rake=rake';
end

%this is CW from north, just like the rake
strike=strike-180;
alpha = strike-rake; 
ss_fact = -cosd(alpha);
ds_fact = -sind(alpha);

   Kern_GPS = G1.*repmat(ss_fact, Nsites*3, 1) + G2.*repmat(ds_fact, Nsites*3, 1);

%   Kern_GPS = G_par_corr;
%fix directions
Kern=Kern_GPS(ind_Kern, :);                                                                                                                                                                                                                                                                                                           



 
%solve for slip direction
%Kern=[-G1(ind_Kern, :), G2(ind_Kern, :)];

%weight d and Kern by inverse variances
  sig_w=sig/mean(sig);
   
   %normalize weights so we don't have to worry about smoothing changes
   
   sig_w=sig_w./mean(sig_w);
   
   
% 
% %weighting of only horizontal vs. vertical
%   sig_w=ones(size(sig))/3;
%   sig_w(3:3:end)=1; %factor of 3 difference between horiz and vert
   
%mixed weighting   
%   sig_w2=sig.^2/mean(sig.^2);
%    sig_w2(sig_w2<.05)=.05;
%    sig_w2(sig_w2>5)=5;
%    sig_w=ones(size(sig))/9;
%    sig_w(3:3:end)=1; %factor of 9 (=3^2) difference between horiz and vert
% 
% sig_w=sqrt(sig_w.^2+sig_w2.^2);
   
%weight all components equally
 %  sig_w=ones(size(sig));




%triangle areas
p12 = nd(el(:, 2), :) - nd(el(:, 1), :);
p13 = nd(el(:, 3), :) - nd(el(:, 1), :);
a = cross(p12',p13');
a = sqrt(sum(a.^2))/2;

  
smooth_mat=mean(abs(Kern_GPS(ind_Kern, :)))/mean(mean(abs(Kern_GPS(ind_Kern, :))));
smooth_mat=smooth_mat'./a'*mean(a);
smooth_mat=repmat(smooth_mat, 1, length(smooth_mat));
%smooth_mat=smoothing*smooth_mat.^-1;

d=d./sig_w;
Kern_w=repmat(sig_w, 1, size(Kern, 2)); %weighting by data uncertainties
%Kern_w=ones(size(Kern)); %no weighting
Kern=Kern./Kern_w;

%restrict to triangles in ETS zone
if(skinny==1)
    ind_ETStri=find(a<100 & trilat'<51.2); %just triangles near tremor;
else
    ind_ETStri=find(a>0); %all triangles
   % ind_ETStri=find(tridepth>-65); %traingles above 65 km depth
end
Kern=Kern(:, ind_ETStri);
Lap=Lap(ind_ETStri, ind_ETStri);
smooth_mat=smooth_mat(ind_ETStri, ind_ETStri);
el=el(ind_ETStri, :);
Kern_w=Kern_w(:, ind_ETStri);


%% Do the inversion

%single slip direction, unweighted smoothing
%slip_rate=lsqnonneg([Kern; smoothing*Lap], [d; zeros(length(Lap), 1)]);
%single slip direction, GF weighted smoothing

%slip_rate=lsqnonneg([Kern; smooth_mat.*Lap], [d; zeros(length(Lap), 1)]);

smoothing_tot=[1, 2, 4, 8, 16, 32, 64, 128, 256];

   smoothing=16;
    smooth_mat2=smoothing*smooth_mat.^-1;
    %single slip direction, GF weighted smoothing
    slip_rate=lsqnonneg([Kern; smooth_mat2.*Lap], [d; zeros(length(Lap), 1)]);
  


%L-curve analysis
% for n=1:length(smoothing_tot)
%     
%     disp('start  ')
%     disp(n)
%     smoothing=smoothing_tot(n);
%     smooth_mat2=smoothing*smooth_mat.^-1;
%     %single slip direction, GF weighted smoothing
%     slip_rate(:, n)=lsqnonneg([Kern; smooth_mat2.*Lap], [d; zeros(length(Lap), 1)]);
%     %slip_rate=lsqlin([Kern; smooth_mat.*Lap], [d; zeros(length(Lap), 1)]);
%     
%     Lap_tot(n)=norm(smooth_mat.^-1.*Lap*slip_rate(:, n));
%     norm_tot(n)=norm(d-Kern*slip_rate(:, n));
%     
%     
%     disp(n)
%     disp('  finished')
% end
% slip_rate2=slip_rate;


%variable slip direction, GF weighted smoothing
%Lap2=[smooth_mat.*Lap, 0.1*smooth_mat.*Lap; 0.1*smooth_mat.*Lap, smooth_mat.*Lap];
%slip_rate=lsqnonneg([Kern; Lap2], [d; zeros(length(Lap2), 1)]);

%slip_rate=ones(size(slip_rate));
RMS_Error=rms(d-Kern*slip_rate); %weighted RMS error
%undo data and Kern weighting
d_plot=d.*sig_w;
Kern=Kern.*Kern_w;
%% make plots

%for L-curve
%slip_rate=slip_rate2(:, 4);

%two component
%total_slip=sqrt(slip_rate(1:length(el)).^2+slip_rate(length(el)+1:end).^2);
%figure; trisurf(el,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),total_slip, 'EdgeAlpha', 0.1);

% two component
comps=Kern*slip_rate;
Ecomp=comps(1:3:end);
Ncomp=comps(2:3:end);
Ucomp=comps(3:3:end);




load pos_map
%one component
figure; trisurf(el,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),slip_rate, 'EdgeAlpha', 0.1);
figure; trisurf(el,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),slip_rate, 'EdgeAlpha', 0.1);


hold on
plot(lands(4).Lon, lands(4).Lat, 'k', 'Linewidth', 2)
sc=70;
quiver([Lons(ind_GPS), -128], [Lats(ind_GPS), 43.5], sc*[d_plot(1:3:end)', 0.01], sc*[d_plot(2:3:end)', 0], 0, 'k', 'linewidth', 2)

 quiver([Lons(ind_GPS), -128], [Lats(ind_GPS), 44], ...
     sc*[Ecomp', 0.01], ...
     sc*[Ncomp', 0], 0, 'r', 'linewidth', 2);
 
 view(2)
colormap(pos_map)
caxis([0 .05])
colorbar
daspect([1/cosd(abs(origin(1))),1,110]);
set(gca, 'fontsize', 16);
plot( -97.653, 16.646, 's', 'markerfacecolor', 'green')
%text(lons, lats, sta_names, 'fontsize', 16)
xlim([-130 -120])
ylim([40 50])
%xlim([-125 -121])
%ylim([39.5 42.5])
title('Total Slip')
 
 figure; trisurf(el,nd_ll(:,1),nd_ll(:,2),nd_ll(:,3),slip_rate, 'EdgeAlpha', 0.1);
 hold on
plot(lands(4).Lon, lands(4).Lat, 'k', 'Linewidth', 2)
sc=50;
quiver([Lons(ind_GPS), -128], [Lats(ind_GPS), 43.5], 0*[d_plot(1:3:end)', 0.01], sc*[d_plot(3:3:end)', 0], 0, 'k', 'linewidth', 2)
% two component
comps=Kern.*Kern_w*slip_rate;
Ecomp=comps(1:3:end);
Ncomp=comps(2:3:end);
Ucomp=comps(3:3:end);
 quiver([Lons(ind_GPS), -128], [Lats(ind_GPS), 44], ...
     0*[Ecomp', 0.01], ...
     sc*[Ucomp', 0], 0, 'r', 'linewidth', 2);

%one component
% quiver([Lons(ind_GPS)'; -128], [Lats(ind_GPS)'; 44], ...
%    sc*[Kern(1:3:end, :)*slip_rate; 0.01], ...
%    sc*[Kern(2:3:end, :)*slip_rate; 0], 0, 'r', 'linewidth', 2)

view(2)
colormap(pos_map)
caxis([0 .05])
colorbar
daspect([1/cosd(abs(origin(1))),1,110]);
set(gca, 'fontsize', 16);
plot( -97.653, 16.646, 's', 'markerfacecolor', 'green')
%text(lons, lats, sta_names, 'fontsize', 16)
xlim([-130 -120])
ylim([40 50])
%xlim([-125 -121])
%ylim([39.5 42.5])
title('Total Slip')



% ss=-slip_rate(1:length(el))';
% ds=-slip_rate(length(el)+1:end)';
% figure;
% quiver(trilon, trilat, -ss.*sind(strike)+ds.*sind(strike+90),  -ss.*cosd(strike)+ds.*cosd(strike+90), 1, 'Linewidth', 1, 'MaxHeadSize', 2)
% xlim([-130 -120])
% ylim([40 50])
%xlim([-125 -121])
%ylim([39.5 42.5])

%calculate likelihood for MLE choosing of smoothing
%following least-squares equations without nonneg constraint
%and following Segall and Matthews 1997 for log likelihood equation
%dropping constant term in Segall & Matthews eqn 38




%% 


save(save_file)


