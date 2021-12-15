clc; clear all; close all;

rng('default'); 

l = qd_layout;  %network layout of the simulation run    
l.name = 'Simple Scenario';
l.simpar.center_frequency = 5.9e9; %carrier frequency is 5.9 GHz    
l.tx_name = {'LTE BS'}; %the transmitter is an LTE BS - omni directional antenna
l.rx_name = {'CV2X UE'}; %the receiver is a CV2X UE BS - omni directional antenna
l.tx_array.center_frequency=l.simpar.center_frequency; %antennas operate at the carrier frequency
l.simpar.use_absolute_delays = 1; %include delay of the LOS path                 
l.simpar.show_progress_bars = 0; %don't print simulation stats                    
l.tx_position = [ 5,5,10 ]'; %BS fixed position                     
l.rx_track = qd_track( 'linear' , 40, pi ); %linear 40m track going west 
l.rx_track.name = 'CV2X UE'; %name of the trajectory of the car
l.rx_track.initial_position = [25 ; 25 ; 1.5 ]; %MT initial postion      
l.rx_track.interpolate_positions(0.2); % One channel sample every 0.02 cm
l.rx_track.scenario = {'3GPP_38.901_UMi_NLOS'}; %3GPP densely populated urban area       
l.rx_track.no_segments = l.rx_track.no_snapshots;  % Use spatial consisteny for mobility
%the track is divided in segments. For each segment we take a snapshot.
%This is what the 3GPP model does (p.198 quadriga manual)
%with the parameter 0.2 in the interpolate_positions function and the 40
%length track I can create 9 snapshots, i.e. 9 segments. The track is
%divided in 9 segments. For the scenario used please refer to the p.198 of
%the quadriga manual. For the description of the scenario please refer to
%p.101. I did not use berlin, because the berlin scenario (if you open
%the config file) is for frequencies around 2-3Ghz.

%TAP VALUES-> given by the equation (141) in quadriga manual p.142, they
%are stored in pow. For the 9 segments of mobility and the 20
%clusters(taps)
%size(pow)
%
%ans =
%
%     9    20


%build the coefficients
%for lines 40-44 please refer to p.198 of the documentation
b = l.init_builder;                                     
b.scenpar.PerClusterDS = 0;                     % Disable per-cluster delay spread
b.scenpar.KF_mu = -3;                           % Set los power to 33 % of the total power
b.gen_parameters;                               % Generate small-scale-fading parameters

c = get_channels( b );                          % Generate channel coefficients
c = merge( c, [], 0 );                          % Combine output channels
c.individual_delays = 0;                        % Remove per-antenna delays

dl = c.delay.';                                 % path delays

%these are the channel gains
pow = squeeze( abs(c.coeff).^2 )';              % path powers (quadriga equation)

%for plotting--------------------------------------------------------------
[len,dist] = get_length( l.rx_track ); %Store the length and distances from start point

%also some plotting helpers------------------------------------------------
snapshots=c.no_snap; %save the number of snapshots/segments
paths=c.no_path; %save the number of paths/taps
dist = transpose(dist);
dist_rows = snapshots;
dist_cols = paths;
dist = zeros(dist_rows,dist_cols)+dist;

figure(1);       

h=zeros(paths,snapshots);                                 
for s=1:snapshots
    for p=1:paths
        h(p,s)=stem3(dist(s,:),dl(s,:).*1e6,10.*log10(pow(s,:)));
    end
 if s==1,hold on,end
end

set(gca,'Ydir','reverse', 'Zdir', 'reverse');

%plot the CIR

%--colors in the RGB spectrum for plotting (https://www.tydac.ch/color/)---
colorspec=cell(1, paths);
for tt = 1:paths
  colorspec{tt}=[rand(1,1) rand(1,1) rand(1,1)]; %random RGB value
end
%--------------------------------------------------------------------------
%preallocation-------------------------------------------------------------
Legend=cell(1,paths); 
L=zeros(1,paths);
%--------------------------------------------------------------------------
for i=1:paths
    stem3(dist(:,i),dl(:,i).*1e6,10.*log10(pow(:,i)),'Color',colorspec{i},'Color',colorspec{i},'LineWidth',1.5);
    if i==1
        Legend{i}=strcat('LOS');
        L(i) = plot(nan, nan, 'Color',colorspec{i});
    else
        Legend{i}=strcat('NLOS', num2str(i-1));
        L(i) = plot(nan, nan, 'Color',colorspec{i});
    end
    hold on;
end

title('Path Power with mobility');
xlabel('Distance from start point [m]');
ylabel('Delay [\mus]');
zlabel('Path Gain [dB]'); %path loss
legend(L,Legend);

figure(2);

%P_t=zeros(snapshots,paths)+30; %dbm
pow2=10.*log10(pow)+30;

h2=zeros(paths,snapshots);        
for s=1:snapshots
    for p=1:paths
        h2(p,s)=stem3(dist(s,:),dl(s,:).*1e6,pow2(s,:));
    end
 if s==1,hold on,end
end

set(gca,'Ydir','reverse', 'Zdir', 'reverse');

for i=1:paths
    stem3(dist(:,i),dl(:,i).*1e6,pow2(:,i),'Color',colorspec{i},'Color',colorspec{i},'LineWidth',1.5);
    if i==1
        Legend{i}=strcat('LOS');
        L(i) = plot(nan, nan, 'Color',colorspec{i});
    else
        Legend{i}=strcat('NLOS', num2str(i-1));
        L(i) = plot(nan, nan, 'Color',colorspec{i});
    end
    hold on;
end

title('Power Delay Profile (PDP) with mobility');
xlabel('Distance from start point [m]');
ylabel('Delay [\mus]');
zlabel('Received Power [dBm]');
legend(L, Legend);


