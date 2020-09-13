close all
clear all
format long


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
% Read in recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

filenames = input('Enter recording filename (.wav format) ','s');
[ya,fs] = audioread(filenames);

% Manually enter the filename here
%[ya,fs]=audioread('Experiment.wav');

dt = 1.d0/fs; % Time between recordings
n = length(ya); % Number of recordings

totaltime = n*dt;
fprintf('Press Ctrl C together anytime to end analysis \n'); 
fprintf('Length of recording in seconds %8.7f \n',totaltime);
fprintf('Total number of data points %d \n',n);
fprintf('Sampling rate in Hertz %8.7f \n',1.0/dt);
fprintf('Time step in milliseconds %8.7f \n',1000*dt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
% End: Read in recording
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yamax = max(abs(ya));

fprintf('Maximum absolute height of spikes %d \n',yamax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
% Set parameters for analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
num_selections = 2; % Number of features to use for spike sorting
%totalsignal = 0; % Set totalsignal = 1 to plot the entire recording

% Parameter which manages addition of points to cluster, Suggested range: [.25,1.5]
fraction_radius_group_orig = .6;
fraction_radius_group = fraction_radius_group_orig;

% Parameter which manages distance between clusters, Suggested range: [.25,1]
fraction_radius_all_scores_orig = 0.5;
fraction_radius_all_scores = fraction_radius_all_scores_orig;

% Parameter which regulates the minimum distance the points of one cluster
% can be to another cluster
minimum_distance_cluster_parameter_orig = .125;
minimum_distance_cluster_parameter_orig = .25;
minimum_distance_cluster_parameter = minimum_distance_cluster_parameter_orig;

ijump = 10;

% Spike features to select from to do analysis
fprintf('\n')
fprintf('*************   Selection parameters *********************')
fprintf('\n')
fprintf('1: Height 2: Negative height \n');
fprintf('3: Half positive width 4: Full positive width \n');
fprintf('5: Half negative width 6: Full negative width  \n');
fprintf('7: Positive area 8: Negative area 9: Total area\n');
fprintf('10: Distance between positive and negative peaks \n');
fprintf('********************************************************** \n')
fprintf('\n')

stra = cellstr(char('Height','Negative height','Width (ms)',...
                    'Full positive width','Half negative width','Full negative width',...
                    'Positive area','Negative area','Total area',...
                    'Distance between positive and negative peaks',...
                    'PCA Score 1','PCA Score 2'));
pca_temp = -1;
nit = 0;
while (pca_temp ~= 0 && pca_temp ~= 1)
   nit = nit + 1;
   if (nit > 1)
      fprintf('Please enter 0 or 1 \n') 
   end
   pca_temp = input('Enter 0 for manual selection or 1 to use Principal Component Analysis '); 
end

pca = pca_temp;

if (pca == 0)
    for ijk = 1:2
        if (ijk == 1)
           iii = 0;
           while (iii > 10 || iii < 1)
              iii = input('Enter number for selection criteria 1: ');
              fprintf('Value needs to be between 1 and 10 \n');
           end
           selection(ijk) = iii;
        else
           iii = 0;
           while (iii > 10 || iii < 1)  
              iii = input('Enter number for selection criteria 2: '); 
              fprintf('Value needs to be between 1 and 10 \n');
           end
           selection(ijk) = iii;
        end
    end 
    %selection(1) = 3;
    %selection(2) = 1;
end


% Set pca = 1 to use Principal Component Analysis to perform clustering
% pca = 0;
% pca = 1;

if (pca == 1)
   selection(1) = 11;
   selection(2) = 12;
end

splt = cellstr(char('r-','g-','m-','b-','y-','c-'));
splto = cellstr(char('ro','go','mo','bo','yo','co'));
spltb = cellstr(char('r','g','m','b','y','c'));

createplot = 1; % Set createplot to create high resolution plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
% End: Set parameters for analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%   Plot recording, Select start and end times for analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%fprintf('Enter 1 to plot total signal, 0 to continue without plotting total signal \n');
%fprintf('If a 1 is entered, the program will end due to time required to plot entire signal. \n');
%fprintf('Program can then be rerun and 0 can be entered \n');
totalsignal = input('Enter 1 to plot total signal, 0 to continue '); 
if (totalsignal == 1) 
   ifirst = 1;
   ilast = n;
   if (n > 1000000)
      fprintf('Since the recording is very large %d \n',n)
      fprintf('Plot every %d point of recording \n',ijump);
      ijump = 10;
   end
   xpkt = [ifirst:ijump:ilast];
   xpkt = xpkt*dt;
   ypkt = ya(ifirst:ijump:ilast);
   axmin = double(ifirst)*dt;
   axmax = double(ilast)*dt;
   ypkmin = min(ypkt);
   ypkmax = max(ypkt);
   subplot(2,1,1)
   plot(xpkt,ypkt)
   axis([axmin axmax ypkmin ypkmax])
   fprintf('Plotting total signal \n');
   xlabel('Time (seconds)','FontSize',16)
   ylabel('Amplitude','FontSize',16)
   set(gca,'linewidth',2)
   set(gca,'FontSize',12)
   %fprintf('Plotting entire signal for 5 seconds \n')
   hold on
   %pause(5)
   %close all
   %fprintf('Press enter to continue \n');
   %pause
   %fprintf('Ending program \n')
   %return
end

starttime_temp = -1;
while (starttime_temp < 0 || starttime_temp >= totaltime)
    fprintf('Start time needs to be greater than 0 and less than total time %8.7f \n',totaltime);
    starttime_temp = input('Enter start time (in seconds) to begin analysis ');
end
starttime = starttime_temp;
endtime_temp = starttime;
while(endtime_temp <= starttime || endtime_temp > totaltime)
    fprintf('End time needs to be greater than start time %d and less than total time %8.7f \n',starttime,totaltime);
    endtime_temp = input('Enter end time (in seconds) at which to end analysis ');
end
endtime = endtime_temp;
threshold_temp = -1.;
while (threshold_temp < 0 || threshold_temp > yamax)
    fprintf('Spike threshold needs to be greater than 0 and less than the maximum spike height %8.7f \n',yamax);
    threshold_temp = input('Enter minimum threshold for spikes to be considered for analysis '); 
end
threshold = threshold_temp;
    
fprintf('Threshold in spike sorting algorithm %4.4f \n',threshold);
  
mstart = floor(max(floor(starttime/dt)+1,2)); % Beginning integer recording value corresponding to start time
mfinish = floor(max(floor(endtime/dt))); % End integer recording value corresponding to end time
        
ioverlay = 1; % Set overlay = 1 to create plot of overlayed spikes

fprintf('Extracting recording from starttime to endtime \n')
minypk = 1.d+20;
maxypk = -1.d+20;
xpk(:,1) = [mstart*dt:dt:mfinish*dt];
ypk = ya(mstart:mfinish,1);
minypk = min(ypk);
maxypk = max(ypk);

nelem = mfinish - mstart + 1;
if (nelem > 1000000)
  ijump_selection = 10;
else
  ijump_selection = 1;
end

xpk_plot = xpk(1:ijump_selection:nelem);
ypk_plot = ypk(1:ijump_selection:nelem);

if (totalsignal == 1) 
   plot(xpk_plot,ypk_plot,'c-')
   axis([axmin axmax ypkmin ypkmax])
   hold off
end

xp_threshold(1) = xpk(1);
yp_threshold(1) = threshold;
xp_threshold(2) = xpk(mfinish-mstart+1,1);
yp_threshold(2) = threshold;
if (totalsignal == 1)
   subplot(2,1,2)
end
plot(xpk_plot,ypk_plot,'c-')
if (totalsignal == 1)
   axis([min(xpk_plot) max(xpk_plot) ypkmin ypkmax])
end
hold on
plot(xp_threshold,yp_threshold,'r-','Linewidth',2)
%if (totalsignal == 1)
%   legend('Total signal','Selected signal','Threshold')
%else
%   legend('Selected signal','Threshold')
%end
legend('Selected signal            ','Threshold','Orientation','horizontal','Location','South')
%legend('','Threshold')
legend boxoff
xlabel('Time (seconds)','FontSize',16)
ylabel('Amplitude','FontSize',16)
set(gca,'linewidth',2)
set(gca,'FontSize',12)
%fprintf('Plotting recording from selected start to end times for 5 seconds\n');
if (createplot == 1)
   highres('totalsignal')
end


fprintf('Are you satisfied with your choice of the spike threshold and start and end times for analysis? \n')
ithreshold = input('Enter 0: No (program will end and can be rerun) OR 1: Yes ');
if (ithreshold == 0)
   return
end

hold off
close all
%fprintf('Press enter to continue \n');
%pause
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%   End: Plot recording, Select start and end times for analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%           Find locations of peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
fprintf('Finding location of peaks \n') 
istart = 0;
ilocation = 0;
ypkmax = 0.d0;
ijk = 0;
peak = zeros(mfinish-mstart+1,1); % Stores the recordings where the peaks occur

for i = mstart:mfinish
    if (abs(ya(i,1)) < threshold)
       istart = 1;
    end
    if (abs(ya(i,1)) > threshold && istart == 1)
       if (abs(ya(i,1)) > ypkmax)
          ypkmax = abs(ya(i,1));
          ilocation = i;
       end
    end
    if (ilocation > 0)
       if (abs(ya(i-1,1)) >= .5d0*threshold && abs(ya(i,1)) <= .5d0*threshold && istart == 1)
          ijk = ijk + 1;
          
          if (mod(ijk,400) == 0)
             fprintf('Number of peaks (positive or negative) found above threshold %d \n',ijk)
          end
          peak_pre(ijk) = ilocation;
          height_pre(ijk) = ya(ilocation,1);

          facya = .5*(ya(ilocation+1,1) - ya(ilocation-1,1))/(ya(ilocation+1,1)-2.0*ya(ilocation,1)+ya(ilocation-1,1));
          peak_int(ijk) = dt*(double(ilocation) - facya);
                    
          ypkmax = 0.d0;
          ilocation = 0;
          istart = 0;
       end
    end
    
end

npeaks = ijk;

ifilter = 0;
for ijk = 1:npeaks
    if (height_pre(ijk) > 0.0)
       ifilter = ifilter + 1;
       peak(ifilter) = peak_pre(ijk);
       height(ifilter) = height_pre(ijk);
    else
       iloc = peak_pre(ijk);
       if (ijk + 1 <= npeaks)
          iloc_ahead = peak_pre(ijk+1);
       else
          iloc_ahead = -100000;
       end
       if (ijk - 1 >= 1)
          iloc_before = peak_pre(ijk-1);
       else
          iloc_before = -100000;
       end
       time_ahead = double(iloc_ahead - iloc)*dt/1000.; % time to next spike in ms
       time_before = double(iloc - iloc_before)*dt/1000.; % time to previous spike in ms
       if (time_ahead <= 2.0 || time_before <= 2.0)
          too_close = 1; % Spike is within 2 ms of another spike
       else
          too_close = 0;
       end
       max_allowed = floor(1.5/(dt*1000)); % Positive peak needs to be within 1.5 ms
       if (too_close == 0)
          % Search for positive peak
          ii = iloc;
          ypkmax = 0.0;
          ilocation = 0;
          while (ii <= max_allowed)
              if (ya(ii,1) > ypkmax && ya(ii,1) >= .5*threshold)
                 ypkmax = ya(ii,1);
                 ilocation = ii;
              end
              ii = ii + 1;
          end
          ilocation
          
          if (location > 0) 
             ifilter = ifilter + 1;
             peak(ifilter) = ilocation;
             height(ifilter) = ypkmax;
          end
       end
    end
end

npeaks = ifilter;
       
fprintf('Number of peaks above threshold %d \n \n',npeaks);

if (npeaks < 10)
   fprintf('Number of peaks may be too few to do a spike sorting analysis \n')
   fprintf('Try decreasing threshold or increasing time range \n');
   fprintf('Ending progam \n');
   return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%           End: Find locations of peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%  Determining other features of action potential
%  For example: negpeak is location of negative peak, 
%               peak is location of positive peak
%  beginpeak firstneg negpeak secondneg firstzero beghalfheight peak endhalfheight endpeak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
firstzero = zeros(npeaks,1);
beginpeak = zeros(npeaks,1);
endpeak = zeros(npeaks,1);
           
for kk = 1:npeaks
    ijk = 0;
    firstzerol = 0;
    secondpeak = 0;
    heightnegative(kk) = 1.d+10;      
          
    for j = peak(kk):-1:peak(kk)-200
        if (ya(j,1) >= 0.d0 && ya(j-1,1) <= 0.d0 && firstzerol == 0 && secondpeak == 0)
           firstzerol = 1;
           firstzero(kk) = j; % firstzero stores the location of first zero to left of peak
           firstzero_int(kk) = dt*( double(j-1) + (0.d0 - ya(j-1,1))/(ya(j,1) - ya(j-1,1)) );
        end
        if (ya(j,1) <= 0.d0 && ya(j-1,1) >= 0.d0 && firstzerol == 1 && secondpeak == 0)
           beginpeak(kk) = j; % beginpeak stores the second zero to left of peak
           beginpeak_int(kk) = dt*( double(j-1) + (0.d0 - ya(j-1,1))/(ya(j,1) - ya(j-1,1)) );
           secondpeak = 1;
        end
    end
    beginpeak(kk) = max(beginpeak(kk),peak(kk)-200);

    heightnegative(kk) = 1.d+10;
    areaneg(kk) = 0.d0;
    for j = beginpeak(kk):firstzero(kk)
        if (ya(j,1) <= heightnegative(kk)) 
           heightnegative(kk) = ya(j,1); % stores location of negative peak
           negpeak(kk) = j;                    
        end
        areaneg(kk) = areaneg(kk) + abs(ya(j,1));
    end
                
    jlo = negpeak(kk);
    facya = .5*(ya(jlo+1,1) - ya(jlo-1,1))/(ya(jlo+1,1)-2.0*ya(jlo,1)+ya(jlo-1,1));
    negpeak_int(kk) = dt*(double(jlo) - facya);
                                           
    firstzeror = 0;
    for j = peak(kk)+5:peak(kk)+100
        %if (ya(j-1,1) >= 0.d0 && ya(j,1) <= 0.d0 && firstzeror == 0)
        if (ya(j,1) > ya(j-1,1) && firstzeror == 0)
           endpeak(kk) = j; % endpeak stores location of zero to right of peak
           facya = .5*(ya(j,1) - ya(j-2,1))/(ya(j,1)-2.0*ya(j-1,1)+ya(j-2,1));
           endpeak_int(kk) = dt*(double(j-1) - facya);
           firstzeror = 1;
        end
    end
    endpeak(kk) = min(endpeak(kk),peak(kk)+100);

    areapos(kk) = 0.d0;
    for j = firstzero(kk):endpeak(kk)
        areapos(kk) = areapos(kk) + abs(ya(j,1));
    end
          
    halfheightnegative = .5d0*heightnegative(kk);
              
    for j = beginpeak(kk):firstzero(kk)
        if (ya(j-1,1) >= halfheightnegative && ya(j,1) <= halfheightnegative)
           firstneg(kk) = j; % firstneg and secondneg store locations of negative half height
           firstneg_int(kk) = dt*( double(j-1) + (halfheightnegative - ya(j-1,1))/(ya(j,1) - ya(j-1,1)) );
        end 
        if (ya(j-1,1) <= halfheightnegative && ya(j,1) >= halfheightnegative)
           secondneg(kk) = j;
           secondneg_int(kk) = dt*( double(j-1) + (halfheightnegative - ya(j-1,1))/(ya(j,1) - ya(j-1,1)) );
        end 
    end
          
    halfheight = .5d0*height(kk);
    for j = firstzero(kk):peak(kk)
        if (ya(j-1,1) <= halfheight && ya(j,1) >= halfheight)
           beghalfheight(kk) = j; % beghalfheight and endhalfheight store locations of positive half height
           beghalfheight_int(kk) = dt*( double(j-1) + (halfheight - ya(j-1,1))/(ya(j,1) - ya(j-1,1)) );
        end
    end
          
    for j = peak(kk)+1:endpeak(kk)
        if (ya(j-1,1) >= halfheight && ya(j,1) <= halfheight)
           endhalfheight(kk) = j;  
           endhalfheight_int(kk) = dt*( double(j-1) + (halfheight - ya(j-1,1))/(ya(j,1) - ya(j-1,1)) ); 
        end
    end
          
end
            
% beginpeak firstneg negpeak secondneg firstzero beghalfheight peak endhalfheight endpeak  

maxdistb = -1.d0;
for kk = 1:npeaks 
    distancebeginpeak = peak(kk) - beginpeak(kk);
    maxdistb = max(maxdistb,distancebeginpeak);
end
ap_before = .5;
ap_after = .2;

maxdistb = floor(ap_before/(1000*dt));
   
maxdiste = -1.d0;
for kk = 1:npeaks 
    distanceendpeak = endpeak(kk) - peak(kk);
    maxdiste = max(maxdiste,distanceendpeak);
end
maxdiste = floor(ap_after/(1000*dt));
              
for kk = 1:npeaks
    halfposheightdur(kk) = double(endhalfheight(kk)-beghalfheight(kk))*dt;
    halfposheightdur(kk) = endhalfheight_int(kk) - beghalfheight_int(kk);

    fullposheightdur(kk) = double(endpeak(kk)-firstzero(kk))*dt;
    fullposheightdur(kk) = endpeak_int(kk) - firstzero_int(kk);

    halfnegheightdur(kk) = double(secondneg(kk)-firstneg(kk))*dt;
    halfnegheightdur(kk) = secondneg_int(kk) - firstneg_int(kk);

    fullnegheightdur(kk) = double(firstzero(kk)-beginpeak(kk))*dt;
    fullnegheightdur(kk) = firstzero_int(kk) - beginpeak_int(kk);
                
    distancepeaks(kk) = double(peak(kk) - negpeak(kk))*dt;
    distancepeaks(kk) = peak_int(kk) - negpeak_int(kk);

    areapos(kk) = areapos(kk)*dt;
    areaneg(kk) = areaneg(kk)*dt;
    areatot(kk) = areapos(kk) + areaneg(kk);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%  End: Determining other features of action potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%  Normalize features of action potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
heightmean = mean(height);
heightsd = std(height);
   
heightnegativemean = mean(heightnegative);
heightnegativesd = std(heightnegative);
   
areaposmean = mean(areapos);
areapossd = std(areapos);
   
areanegmean = mean(areaneg);
areanegsd = std(areaneg);
   
areatotmean = mean(areatot);
areatotsd = std(areatot);
   
halfposheightdurmean = mean(halfposheightdur);
halfposheightdursd = std(halfposheightdur);
   
fullposheightdurmean = mean(fullposheightdur);
fullposheightdursd = std(fullposheightdur);
   
halfnegheightdurmean = mean(halfnegheightdur);
halfnegheightdursd = std(halfnegheightdur);
    
fullnegheightdurmean = mean(fullnegheightdur);
fullnegheightdursd = std(fullnegheightdur);
   
distancepeaksmean = mean(distancepeaks);
distancepeakssd = std(distancepeaks);

% Normalize features of action potential
for kk = 1:npeaks
    heightz(kk) = (height(kk) - heightmean)/heightsd;
    heightnegativez(kk) = (heightnegative(kk) - heightnegativemean)/heightnegativesd;
    areaposz(kk) = (areapos(kk) - areaposmean)/areapossd;
    areanegz(kk) = (areaneg(kk) - areanegmean)/areanegsd;
    areatotz(kk) = (areatot(kk) - areatotmean)/areatotsd;
    halfposheightdurz(kk) = (halfposheightdur(kk) - halfposheightdurmean)/halfposheightdursd;
    fullposheightdurz(kk) = (fullposheightdur(kk) - fullposheightdurmean)/fullposheightdursd;
    halfnegheightdurz(kk) = (halfnegheightdur(kk) - halfnegheightdurmean)/halfnegheightdursd;
    fullnegheightdurz(kk) = (fullnegheightdur(kk) - fullnegheightdurmean)/fullnegheightdursd;
    distancepeaksz(kk) = (distancepeaks(kk) - distancepeaksmean)/distancepeakssd;
end

stra = cellstr(char('Height','Negative height','Half positive width ',...
                    'Full positive width','Half negative width','Full negative width',...
                    'Positive area','Negative area','Total area',...
                    'Distance between positive and negative peaks',...
                    'PCA Score 1','PCA Score 2'));

coefficient_variation(1) = heightsd/heightmean;
coefficient_variation(2) = -heightnegativesd/heightnegativemean;
coefficient_variation(3) = abs(halfposheightdursd/halfposheightdurmean);
coefficient_variation(4) = fullposheightdursd/fullposheightdurmean;
coefficient_variation(5) = halfnegheightdursd/halfnegheightdurmean;
coefficient_variation(6) = fullnegheightdursd/fullnegheightdurmean;
coefficient_variation(7) = areapossd/areaposmean;
coefficient_variation(8) = areanegsd/areanegmean;
coefficient_variation(9) = areatotsd/areatotmean;
coefficient_variation(10) = abs(distancepeakssd/distancepeaksmean);

[coef_variation_sort,points_sort_coef] = sort(coefficient_variation,2,'descend'); 

% Print out coefficient of variation from largest to smallest
for i = 1:10
    jj = points_sort_coef(i);
    fprintf('Coefficient of variation of %s is %d \n',stra{jj},coefficient_variation(jj));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%  End: Normalize features of action potential
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
fprintf('\n');   
fprintf('First selection criteria %s \n', stra{selection(1)});
fprintf('Second selection criteria %s \n', stra{selection(2)});
    
xscore = zeros(npeaks,1);
yscore = zeros(npeaks,1);
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%           Manual Selection of Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
if (pca == 0) 
   for ijk = 1:num_selections  
       
       if (selection(ijk) == 1)
          for kk = 1:npeaks
              score(kk,ijk) = heightz(kk);
              scoreo(kk,ijk) = height(kk);
           end 
       end
                    
       if (selection(ijk) == 2)
          for kk = 1:npeaks
              score(kk,ijk) = heightnegativez(kk);
              scoreo(kk,ijk) = heightnegative(kk);
          end 
       end
                    
       if (selection(ijk) == 3)
          for kk = 1:npeaks
              score(kk,ijk) = halfposheightdurz(kk);
              scoreo(kk,ijk) = halfposheightdur(kk);
          end 
       end
          
       if (selection(ijk) == 4)
          for kk = 1:npeaks
              score(kk,ijk) = fullposheightdurz(kk);
              scoreo(kk,ijk) = fullposheightdur(kk);
          end 
       end
          
       if (selection(ijk) == 5)
          for kk = 1:npeaks
              score(kk,ijk) = halfnegheightdurz(kk);
              scoreo(kk,ijk) = halfnegheightdur(kk);
          end 
       end
          
       if (selection(ijk) == 6)
          for kk = 1:npeaks
              score(kk,ijk) = fullnegheightdurz(kk);
              scoreo(kk,ijk) = fullnegheightdur(kk);
          end 
       end
                    
       if (selection(ijk) == 7)
          for kk = 1:npeaks
              score(kk,ijk) = areaposz(kk);
              scoreo(kk,ijk) = areapos(kk);
          end 
       end
          
       if (selection(ijk) == 8)
          for kk = 1:npeaks
              score(kk,ijk) = areanegz(kk);
              scoreo(kk,ijk) = areaneg(kk);
          end 
       end
          
       if (selection(ijk) == 9)
          for kk = 1:npeaks
              score(kk,ijk) = areatotz(kk);
              scoreo(kk,ijk) = areatot(kk);
          end 
       end
                    
       if (selection(ijk) == 10)
          for kk = 1:npeaks
              score(kk,ijk) = distancepeaksz(kk);
              scoreo(kk,ijk) = distancepeaks(kk);
          end 
       end  
                    
   end
          
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%           End: Manual Selection of Criteria
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%           Principal Component Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
range = maxdiste+maxdistb+1;
xpkk = zeros(range,1);
ypkk = zeros(range,1);
xmat = zeros(npeaks,range);
for kk = 1:npeaks
       
    % beginpeak firstneg negpeak secondneg firstzero beghalfheight peak endhalfheight endpeak 
       
    ijk = 0;
    if (peak(kk)-maxdistb < beginpeak(kk)-1)
       for j = peak(kk)-maxdistb:beginpeak(kk)-1
           ijk = ijk + 1;
           xpkk(ijk) = ijk;
           ypkk(ijk) = 0.d0;
       end  
    end
       
    for j = firstneg(kk):peak(kk)
        ijk = ijk + 1;
        xpkk(ijk) = ijk;
        ypkk(ijk) = ya(j,1);
    end
       
    for j = peak(kk)+1:endpeak(kk)
        ijk = ijk + 1;
        xpkk(ijk) = ijk;
        ypkk(ijk) = ya(j,1);
    end
       
    if (peak(kk)+maxdiste > endpeak(kk)+1)
       for j = endpeak(kk)+1:peak(kk)+maxdiste
           ijk = ijk + 1;
           xpkk(ijk) = ijk;
           ypkk(ijk) = 0.d0;
       end  
    end
       
    ijk = 0;
    for j = peak(kk)-maxdistb:peak(kk)+maxdiste
        ijk = ijk + 1;
        xmat(kk,ijk) = ypkk(ijk);
    end
       
end

   
if (pca == 1)
      
    sigma = zeros(range,range);
    U = zeros(range,range);
    V = zeros(range,range);
    S = zeros(range,range);
    Vr = zeros(npeaks,npeaks);

    xmatt = zeros(range,npeaks);
    mean = zeros(npeaks,1);
      
    xmatt = xmat';
    xmatsave = xmatt;
     
    for i = 1:range
        mean(i) = 0.d0;
        for k = 1:npeaks
            mean(i) = mean(i) + xmatt(i,k);
        end
        mean(i) = mean(i)/double(npeaks);
    end
      
    for i = 1:range
        for k = 1:npeaks
            xmatt(i,k) = xmatt(i,k) - mean(i);
        end
    end
      
    sig = xmatt; % sig is a ilength by ilength vector

    % xmat is a npeak (m) BY number_of_data_points_action_potential (n)
    %%sigma = xmat'*xmat;
    %%sigma = sigma/double(npeaks); %Covariance matrix

    %%[U,S,V] = svd(sigma); % U is a n BY n matrix, U' is a n BY n matrix
      
    [U,S,Vr] = svd(sig);
      
    pcascores = zeros(range,npeaks);
                

    pcascores = U'*xmatt; % scores is a n BY m matrix
      
    for ijk = 1:num_selections
        for i = 1:npeaks
            score(i,ijk) = pcascores(ijk,i);
            scoreo(i,ijk) = score(i,ijk);
        end
    end
       
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%           End: Principal Component Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
      
fprintf('\n');

normalize_scores = 1;
% Set normalize_scores = 0 to use unnormalized scores
if (normalize_scores == 1)
   fprintf('Normalizing scores by subtracting mean and dividing by standard deviation \n')
end
for i = 1:npeaks
    xscore(i) = score(i,1);
    yscore(i) = score(i,2);
    xscoreo(i) = scoreo(i,1);
    yscoreo(i) = scoreo(i,2);
    if (normalize_scores == 0)                 
       xscore(i) = scoreo(i,1);
       yscore(i) = scoreo(i,2);
    end
end

                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clustering algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
fprintf('Begin clustering algorithm \n')
%fprintf('If points are being added too quickly, reduce fraction_radius_group \n')
%fprintf('If points are being skipped or added too slowly or have stopped, increase fraction_radius_group \n')
%fprintf('\n')
   
% The radius of each point is defined to be distance to the "number_neighbors" 
% point which are ordered from least to greatest in terms of distance.  
% The smaller the radius, the heigher the density of each point.
% The density of each point is used to select centers of clusters. 

number_neighbors = max(2,floor(.1*npeaks)); 

points_neighbor_sorted = zeros(npeaks,npeaks);
dist_neighbor_sorted = zeros(npeaks,npeaks);
cluster_group = zeros(npeaks,1);
in_cluster = zeros(npeaks,1);
center_group_element = zeros(npeaks,1);
dist_veca = zeros(npeaks,npeaks);
for i = 1:npeaks
    radius_group(i) = 1.e+20;
end
  
fprintf('Finding distances between scores \n')  
for i = 1:npeaks 
    if (mod(i,400) == 0)
      fprintf('Finished calculating distance up to peak %d out of %d peaks \n',i,npeaks);
    end
    for j = 1:npeaks
        dist_vec(j) = 0.d0;
        for ijk = 1:num_selections
            dist_vec(j) = dist_vec(j) + (score(i,ijk)-score(j,ijk))^2;
        end
        dist_vec(j) = sqrt(dist_vec(j));
        % dist_veca(i,j) is the distance from score i to score j
        dist_veca(i,j) = dist_vec(j);
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        
    [dist_sort,points_sort] = sort(dist_vec,2,'ascend');          
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                        

    % Define the distance for point i to the number_neighbors closest point
    % number_neighbors is defined to be the 1/5 of the total number of
    % peaks
    dist_neighbor(i) = dist_sort(number_neighbors); 

    for j = 1:npeaks
        points_neighbor_sorted(i,j) = points_sort(j);
        dist_neighbor_sorted(i,j) = dist_sort(j);
    end
      
end
  
number_groups = 0;
iretry_groups = 1;
fprintf('End finding distances between scores \n')  
   
while (iretry_groups == 1 && number_groups < 4)
   
   itry_again_cluster = 1;
   nit_center = 0;
   
   while (itry_again_cluster == 1)
       
      nit_center = nit_center + 1;
      % Find the point with the highest density
      number_groups = number_groups + 1;
      center_group_element(number_groups) = 0;  

      for i = 1:npeaks    

          % Check to see if the peak i has already been included in a cluster
          % or if the peak i is too close to a cluster
          ifound = 0; 

          if (in_cluster(i) == 0)
             for j = 1:number_groups-1
                 for ii = 1:num_group_elements(j)
                     jj = group_elements(j,ii);
                     dd = dist_veca(i,jj);
                     % radius_all_scores is a measure of the size of the entire set of
                     % scores
                     if ((i == jj) || (dd < fraction_radius_all_scores*radius_all_scores)  ) 
                        ifound = 1;
                     end
                 end
             end
          else
             ifound = 1;
          end
          
          if (ifound == 0) 
             % Find the point i with the smallest radius to the nearest number_neighbors
             % or equivalently the point i with the highest surrounding
             % density of scores
             if (dist_neighbor(i) < radius_group(number_groups))
                center_group_element(number_groups) = i;
                radius_group(number_groups) = dist_neighbor(i);
             end 
          end
       
      end
               
      kc = center_group_element(number_groups);

      if (kc == 0)
         fprintf('\n')
         fprintf('No new cluster center was not found \n')
         fprintf('Try decreasing value of variable fraction_radius_all_scores \n')
         fprintf('Current value of fraction_radius_all_scores %d \n',fraction_radius_all_scores)
         fprintf('Enter 0 to stop finding new clusters \n')
         itry_again_cluster = input('Enter 1 to continue finding cluster by decreasing fraction_radius_all_scores ');
         if (itry_again_cluster == 1)
            fraction_radius_all_scores = input('Enter new value of fraction_radius_all_scores: Suggested range [.25,1] ');
         end
         number_groups = number_groups - 1;
      else
         itry_again_cluster = 0;
      end
      
   end
 
   if (kc > 0)
      number_in_cluster = 1;
      cluster_group = 0;
      cluster_group(1) = kc;
      in_cluster(kc) = 1;
      
      rad_search = fraction_radius_group*radius_group(number_groups);
      iretry = 1;
      nit = 0;
      nit_progress = 0;
      
      fprintf('Adding points to cluster %d \n',number_groups);
      while (iretry == 1) 
         
         nit = nit + 1;
       
         close_points = 0;
         number_in_cluster_new = number_in_cluster;
         for i = 1:number_in_cluster
             ii = cluster_group(i); 

             jj = 0;
             radius_exceeded = 0;
             while (jj < npeaks && radius_exceeded == 0)
                 jj = jj + 1;
                 ic = points_neighbor_sorted(ii,jj);
                 if (dist_neighbor_sorted(ii,jj) < rad_search)      
                    
                    % Determine if candidate point ic is far enough way
                    % from other clusters
                    ifound = 0;
                    if (number_groups > 1)
                       min_distance_neighboring_cluster = minimum_distance_cluster_parameter*radius_all_scores;
                    end
                    if (in_cluster(ic) == 0)
                       close_points = close_points + 1; 
                       for ja = 1:number_groups-1
                           for iia = 1:num_group_elements(ja)
                               jja = group_elements(ja,iia);
                               dd = dist_veca(ic,jja);
                               % radius_all_scores is a measure of the size
                               % of the entire set of scores
                               
                               % fraction_radius_all_score*radius_all_scores
                               % is the minimum distance that must be
                               % maintained between cluster centers
                               
                               % min_distance_neighboring_cluster is the
                               % minimum distance that must be maintainrf
                               % between any point in a cluster and points
                               % in a neighboring cluster
                               
                               if ( dd <  min_distance_neighboring_cluster ) 
                                  ifound = 1;
                               end                  
                           end
                       end
                    else
                        ifound = 1;
                    end

                    if (ifound == 0)
                       number_in_cluster_new = number_in_cluster_new + 1;
                       cluster_group(number_in_cluster_new) = ic;
                       in_cluster(ic) = 1;
                    end
          
                 else                     
                    radius_exceeded = 1;
                 end
             end           
         end
         no_change = 0;
         if (number_in_cluster == number_in_cluster_new  && nit > 1)
             nit_progress = nit_progress + 1;
             fprintf('No additional points were added to cluster \n');
             if (close_points == 0)
                fprintf('Suggest increasing parameter fraction_radius_group \n');
                fprintf('to retry adding additional points to cluster \n');
                fprintf('Current value of fraction_radius_group is %d \n',fraction_radius_group);
                fprintf('Enter 0 to keep value of fraction_radius_group \n')
                itry_again = 0;
                itry_again = input('Enter 1 to change value of fraction_radius_group ');
                if (itry_again == 1)
                   fprintf('Suggested values for fraction_radius_group [.6 1.5] \n');
                   fraction_radius_group = input('Enter value of fraction_radius_group ');
                   rad_search = fraction_radius_group*radius_group(number_groups);
                end
             else
                fprintf('Suggest decreasing minimum distance parameter minimum_distance_cluster_parameter \n')
                fprintf('to be maintained to nearest cluster \n')
                fprintf('Current value of minimum_distance_cluster_parameter is %d \n',minimum_distance_cluster_parameter);
                fprintf('Enter 0 to keep value of minimum_distance_cluster_parameter \n')
                itry_again = 0;
                itry_again = input('Enter 1 to change value of minimum_distance_cluster_parameter ');
                if (itry_again == 1)
                   fprintf('Suggested values for minimum_distance_cluster_parameter [.05 .5] \n');
                   minimum_distance_cluster_parameter = input('Enter value of minimum_distance_cluster_parameter ');
                end
             end
             no_change = 1;
         end
          
         number_in_cluster = number_in_cluster_new;
            
         if (number_groups == 1)
            % radius_all_scores is a measure of the size of the entire set of
            % scores
            radius_all_scores = dist_neighbor_sorted(kc,npeaks);
         end
        
         if (no_change == 0)     
            if (number_groups == 1)
               for i = 1:number_in_cluster
                   j = cluster_group(i);
                   xscoreplota(i) = xscore(j);
                   yscoreplota(i) = yscore(j);
               end
               subplot(1,2,1)
               plot(xscore,yscore,'co',xscoreplota,yscoreplota,char(splto(1)),'Markersize',10)
               xlabel(stra(selection(1)),'Fontsize',14);
               ylabel(stra(selection(2)),'Fontsize',14);
            end
            if (number_groups == 2)
               for i = 1:number_in_cluster
                   j = cluster_group(i);
                   xscoreplotb(i) = xscore(j);
                   yscoreplotb(i) = yscore(j);
               end
               plot(xscore,yscore,'co',xscoreplota,yscoreplota,char(splto(1)),xscoreplotb,yscoreplotb,char(splto(2)),'Markersize',10)
               xlabel(stra(selection(1)),'Fontsize',14);
               ylabel(stra(selection(2)),'Fontsize',14);
            end
            if (number_groups == 3)
               for i = 1:number_in_cluster
                   j = cluster_group(i);
                   xscoreplotc(i) = xscore(j);
                   yscoreplotc(i) = yscore(j);
               end
               plot(xscore,yscore,'co',xscoreplota,yscoreplota,char(splto(1)),xscoreplotb,yscoreplotb,char(splto(2)),xscoreplotc,yscoreplotc,char(splto(3)),'Markersize',10)
               xlabel(stra(selection(1)),'Fontsize',12);
               ylabel(stra(selection(2)),'Fontsize',12);
            end
            if (number_groups == 4)
               for i = 1:number_in_cluster
                   j = cluster_group(i);
                   xscoreplotd(i) = xscore(j);
                   yscoreplotd(i) = yscore(j);
               end
               plot(xscore,yscore,'co',xscoreplota,yscoreplota,char(splto(1)),xscoreplotb,yscoreplotb,char(splto(2)),xscoreplotc,yscoreplotc,char(splto(3)),xscoreplotd,yscoreplotd,char(splto(4)),'Markersize',10)
               xlabel(stra(selection(1)),'Fontsize',12);
               ylabel(stra(selection(2)),'Fontsize',12);
            end
         end            
         fprintf('Number of points in cluster %d is %d \n',number_groups,number_in_cluster)
         iretry = input('Refer to figure.  Add to cluster?  Yes: 1, No: 0 to stop adding points to cluster ');
      end 
      
      if (nit_progress > 0)
          
         if (fraction_radius_group ~= fraction_radius_group_orig)
            fprintf('Return fraction_radius_group to its original value? %d \n',fraction_radius_group_orig);
            fprintf('Enter 0 to keep value of fraction_radius_group \n')
            iorig = input('Enter 1 to return fraction_radius_group to its original value ');
            if (iorig == 1)
               fraction_radius_group = fraction_radius_group_orig;
            end
         end
         
         if (minimum_distance_cluster_parameter ~= minimum_distance_cluster_parameter_orig)
            fprintf('Return minimum_distance_parameter_cluster to its original value? %d ',minimum_distance_cluster_parameter_orig);
            fprintf('Enter 0 to keep value of minimum_distance_parameter_cluster \n')
            iorig = input('Enter 1 to return minimum_distance_parameter_cluster to its original value ');
            if (iorig == 1)
               minimum_distance_cluster_parameter = minimum_distance_cluster_parameter_orig;
            end
         end
      end
      
      
      [cluster_sort,cpoints_sort] = sort(cluster_group,2,'ascend');  
      
      num_group_elements(number_groups) = number_in_cluster;
      for i = 1:num_group_elements(number_groups)
          %group_elements(number_groups,i) = cluster_group(i);
          group_elements(number_groups,i) = cluster_sort(i);
      end
      
      fprintf('\n')
      if (number_groups < 4)
         iretry_groups = input('Enter 1 to construct a new cluster, 0 to stop adding clusters ');
      else
         iretry_groups = 0;
      end
   else
      iretry_groups = 0;
   end
                         
end                 

if (createplot == 1)
   highres('cluster_plot') 
end
  
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End: Clustering algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot similar spikes.  The spikes are overlayed to match the peak to view
% how similar they are.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%hold off
ngroupsplot = min(6,number_groups);
if (ioverlay == 1)
   fprintf('Overlaying spikes \n');    
   mingx = 1.e+20;
   maxgx = -1.e+20;
   for i = 1:ngroupsplot
       ypkk_average = zeros(maxdiste+maxdistb+1,1);
       for k = 1:num_group_elements(i)
           kk = group_elements(i,k);
           ijk = 0;
           xpkk = zeros(maxdiste+maxdistb+1,1);
           ypkk = zeros(maxdiste+maxdistb+1,1);

           for j = peak(kk)-maxdistb:peak(kk)+maxdiste  
               ijk = ijk +1;
               xpkk(ijk) = ijk;
               ypkk(ijk) = ya(j,1);
               ypkk_average(ijk) = ypkk_average(ijk) + ypkk(ijk);
           end
           if (mod(i,7) == 1)
               subplot(1,2,2);
               mingx = min(mingx,1.1*min(ypkk));
               maxgx = max(maxgx,1.1*max(ypkk));
               plot(1000.*xpkk*dt,ypkk,char(splt(1)),'LineWidth',2)
               axis([0 ap_before+ap_after mingx maxgx])
               if (k == num_group_elements(i))
                  ypkk_average = ypkk_average/double(num_group_elements(i));
                  plot(1000.*xpkk*dt,ypkk_average,'k-','LineWidth',2)
               end
               xlabel('Time (ms)','FontSize',14)
               ylabel('Amplitude','FontSize',14)
               set(gca,'linewidth',2)
               set(gca,'FontSize',12)
               hold on
            end
            if (mod(i,7) == 2)
               mingx = min(mingx,1.1*min(ypkk));
               maxgx = max(maxgx,1.1*max(ypkk)); 
               plot(1000.*xpkk*dt,ypkk,char(splt(2)),'LineWidth',2)
               axis([0 ap_before+ap_after mingx maxgx])
               if (k == num_group_elements(i))
                  ypkk_average = ypkk_average/double(num_group_elements(i));
                  plot(1000.*xpkk*dt,ypkk_average,'k-','LineWidth',2)
               end
            end
            if (mod(i,7) == 3)
               mingx = min(mingx,1.1*min(ypkk));
               maxgx = max(maxgx,1.1*max(ypkk)); 
               plot(1000.*xpkk*dt,ypkk,char(splt(3)),'LineWidth',2)
               axis([0 ap_before+ap_after mingx maxgx])
               if (k == num_group_elements(i))
                  ypkk_average = ypkk_average/double(num_group_elements(i));
                  plot(1000.*xpkk*dt,ypkk_average,'k-','LineWidth',2)
               end
            end
            if (mod(i,7) == 4)
               mingx = min(mingx,1.1*min(ypkk));
               maxgx = max(maxgx,1.1*max(ypkk)); 
               plot(1000.*xpkk*dt,ypkk,char(splt(4)),'LineWidth',2)
               axis([0 ap_before+ap_after mingx maxgx])
               if (k == num_group_elements(i))
                  ypkk_average = ypkk_average/double(num_group_elements(i));
                  plot(1000.*xpkk*dt,ypkk_average,'k-','LineWidth',2)
               end
            end
            if (mod(i,7) == 5)
               mingx = min(mingx,1.1*min(ypkk));
               maxgx = max(maxgx,1.1*max(ypkk)); 
               plot(1000.*xpkk*dt,ypkk,char(splt(5)),'LineWidth',2)
               axis([0 ap_before+ap_after mingx maxgx])
               if (k == num_group_elements(i))
                  ypkk_average = ypkk_average/double(num_group_elements(i));
                  plot(1000.*xpkk*dt,ypkk_average,'k-','LineWidth',2)
               end
            end
            if (mod(i,7) == 6)
               mingx = min(mingx,1.1*min(ypkk));
               maxgx = max(maxgx,1.1*max(ypkk)); 
               plot(1000.*xpkk*dt,ypkk,char(splt(6)),'LineWidth',2)
               axis([0 ap_before+ap_after mingx maxgx])
               if (k == num_group_elements(i))
                  ypkk_average = ypkk_average/double(num_group_elements(i));
                  plot(1000.*xpkk*dt,ypkk_average,'k-','LineWidth',2)
               end
            end
 
       end
   end
      
end
   
if (createplot == 1)
   highres('overlay') 
end
fprintf('Press enter to continue \n'); 
pause
hold off
subplot(1,1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End: Plot similar spikes. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting selected colored-coded spikes interspersed among time range
% of spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
fprintf('Plotting time range and spike types \n');
fprintf('Number of types of spikes %d \n',number_groups)

distmax_min = max(xpk) - min(xpk);            
plot(xpk,ypk,'c-','LineWidth',2)
%axis([min(xpk) max(xpk) min(ypk) max(ypk)])
%axis([min(xpk) min(xpk)+.02*distmax_min min(ypk) max(ypk)])
axis([min(xpk) min(xpk)+2.0 min(ypk) max(ypk)])
xlabel('Time (seconds)','FontSize',16)
ylabel('Amplitude','FontSize',16)
set(gca,'linewidth',2)
set(gca,'FontSize',12)
hold on
          


freq = zeros(1000,number_groups);
    
ijkt = 0;

for i = 1:number_groups

    for k = 1:num_group_elements(i)
              
        kk = group_elements(i,k);
        % beginpeak firstneg negpeak secondneg firstzero beghalfheight peak endhalfheight endpeak      
        distancebeginpeak = peak(kk) - beginpeak(kk);
               
        ijk = 0;
        xpkk = zeros(endpeak(kk)-beginpeak(kk)+1,1);
        ypkk = zeros(endpeak(kk)-beginpeak(kk)+1,1);

        for j = beginpeak(kk):endpeak(kk)
            ijk = ijk +1;
            xpkk(ijk) = (ijk+beginpeak(kk))*dt;
            ypkk(ijk) = ya(j,1);
            if (i == 1)
               ijkt = ijkt + 1;
               xpkks(ijkt) = xpkk(ijk);
               ypkks(ijkt) = ypkk(ijk);
            end
        end
              
        if (mod(i,7) == 1)
           plot(xpkk,ypkk,char(splt(1)),'LineWidth',2)
           hold on
        end
        if (mod(i,7) == 2)
           plot(xpkk,ypkk,char(splt(2)),'LineWidth',2)
           hold on
        end
        if (mod(i,7) == 3)
           plot(xpkk,ypkk,char(splt(3)),'LineWidth',2)
           hold on
        end
        if (mod(i,7) == 4)
           plot(xpkk,ypkk,char(splt(4)),'LineWidth',2)
           hold on
        end
        if (mod(i,7) == 5)
           plot(xpkk,ypkk,char(splt(5)),'LineWidth',2)
           hold on
        end
        if (mod(i,7) == 6)
           plot(xpkk,ypkk,char(splt(6)),'LineWidth',2)
           hold on
        end 
 
    end

end
   
if (createplot == 1)
   highres('color_coded_spikes') 
end

fprintf('Press enter to continue \n');
pause
fprintf('\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End: Plotting selected colored-coded spikes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating interspike intervals and frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n')
fprintf('Calculating interspike intervals and frequencies \n')
freqmin = 10.; % Minimum frequency
freqmax = 120.; % Maximum frequency
fprintf('Minimum frequency used %d \n',freqmin);
fprintf('Maximum frequency used %d \n',freqmax);
for i = 1:number_groups

    numfreq = 0;
    kkij = 0;
    number_isi(i) = 0;
    for k = 1:num_group_elements(i)-1
        kk = group_elements(i,k);
        kkp = group_elements(i,k+1);
        idif = peak(kkp)-peak(kk);
        spiketime = double(idif)*dt;
        freq(k,i) = 1.d0/spiketime;
        numfreq = numfreq + 1;
           
        if (spiketime >= 1./freqmax && spiketime <= 1./freqmin )
           kkij = kkij + 1;
           xisi1(kkij,i) = .5*(peak(kkp)+peak(kk))*dt;
           yisi1(kkij,i) = spiketime;
           yisif(kkij,i) = 1.0/spiketime;
           number_isi(i) = kkij;
        end
    end 
       
    fprintf('Total number of frequencies %d for group %d \n',numfreq,i);
         
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End: Calculating interspike intervals and frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting four window frame for each spike type
% Upper left window pane: Plot of Spike type interspersed with time interval
% Lower left window pane: Plot of interspike interval over time
% Upper right window pane: Distribution of firing frequencies
% Lower right window pane: Distribution of interspike intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interval_freq = 10;
   
xfreq = zeros(interval_freq,1);
hold off
num_skip = 3;  
fprintf('\n')
fileID = fopen('freq_dist','w');
for kk = 1:number_groups
    
    if (number_isi(kk) > 1)
       fprintf('Creating plots for group %d \n',kk); 
    
       isimin = 1.d+20;
       isimax = -1.d+20;
       isimean = 0.0;
       isimean2 = 0.0;
       for ij = 1:number_isi(kk)
           isimin = min(isimin,yisi1(ij,kk));
           isimax = max(isimax,yisi1(ij,kk));
           isimean = isimean + yisi1(ij,kk);
           isimean2 = isimean2 + yisi1(ij,kk)^2;
       end
       fnum = double(number_isi(kk));
       isimean = isimean/fnum;
       standard_deviation = sqrt((isimean2 - (isimean^2)/fnum)/(fnum-1));
       fprintf('Standard deviation of Interspike Interval scores for Group %d is %5.2f ms \n',kk,1000.*standard_deviation);
       disi = (isimax - isimin)/double(interval_freq);
       binisi = zeros(interval_freq,1);
       ijbinisi = 0;
       for ijk = 1:interval_freq
           isi_low = isimin + double(ijk-1)*disi;
           isi_high = isimin + double(ijk)*disi;
           xisi(ijk) = isimin + (double(ijk)-.5d0)*disi;
           for ij = 1:number_isi(kk)    
               if (yisi1(ij,kk) >= isi_low && yisi1(ij,kk) < isi_high)
                  binisi(ijk) = binisi(ijk) + 1.0;
                  ijbinisi = ijbinisi + 1;
               end
           end 
       end
       fprintf('Total number of interspike intervals %d for group %d \n',ijbinisi,kk)

       freqminc = 1.d+10;
       freqmaxc = -1.d+10;
       for ij = 1:num_group_elements(kk)-1
           ff = freq(ij,kk);
           freqminc = min(freqminc,freq(ij,kk));
           freqmaxc = max(freqmaxc,freq(ij,kk));
       end
       freqmin = max(freqmin,freqminc);
       freqmax = min(freqmax,freqmaxc);

       
       dfreq = (freqmax-freqmin)/double(interval_freq);
       ijbin = 0;
       binfreq = zeros(interval_freq,1);
       for ijk = 1:interval_freq
           freq_low = freqmin + double(ijk-1)*dfreq;
           freq_high = freqmin + double(ijk)*dfreq;
           xfreq(ijk) = freqmin + (double(ijk)-.5d0)*dfreq;
           for ij = 1:num_group_elements(kk)-1      
               if (freq(ij,kk) >= freq_low && freq(ij,kk) < freq_high)
                   binfreq(ijk) = binfreq(ijk) + 1;
                   ijbin = ijbin + 1;
               end
           end 
       end
    
       sts = strcat(char(splt(kk)),'o');

       xisiif = zeros(number_isi(kk),1);
       yisiif = zeros(number_isi(kk),1);   
       for ijka = 1:number_isi(kk)
         xisiif(ijka) = xisi1(ijka,kk);
         yisiif(ijka) = yisif(ijka,kk);
         
         %xisiia = xisii(ijka);
         %yisiia = yisii(ijka);
         
       end    
       %plot(xisii,yisii,char(splt(kk)))
       aaa = min(xisiif);
       bbb = max(xisiif);
       ccc = min(yisiif);
       ddd = max(yisiif);
      
       subplot(2,2,2)
       minxfreq = min(xfreq);
       maxxfreq = max(xfreq);
       distxfreq = maxxfreq - minxfreq;
       maxbin = max(binfreq);
       bar(xfreq,binfreq,char(spltb(kk)))
       axis([minxfreq-.1*distxfreq maxxfreq+.1*distxfreq 0 1.1*maxbin])
       xlabel('Frequency (Hz)','FontSize',14)
       ylabel('Counts','FontSize',14)
       hold off
       
       if (kk == 2)
           for ijk = 1:interval_freq
               xxa = xfreq(ijk);
               yya = binfreq(ijk);
               fprintf(fileID,'%5.1f %4.0f \n',xxa,yya);
           end
       end
       
       xisii = zeros(number_isi(kk),1);
       yisii = zeros(number_isi(kk),1);
       

       xisii_avg = zeros(number_isi(kk),1);
       yisii_avg = zeros(number_isi(kk),1); 
       freq_avg = zeros(number_isi(kk),1);
       
       for ijka = 1:number_isi(kk)
         xisii(ijka) = xisi1(ijka,kk);
         yisii(ijka) = yisi1(ijka,kk);

         num_avg = 0;
         xisii_avg(ijka) = 0.0;
         yisii_avg(ijka) = 0.0;
         for ijkb = max(ijka-num_skip,1):min(ijka+num_skip,number_isi(kk))
             xisii_avg(ijka) = xisii_avg(ijka) + xisi1(ijkb,kk);            
             yisii_avg(ijka) = yisii_avg(ijka) + yisi1(ijkb,kk);
             num_avg = num_avg + 1;
         end
         xisii_avg(ijka) = xisii_avg(ijka)/double(num_avg);
         yisii_avg(ijka) = yisii_avg(ijka)/double(num_avg);
         freq_avg(ijka) = 1.d0/yisii_avg(ijka);

       end
       
       ijkc = 0;
       for ijka = num_skip+1:num_skip:number_isi(kk)-num_skip
           ijkc = ijkc + 1;
       end
       
       xisiib_avg = zeros(ijkc,1);
       yisiib_avg = zeros(ijkc,1);
       freqb_avg = zeros(ijkc,1);
       
       ijkc = 0;
       for ijka = num_skip+1:num_skip:number_isi(kk)-num_skip
           ijkc = ijkc + 1;
           xisiib_avg(ijkc) = xisii_avg(ijka);
           yisiib_avg(ijkc) = yisii_avg(ijka);
           freqb_avg(ijkc) = freq_avg(ijka);
       end
       
       %plot(xisii,1000*yisii,char(splt(kk)))
       subplot(2,2,3)
       plot(xisii,1000*yisii,sts)
       axis([min(xisii) max(xisii) 1000.*min(yisii) 1000.*max(yisii)])
       xlabel('Time (seconds)','FontSize',14)
       ylabel('Interspike Interval (ms)','FontSize',12)
       hold on
       if (ijkc > 1)
          plot(xisiib_avg,1000*yisiib_avg,'k-','Linewidth',3)
          axis([min(xisii) max(xisii) 1000.*min(yisii) 1000.*max(yisii)])
          %legend('','Average')
          %legend boxoff
       end
      
       hold off
       
       subplot(2,2,1)     
       plot(xisiif,yisiif,sts)
       axis([aaa bbb ccc ddd])
       xlabel('Time (seconds)','FontSize',14)
       ylabel('Frequency (Hz)','FontSize',12)
       hold on
       
       if (ijkc > 1)
          plot(xisiib_avg,freqb_avg,'k-','Linewidth',3)
          %legend('','Average')
          %legend boxoff
          axis([aaa bbb ccc ddd])
       end
       hold off
            
       subplot(2,2,4)
       bar(1000*xisi,binisi,char(spltb(kk)))
       distxim = max(1000*xisi) - min(1000*xisi);
       axis([min(1000*xisi)-.1*distxim max(1000*xisi)+.1*distxim 0 1.1*max(binisi)])
       xlabel('Interspike Interval (ms)','FontSize',14)
       ylabel('Counts','FontSize',14)
       hold off
       
       fprintf('Plotting information for group %d \n',kk);
       
       if (createplot == 1) 
          if (kk == 1)
             highres('windowpane1')
          end
          if (kk == 2)
             highres('windowpane2')
          end
          if (kk == 3)
             highres('windowpane3')
          end
          if (kk == 4)
             highres('windowpane4')
          end
       end
    
       if (kk < number_groups)
          fprintf('Press enter to continue \n');
          pause
       end
       
    else
        
        fprintf('Too few spikes to construct interspike interval and frequency histograms for group %d \n',number_groups)
        
    end  
        
        
end  


