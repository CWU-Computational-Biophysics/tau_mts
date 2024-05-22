% Tau/mpa6 project coarse grain model
% Last updated: 12/11/23
clear

% define the name of the .mat file to save
save_str = "T0_5_2";


% Grid parameters
N=1000; % Total number of time steps
M=20; % Number of spatial grid points
dt=0.001; % Time step (s)
L=20; % Domain length
dx=L/(M-1); % Grid spacing (conceptualize as 36nm corresponding to length occupied by tau or map6)

% Initialize variables. 
MTgrid=zeros(N,M); % time, distance, Values: 0=empty, 1=tau, 2=map6

t=0:dt:(N-1)*dt;
x=0:dx:L;



growthstate=0; % 0=static,-1=shrinking, 1=growth

MTlength=zeros(numel(t),3);
MTlength(:,1)=t;
MTlength(1,2)=L;
MTlength(1,3)=growthstate;

% Tunable parameters
T0=5.0; %base binding rate for tau
Toff=5.0; %undbinding rate for tau
M0=0.1; %base binding rate for map6
Moff=0.1; %unbinding rate for map6
alphaT=0.0; % cooperativity for tau
alphaM=0.0; %cooperativity for map6
%Tinit=0.5; % initial fraction of binding sites on labile MT segment with tau bound
%Minit=0.5; % initial fraction of binding sites on stable MT segment with map6 bound
fmp=500; % rescue frequency
fpm=300; % catastrophe frequency

if (T0+M0)*dt>1MT
    'on rates too high'
end


for j=1:N-1 %loop over time

 
    
    % enforce tau and map cooperativity
    taunext=0;
    mapnext=0;
    if(MTgrid(j,2)==1)
        taunext=taunext+1;
    elseif(MTgrid(j,2)==2)
        mapnext=mapnext+1;
    end

    Toffeff=Toff-alphaT*Toff*taunext;
    Moffeff=Moff-alphaM*Moff*mapnext;

    % Update occupancy for i=1

    if(MTgrid(j,1)==0) % site un-occupied
        if rand<T0*dt % chance for tau binding
            MTgrid(j+1,1)=1; 
        elseif T0*dt<rand<(T0+M0)*dt; % chance for map binding
            MTgrid(j+1,1)=2;
        else
            MTgrid(j+1,1)=0; % stays unoccupied
        end
    end

    if(MTgrid(j,1)==1)
        if rand<Toffeff*dt
            MTgrid(j+1,1)=0;
        else
            MTgrid(j+1,1)=1;
        end
    end

    if(MTgrid(j,1)==2)
        if rand<Moffeff*dt
            MTgrid(j+1,1)=0;
        else
            MTgrid(j+1,1)=2;
        end
    end


    for i=2:M-1 %loop over x

       % determine Ton based on next neigbor occupancy
       
       taunext=0;
       mapnext=0;
            
       if(MTgrid(j,i-1)==1)
           taunext=taunext+1;
       elseif(MTgrid(j,i-1)==2)
           mapnext=mapnext+1;
       end
      
       if(MTgrid(j,i+1)==1)
           taunext=taunext+1;
       elseif(MTgrid(j,i+1)==2)
           mapnext=mapnext+1;
       end

       %Ton=T0+alphaT*taunext;
       Toffeff=Toff-alphaT*Toff*taunext;
       Moffeff=Moff-alphaM*Moff*mapnext;
       %%%%%%

       % Update occupancy numbers
        
       if(MTgrid(j,i)==0) % site un-occupied
           if rand<T0*dt % chance for tau binding
               MTgrid(j+1,i)=1; 
           elseif T0*dt<rand<(T0+M0)*dt; % chance for map binding
               MTgrid(j+1,i)=2;
           else
               MTgrid(j+1,i)=0; % stays unoccupied
           end
       end

       if(MTgrid(j,i)==1)
           if rand<Toffeff*dt
               MTgrid(j+1,i)=0;
           else
               MTgrid(j+1,i)=1;
           end
       end

       if(MTgrid(j,i)==2)
           if rand<Moffeff*dt
               MTgrid(j+1,i)=0;
           else
               MTgrid(j+1,i)=2;
           end
       end

    end


    % Ton for i=M
    taunext=0;
    mapnext=0;
    if(MTgrid(j,M-1)==1)
        taunext=taunext+1;
    elseif(MTgrid(j,M-1)==2)
        mapnext=mapnext+1;
    end
    %Ton=T0+alphaT*taunext;
    Toffeff=Toff-alphaT*Toff*taunext;
    Moffeff=Moff-alphaM*Moff*mapnext;

    % Update occupancy for i=M

    if(MTgrid(j,M)==0) % site un-occupied
        if rand<T0*dt % chance for tau binding
            MTgrid(j+1,M)=1; 
        elseif T0*dt<rand<(T0+M0)*dt; % chance for map binding
            MTgrid(j+1,M)=2;
        else
            MTgrid(j+1,M)=0; % stays unoccupied
        end
    end

    if(MTgrid(j,M)==1)
        if rand<Toffeff*dt
            MTgrid(j+1,M)=0;
        else
            MTgrid(j+1,M)=1;
        end
    end

    if(MTgrid(j,M)==2)
        if rand<Moffeff*dt
            MTgrid(j+1,M)=0;
        else
            MTgrid(j+1,M)=2;
        end
    end

    % chance to grow

    if(sum(MTgrid(j,M-2:M))>0) % dynamic MTs when tau bound anywhere in last three sub-units
        
        % chance to change growth state
        if growthstate==0
            growthstate=1;
        elseif growthstate==-1
            if rand<fmp*dt
                growthstate=1;
            end
        elseif growthstate==1
            if rand<fpm*dt
                growthstate=-1;
            end
        end
        
        % enact growth or shrinking

        if growthstate==1 % growing
         
            M=M+1;
            L=L+dx;

            dx=L/(M-1); % Grid spacing (conceptualize as 36nm corresponding to length occupied by tau or map6)
            x=0:dx:L;

            x=0:dx:L;

            MTgrid=[MTgrid,zeros(N,1)];
            MTgrid(1:j,M)=-1; % -1 means 'location does not exist at this time'
            MTgrid(j+1,M)=1;

        elseif growthstate==-1 % shrinking
           % M=M-1;
            L=L-dx;
            %x=0:dx:L;

            MTgrid(j+1,1:M-1)=MTgrid(j,1:M-1);
            MTgrid(j+1,M)=-1; % location no longer exists
        end



        %MTgrid(j+1,:)=[MTgrid(j,:),1];
       % MTgrid(j+1,1:M-1)=MTgrid(j,1:M-1);
       % MTgrid(j+1,M)=1; % new site has tau
    else 
        growthstate=0; %static
    end

    MTlength(j+1,2)=L;
    MTlength(j+1,3)=growthstate;
end

% figure
% bar(x',MTgrid(N,:))
% 
% figure
% plot(t',MTlength(:,2))
% title('MT length vs time')
% 
% figure
% plot(t',MTlength(:,3))
% title('growth state vs time')

% Calculate binding fraction asymmetry
% defined as tau fraction, (tau/(tau+map6)), at plus end (final five
% positions) vs. tau fraction along the length

taufractip=sum(MTgrid(N,M-10:M)==1)/numel(MTgrid(N,M-10:M));
taufraclength=sum(MTgrid(N,1:M-11)==1)/numel(MTgrid(N,1:M-11));

tauplusendasym=taufractip/taufraclength;

mapfractip=sum(MTgrid(N,M-10:M)==2)/numel(MTgrid(N,M-10:M));
mapfraclength=sum(MTgrid(N,1:M-11)==2)/numel(MTgrid(N,1:M-11));

mapplusendasym=mapfractip/mapfraclength;

% L


% figure
% fig = figure('units','inch','position',[10,10,15,4]);
% for i=1:M
%         if(MTgrid(N,i)==0) % nothing bound
%             bar(i,1,FaceColor="#77AC30")
%         elseif(MTgrid(N,i)==1) % tau bound
%            bar(i,1,FaceColor="#7E2F8E")
%         elseif(MTgrid(N,i)==2) % map6 bound
%             bar(i,1,FaceColor="#0072BD")
%         elseif(MTgrid(N,i)==-1) % no MT here
%             bar(i,1,'w',EdgeColor='w')
%         end
%         hold on
%  end


% movie=1;
% if movie ==1
% 
% % Initialize video
% myVideo = VideoWriter('MT_tau_map6_121523'); %open video file (change file name to avoid overwriting)
% myVideo.FrameRate = 10;  %can adjust this, 5 - 10 works well for me
% open(myVideo)
% 
% figure
% fig = figure('units','inch','position',[10,10,15,4]);
% for j=1:N
%     for i=1:M
%         if(MTgrid(j,i)==0) % nothing bound
%             bar(i,1,FaceColor="#77AC30")
%         elseif(MTgrid(j,i)==1) % tau bound
%            bar(i,1,FaceColor="#7E2F8E")
%         elseif(MTgrid(j,i)==2) % map6 bound
%             bar(i,1,FaceColor="#0072BD")
%         elseif(MTgrid(j,i)==-1) % no MT here
%             bar(i,1,'w',EdgeColor='w')
%         end
%         hold on
%     end
%     hold off
%     pause(0.1)
%     frame=getframe(gcf);
%     writeVideo(myVideo,frame);
% end
% close(myVideo)
% end

% export the workspace
% define the export location
export_dir = fullfile("data/raw_sim_data");
export_file = fullfile(export_dir, save_str);
[~, ~, ~] = mkdir(export_dir);
save(export_file, "-nocompression", "-v7");