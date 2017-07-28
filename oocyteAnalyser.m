 % This script reads in the tracking results by plusendtrack
% and analyses an oocyte

%output will be saved in roi*/track/OOCYTEoutput.m

%Output comprises the following parameters:
% Angle_X: All frame2frame angles in each compartment (compartments divide the oocyte in areas of 10um)
% AvAngle_X: The average angle for each compartment
% AvAngle_trackX: The end2end angles of each track in each compartment
% AvMagnitude_X: The magnitude of the Av Vector in each compartment
% fracPostAnt: The fraction of posterior oriented MTs in each compartment
% track: information about each individual track


% by Max Jakobs 23/10/2015

clear all

close all

%segmentation parameter in um
Parameters.domainsize=15;

%% load the ROI mask specified during plusTipTrack
topfolder = uigetdir('Select TOP folder containing all projects allready analysed with plusTipTrack (two folders above roi_* folders))'); 
cd(topfolder);
d=dir;
b=[d.isdir];
directories=d(b);
directories=directories(3:end);
% finish loading

%% the following is the main loop and analyzes ALL directories in the folder that have roi_* in their subfolders

%keep in mind that q runs over all folders above the roi and p loops
%through all the roi subfolders
disp('Watch out that the roi subfolders start at roi_1 and increase ony by one!')
clear d b
for q=1:length(directories) 
    p=1;
    path=[topfolder '/' directories(q).name '/roi_' int2str(p) '/'];
    if (p==1 && ~isdir(path))
        disp([path,' Does not exist even though it should. Maybe you selected the wrong topfolder?']);
        if ( input('Omit this folder and continue?(yes/no)','s')=='yes')
            continue;
        else
            error('exiting...');
        end
    end
    %% this while loop analyzes EACH roi individually, i.e. averything from  now on is part of one analyzation 
   while isdir(path)
        
        %read the mask from plusTipTracks
        chosenfile  = [path '/roiMask.tif']; 
        mask        = imread(chosenfile); 
        clear file chosenfile;

        %load the trackFinal struct from plusendtrack
        chosenfile  = [path 'track/trackResults.mat' ]; 
        load(chosenfile);

        %load projData
        chosenfile  = [path 'meta/projData.mat' ]; 
        load(chosenfile);
        
        %load pixelsize and dt
        Parameters.pixelsize=projData.pixSizeNm*0.001;
        Parameters.dt=projData.secPerFrame;
        
         %% prepare track data in a neat struct and get frame2frame velocities and accelerations
        % total number of tracks
        ntrack=length(tracksFinal);
        %load coordinates from projData
        for i=1:ntrack
            j=0;
            while j*8<length(tracksFinal(i).tracksCoordAmpCG)
                track(i).coord(j+1,1)=tracksFinal(i).tracksCoordAmpCG(j*8+1);
                track(i).coord(j+1,2)=tracksFinal(i).tracksCoordAmpCG(j*8+2);
                j=j+1;
            end
        end

        %% load the full image to be cropped by the mask

        imgpath        = [path '../images/']; 
        [listOfImages] = searchFiles('.tif',[],imgpath,0);
        chosenfile     = strcat(imgpath , listOfImages(1,1));
        disp(['Processing: ',path]);
        img            = imadjust(imread(chosenfile{1,1}));
        
        clear file chosenfile listOfImages imgpath;
        
        
        %% Draw oocyte axis to analyze directionality of comets
        
        % crop image according to mask inherited from plusTipTracker
        ROI=img;
        ROI(mask==0)=0;

        %check if drawn axis exists and if not ask user to draw
        file=[path 'track/axis.mat'];
        if exist(file,'file')~=2
            %ask user to draw an axis in the cropped image for 1D anallysis
            axis.coords=ginputExtra(ROI);
            
            if length(axis.coords)>2
                error('only 2 points please!')
            end
            %array of individual length segments of oocyte
            dx=axis.coords(2,1)-axis.coords(1,1);
            dy=axis.coords(2,2)-axis.coords(1,2);
            axis.length=sqrt(dx*dx+dy*dy);
            
            if axis.length<(Parameters.domainsize/Parameters.pixelsize)
                error('drawn axis shorter that domainsize!!! NOT OK!');
            end
            
            clear dx dy;
            
            
            %save oocyte axis data in track/axis.mat
            file=[path 'track/axis.mat'];
            save(file,'axis');
        
        % if axis was drawn previously load this one instead    
        else
            load([path 'track/axis.mat']);
        end
        
        
        
        %theta is the angle of the drawn axis
        theta=-atan((axis.coords(2,2)-axis.coords(1,2))/(axis.coords(2,1)-axis.coords(1,1)));

        %% get position and movement angles of granules along the drawn neurite
        
        % av track pos and angles at the same time and correct them by
        % theta so that the angles are with respect to the drawn axis!!
        for trk=1:length(track)
            track(trk).AvAngl=atan2(-(track(trk).coord(end,2)-track(trk).coord(1,2)),track(trk).coord(end,1)-track(trk).coord(1,1))-theta;
            if track(trk).AvAngl<-pi
                track(trk).AvAngl=track(trk).AvAngl+pi;
            elseif track(trk).AvAngl>pi
                track(trk).AvAngl=track(trk).AvAngl-pi;
            end
            track(trk).angles=atan2(-(track(trk).coord(2:end,2)-track(trk).coord(1:end-1,2)),track(trk).coord(2:end,1)-track(trk).coord(1:end-1,1))-theta;
            for i=1:length(track(trk).angles)
                if track(trk).angles(i)<-pi
                    track(trk).angles(i)=track(trk).angles(i)+pi;
                elseif track(trk).angles(i)>pi
                    track(trk).angles(i)=track(trk).angles(i)-pi;
                end
            end
        end
        clear trk tmp dx dy dist file;
        
         %% angle data
        % first bin the angles according to parameters.Nseg
        dx=Parameters.domainsize/Parameters.pixelsize;
        
        x_f2f=[];
        x_AV=[];
        tot_angle=[];
        z=[];
        directions=[];
        
        Oax=sort((cos(theta)*axis.coords(2,1)-sin(theta)*axis.coords(2,2))-dx/2:-dx:(cos(theta)*axis.coords(1,1)-sin(theta)*axis.coords(1,2)));
        
     
        
        for i=1:ntrack
            
            x_AV=[x_AV; cos(theta)*(track(i).coord(end,1)+track(i).coord(1,1))/2-sin(theta)*(track(i).coord(end,2)+track(i).coord(1,2))/2];
            x_f2f=[x_f2f; cos(theta).*track(i).coord(1:end-1,1)-sin(theta).*track(i).coord(1:end-1,2)];
            z=[z, track(i).AvAngl];
            
            
      
            %AV angle
            if abs(track(i).AvAngl)<pi/2 
                track(i).AVdirect=1;
            elseif abs(track(i).AvAngl)>=pi/2 && abs(track(i).AvAngl)<=pi   
                track(i).AVdirect=0;    
            else
                error()
            end
            
            directions=[directions,track(i).AVdirect];
            
            for j=1:length(track(i).angles)
                tot_angle=[tot_angle,track(i).angles(j)];
                if abs(track(i).angles(j))<pi/2 
                    track(i).direct(j)=1;
                elseif abs(track(i).angles(j))>=pi/2 && abs(track(i).angles(j))<=pi
                    track(i).direct(j)=0;    
                else
                    error()
                end
            end
                
            
            %error messages
            if length(tot_angle)~=length(x_f2f)
                disp(i)
                disp(ntrack)
                disp(length(tot_angle))
                disp(length(x_f2f))
                error();
            elseif length(directions)~=length(x_AV)
                disp(i)
                disp(ntrack)
                disp(length(directions))
                disp(length(x_AV))
                error();
            elseif length(z)~=length(x_AV) 
                disp(i)
                disp(ntrack)
                disp(length(z))
                disp(length(x_AV))
                error();
            end
        end
        
        [~,~,~,Angle_X]=bindata(x_f2f,tot_angle,Oax);
        [~,~,~, Direct]=bindata(x_AV,directions,Oax);
        [~,~,~, AvAngle_trackX]=bindata(x_AV,z,Oax);
        
        
        
        x=cell([length(Angle_X) 1]);
        y=cell([length(Angle_X) 1]);
        AvAngle_X=cell([length(Angle_X) 1]);
        AvMagnitude_X=cell([length(Angle_X) 1]);
        %convert to kartesian
        for i=1:length(Angle_X)
            x{i}=cos(Angle_X(i).dat);
            y{i}=sin(Angle_X(i).dat);
            AvAngle_X{i}=atand(mean(y{i})/mean(x{i}));
            AvMagnitude_X{i}=sqrt(mean(y{i})^2+mean(x{i})^2);
        end
        
         %fraction of Posterior comets in each compartment
         for i=1:length(Direct)
             fracPostAnt{i}=sum(Direct(i).dat)/length(Direct(i).dat);
         end
            
       
        

        %% save to matlab file in the folder where the track struct is from
        save([path 'track/OOCYTEoutput.mat'],'Direct','track','Parameters','Angle_X','AvAngle_trackX','AvAngle_X','AvMagnitude_X','theta','fracPostAnt');
        disp(['Output saved in: ', path, 'track/OOCYTEoutput.mat'])
        clearvars -except p directories topfolder q Parameters 
        p=p+1;
        path=[topfolder '/' directories(q).name '/roi_' int2str(p) '/'];
   end
end

clear all;


