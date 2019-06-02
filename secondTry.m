
%closes all figures that were used befor
close all;

compare=0;
useOwnICP=0;
%user interface
prompt = 'Do you want to use the default config or give your own folder for the scan data?\n Enter 1 for default, anything else for you own input \n';
x = input(prompt);
if(x==1)
    dir_to_search = 'ScanFilesNew\';
    txtpattern = fullfile(dir_to_search, '*.dat');
else
    prompt = 'Provide the folder path that contains the scans.\n The files in the folder need to be in the right order \n Enter to continue';
    x = input(prompt);
    dir_to_search=uigetdir();
    str=input('What is the file ending? (eg. dat). Files need to be tables. Tables in Polarcoordinates and only up to 180 degrees. \n','s');
    s = strcat('*.',str);
    txtpattern = fullfile(dir_to_search,s);
    prompt = 'Compare the matlab icp with own or run icp? \n 1=compare; anything else= run icp \n';
    x= input(prompt);
    if(x==1)
        compare=1;
    else
        prompt = 'Do you want to run my own icp?\n 1= own icp; anything else= matlab icp \n';
        x= input(prompt);
        if(x==1)
            useOwnICP=1;
        end
    end
end


%reads all the fileinfos out of the given directory
files_info=dir(txtpattern);
itterate_over_folder(dir_to_search,files_info,compare,useOwnICP);


%itterates over the files in the given folder
function itterate_over_folder(dir_to_search,files_info,compare,useOwnICP)
    %degree to rad for polar coordinates
    theta=deg2rad(0:180);
    %timer. only 'cause curious.
    tic
    %number indicates with  which can should be started.
    %normaly 1 but can be used if there are any starting issues (eg.
    %calibration of the sensor). Or for debugging purposes
    i=1;
    %If the robot only moves a little bit from scan to scan it is ok to skip
    %some scans to make the calculation faster. But carefull there. Increase
    %only in small steps to ensure still good results.
    skipsize=1;
    %represents the rotation of the actual scan compared with the first scan
    totalRotation=[1 0; 0 1];
    %represents the total translation from the actual scan to the first scan
    totalTranslation=zeros(2,1);

    %translation vector of the matlab icp 
    old_t_mat=[0;0];
    %only used for the comparison between own ipc and matlab ipc
    old_t=[0;0];

    %names the figure
    mainfigure=figure('Name','Robo');
    %only used if the matlab icp is used. 
    oldTform=affine3d;
    while i+skipsize<length(files_info)

       %reads two of the scans
       filename = fullfile(dir_to_search, files_info(i).name);
       X=readtable(filename);
       filename = fullfile(dir_to_search, files_info(i+skipsize).name);
       P=readtable(filename);
       X=table2array(X);
       P=table2array(P);
       %converts the scans from polar to cartesian coordinates
       X_k=get_cartesian_matrix(X,theta);
       P_k=get_cartesian_matrix(P,theta);


       if(compare)
        % compares the own icp with the matlab icp by printing the calculated
        % movement of the robot. cyan is the own calculation and magenta the
        % matlab one.
        [R,t,M]=my_icp(X_k,P_k,0);
        [R_mat,t_mat,matlabM]=matlabicp(X_k,P_k);
        paint_kartesian_image(old_t,'c')
        hold on;
        paint_kartesian_image(old_t_mat,'m')
        hold on;
        old_t=old_t+t;
        old_t_mat=old_t_mat+t_mat;
       else
       %calculates the map 
           if(useOwnICP)
            [R,t,M]=my_icp(X_k,P_k,0);
            %rotation and translation from current scan compared to first scan
            totalRotation=R*totalRotation;
            totalTranslation=totalTranslation+t;
            %coordinates from the actual scan compared  to the first scan
            map=totalRotation.'*P_k+totalTranslation; 
            %paint
            paint_kartesian_image(map,'k')
            hold on;
            paint_kartesian_image(totalTranslation,'c')
            hold on;
           else
            %gives back an homogenous matrix. in it is the rotation and the
            %translation. if used like that we dont have to pay attention if we
            %rotate befor or after the translation
            tform=matlabicp_tform(X_k,P_k);
            oldTform.T=tform.T*oldTform.T;
            %calculates the coordidnates of the current scan compared with the
            %first scan
            new_coordinate=mat_4D_Rot_trans(P_k,oldTform);
            %calculates the position of the robot
            [R_mat,t_mat,matlabM]=matlabicp(X_k,P_k);
            old_t_mat=old_t_mat+t_mat;
            %paint
            paint_kartesian_image(new_coordinate,'k');
            hold on;
            paint_kartesian_image(old_t_mat,'m')
            hold on;

           end
       end

       %if this is on. matlab is fast enough to print the painting and so we
       %see what way the roboter is doing while it is computed.
       %pause(0.1);

       %with this an certaint amout of scans can be skipped. so the calculation gets faster
       %but dangerous if the robot does suddendly big moves -> we gonna miss
       %those.
       i=i+skipsize;
    end
    %timer. I was curious how long it takes. 
    timet=toc; 
    %save the image
    imgname=picture_name(compare,useOwnICP);
    saveas(mainfigure,imgname); 
end

%calculates the filename for the picture
function [str]=picture_name(compare,useOwnICP)
    str='icp.png';
    if(compare)
        str='robot_move_compare.png';
    else
        if(useOwnICP)
                str='ownICP.png';
        else
             str='matlabICP.png';
        end
    end
end

%uses the matlab icp algorithm. gives back the same return values as in the
%own icp.
function [R,t,M_k]=matlabicp(X_k,P_k)
   
    X_3D=[X_k;zeros(1,length(X_k))];
    P_3D=[P_k;zeros(1,length(P_k))];
    ptcloud_x=pointCloud(X_3D');
    ptcloud_p=pointCloud(P_3D.');
    tform=pcregistericp(ptcloud_p,ptcloud_x);
    newcloud=pctransform(ptcloud_p,tform);
    M_k=newcloud.Location(:,1:2);
    M_k=M_k.';
    R=tform.T(1:2,1:2);
    t=tform.T(4,1:2).';
end

%uses matlab icp algorithm. but gives back an tform. this is an homogenous
%matrix with the rotation matrix and the translation vector in it. the
%matlab icp is made for 3 dimensions!
function [tform]=matlabicp_tform(X_k,P_k)  
    %transforms 2D clouds to 3D
    X_3D=[X_k;zeros(1,length(X_k))];
    P_3D=[P_k;zeros(1,length(P_k))];
    %transforms into point clouds cause they are needed for the icp
    ptcloud_x=pointCloud(X_3D');
    ptcloud_p=pointCloud(P_3D.');
    tform=pcregistericp(ptcloud_p,ptcloud_x);
end


%Returns the 2D matrix from P_k (also 2D) that was transformed with the
%tfrom given from the matlabicp
function [M_k]=mat_4D_Rot_trans(P_k,tform)

   P_3D=[P_k;zeros(1,length(P_k))];
   ptcloud_p=pointCloud(P_3D.');
   newcloud=pctransform(ptcloud_p,tform);
   M_k=newcloud.Location(:,1:2);
   M_k=M_k.';
   
end

%transfroms polar coordinates into cartesian coordinates
function X_k=get_cartesian_matrix(X,theta)
    [x,y]=pol2cart(theta.',X);
    %cartesian
    X_k=[x,y];
    X_k=X_k.';
end

%gives the rotation matrix and the translation vector for matching one
%point cloud P_k to another X_k. The Matrices are in cartesian coorinates
%size 2xN. N is the amount of data points.
function [R,t,M]=my_icp(X_k,P_k,paint)

    %rotmatrix [cos -sin; sin cos]

    %M_k is the new Point cloud.
    M_k=P_k;
    %mean value (Mittelwert) for fixed pointcloud
    mue_x=calc_mue(X_k);
    %Delaunay triangulation for faster calculation of closest points
    triang=delaunayn(X_k.');
    %store closest points. if it doesn change anymore calculation can be
    %stopped.
    oldClosest=zeros(length(P_k),1);
    closest=ones(length(P_k),1);
    %variable to store the final rotation an translation
    totalRotation=[1 0; 0 1];
    totalTranslation=zeros(2,1);
    %transpose of X_k is used many times. so store here
    transpX_k=X_k.';
    while(~isequal(oldClosest,closest))
        
        oldClosest=closest;
        %get closest points
        %gives indices of the closest points in X_k for each point in M_k
        closest=dsearchn(transpX_k,triang,M_k.');
        [R,t]=get_optimal_R_and_t(M_k,transpX_k,mue_x,closest);
        %Neue Punktwolke
        M_k=new_point_cloud(M_k,R,t);
        
        if(paint)
            %paint to follow change
            figure('Name','Evolution');
            paint_kartesian_image(X_k,'k');
            hold on;
            paint_kartesian_image(M_k,'m');
            hold on;
        end
       totalRotation=R*totalRotation;
       totalTranslation=totalTranslation+t;
    end

    %retunr values
    R=totalRotation;
    t=totalTranslation;
    M=M_k;
 
end

%returns the rotation matrix and translation vector to transfrom M_k to
%transpX_k
function [R,t]=get_optimal_R_and_t(M_k,transpX_k,mue_x,closest)

        mue_p=calc_mue(M_k);
        %Kreuzterm der Kostenfunktion
        W=zeros(2);
        for i=1:length(closest)
            point=transpX_k(closest(i),:);
            x=point.'-mue_x;
            p=M_k(:,i)-mue_p;
            W=W+x*p.';
        end
        %Singulärwertzerlegung
        [U,S,V]=svd(W);
        %Optimale Rotationsmatrix
        R=U*V.';
        %optimaler Translationsvektor
        t=mue_x-R*mue_p;
end

%calculates the new point cloud of M_k when you use the R and t
function [M]=new_point_cloud(M_k,R,t)
 %Neue Punktwolke
        for i=1:length(M_k)
            point=M_k(:,i);
            rPoint=R*point;
            M_k(:,i)=rPoint+t;
        end
        M=M_k;
end

%calculates mü. for a cartesian 2d point cloud. retunrs a scalar
function [mue]=calc_mue(M_k)
    mue=sum(M_k.')/length(M_k);
    mue=mue.';
end

%paints the matrix to a figure
function paint_kartesian_image(X_k, color)
    x = X_k(1, :);
    y = X_k(2, :);
    scatter(x, y,color);
end


