clc; clear all; close all;

%User Inerface - File path request                                          *Anleitung how to use
prompt = 'Provide user info or default. Enter 1 for default.\n Any other for user info ';
x = input(prompt);
if(x==1)
    dir_to_search = 'ScanFilesNew\';
    txtpattern = fullfile(dir_to_search, '*.dat');
else
    prompt = 'Provide folder path of the files.\n The files in the folder need to be in the right order \n Enter to continue';
    x = input(prompt);
    folder=uigetdir();
    str=input('What is the file ending? (eg .dat). Files need to be tables','s');
%     str = input(prompt);
    s = strcat('*.',str);
    txtpattern = fullfile(folder,s);
end
%init
theta=deg2rad(0:180);
files_info=dir(txtpattern);
tic
%itterate over all datafiles.
i=1;
skipsize=20;                                                               

%ICP main-loop
while i+skipsize<length(files_info)
   filename = fullfile(dir_to_search, files_info(i).name);
   X=readtable(filename);
   filename = fullfile(dir_to_search, files_info(i+skipsize).name);
   P=readtable(filename);
   X=table2array(X);
   P=table2array(P);
   X_k=get_cartesian_matrix(X,theta);
   P_k=get_cartesian_matrix(P,theta);
   [R,t]=icp(X_k,P_k,0);
   i=i+skipsize;
end

timet=toc;

% test1=readtable('ScanFilesNew\Scan000000.dat');
% test2=readtable('ScanFilesNew\Scan000080.dat');
% X=table2array(test1);
% P=table2array(test2);
% 
% 
% 
% X_k=get_cartesian_matrix(X,theta);
% P_k=get_cartesian_matrix(P,theta);
% 
% % % test wolken
% % vector=1:20;
% % X_k=[vector;vector];
% % pz=5*ones(1,20);
% % P_k=[pz;vector]
% 
% % X_k=[-10 0 0 10 0;0 -10 10 0 0];
% % a=60;
% % rot=[cosd(a) -sind(a);sind(a) cosd(a)];
% % P_k=[-10 0 0 10 1;0 -10 10 0 1];
% % 
% %     for i=1:4
% %         point=P_k(:,i);
% %         rPoint=rot*point;
% %         P_k(:,i)=rPoint;
% %     end
% 
%  %paint
%  paint_kartesian_image(X_k,'k');
%  hold on;
%  paint_kartesian_image(P_k,'m');
%  hold on;
% 
% [R,t]=icp(X_k,P_k,1);

%transverts polar to cartesian coordinates of a supplied point-cloud
function X_k=get_cartesian_matrix(X,theta)
    [x,y]=pol2cart(theta.',X);
    %cartesian
    X_k=[x,y];
    X_k=X_k.';
end

%gives the rotation matrix and the translation vector for matching one
%point cloud P_k to another X_k. The Matrices are in cartesian coorinates
%size 2xN. N is the amount of data points.
function [R,t]=icp(X_k,P_k,paint)

    %rotmatrix [cos -sin; sin cos]

    %M_k is the new Point cloud.
    M_k=P_k;
    %mean value (Mittelwert) for fixed pointcloud
    mue_x=calc_mue(X_k);
    %Delaunay triangulation for faster calculation of closest points
    triang=delaunayn(X_k.');
    %store closest points. if it doesn't change anymore, calculation can be
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

        [R,t]=get_optimal_R_and_t(M_k,transpX_k,mue_x,triang);
        %new point cloud
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

    R=totalRotation;
    t=totalTranslation;
    equal=verify(M_k,P_k,R,t,paint);
    if(~equal)
        error('There was a problem while calculating the new point cloud.')
    end
 
end

%gives optimal rotation matrix and translation vector of two point clouds
function [R,t]=get_optimal_R_and_t(M_k,transpX_k,mue_x,triang)

        mue_p=calc_mue(M_k);
        %get closest points
        %gives indices of the closest points in X_k for each point in M_k
        closest=dsearchn(transpX_k,triang,M_k.');
        %cross term of the cost function
        W=zeros(2);
        for i=1:length(M_k)
            point=transpX_k(closest(i),:);                                 
            x=point.'-mue_x;
            p=M_k(:,i)-mue_p;
            W=W+x*p.';
        end
        %singular value decomposition
        [U,S,V]=svd(W);
        %optimal rotation matrix
        R=U*V.';
        %optimal translation vector
        t=mue_x-R*mue_p;
end

%relocates the point cloud according to the supplied rotation matrix and 
%translation vector.
function [M]=new_point_cloud(M_k,R,t)
 %new point cloud
        for i=1:length(M_k)
            point=M_k(:,i);
            rPoint=R*point;
            M_k(:,i)=rPoint+t;
        end
        M=M_k;
end

%get the mean values of a point cloud
function [mue]=calc_mue(M_k)
    mue=sum(M_k.')/length(M_k);
    mue=mue.';
end

%gives map with point cloud on it
function paint_kartesian_image(X_k, color)
    x = X_k(1, :);
    y = X_k(2, :);
    scatter(x, y,color);
end

function [equal]=verify(M_k,P_k,R,t,paint)

    %verify result
    result=P_k;
    for i=1:length(P_k)
        point=P_k(:,i);
        rPoint=R*point;
        result(:,i)=rPoint+t;
    end
    if(paint)
        figure('Name','Solution');
        paint_kartesian_image(M_k,'c');
        hold on;
        paint_kartesian_image(result,'m');
        hold on
    end
    equal=isequal(result,M_k);                                             
end
