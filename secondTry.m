close all;
test1=readtable('ScanFilesNew\Scan000000.dat');
test2=readtable('ScanFilesNew\Scan000080.dat');
X=table2array(test1);
P=table2array(test2);



%get dot cloud
theta=deg2rad(0:180);
[x,y]=pol2cart(theta.',X);
%cartesian
X_k=[x,y];
X_k=X_k.';
[x,y]=pol2cart(theta.',P);
P_k=[x,y];
P_k=P_k.';


% % test wolken
% vector=1:20;
% X_k=[vector;vector];
% pz=5*ones(1,20);
% P_k=[pz;vector]

% X_k=[-10 0 0 10 0;0 -10 10 0 0];
% a=60;
% rot=[cosd(a) -sind(a);sind(a) cosd(a)];
% P_k=[-10 0 0 10 1;0 -10 10 0 1];
% 
%     for i=1:4
%         point=P_k(:,i);
%         rPoint=rot*point;
%         P_k(:,i)=rPoint;
%     end

 %paint
 paint_kartesian_image(X_k,'k');
 hold on;
 paint_kartesian_image(P_k,'m');
 hold on;

[R,t]=icp(X_k,P_k,1);

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
       
        [R,t]=get_optimal_R_and_t(M_k,transpX_k,mue_x,triang);
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

    R=totalRotation;
    t=totalTranslation;
    equal=verify(M_k,P_k,R,t,paint);
    if(~equal)
        error('There was a problem while calculating the new point cloud.')
    end
 
end

function [R,t]=get_optimal_R_and_t(M_k,transpX_k,mue_x,triang)

        mue_p=calc_mue(M_k);
        %get closest points
        %gives indices of the closest points in X_k for each point in M_k
        closest=dsearchn(transpX_k,triang,M_k.');
        %Kreuzterm der Kostenfunktion
        W=zeros(2);
        for i=1:length(M_k)
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

function [M]=new_point_cloud(M_k,R,t)
 %Neue Punktwolke
        for i=1:length(M_k)
            point=M_k(:,i);
            rPoint=R*point;
            M_k(:,i)=rPoint+t;
        end
        M=M_k;
end

function [mue]=calc_mue(M_k)
    mue=sum(M_k.')/length(M_k);
    mue=mue.';
end

function paint_kartesian_image(X_k, color)
    x = X_k(1, :);
    y = X_k(2, :);
    scatter(x, y,color);
end

function [equal]=verify(M_k,P_k,R,t,paint)

    %verify result
    if(paint)
        figure('Name','Solution');
        paint_kartesian_image(M_k,'c');
        hold on;
        result=P_k;
        for i=1:length(P_k)
                point=P_k(:,i);
                rPoint=R*point;
                result(:,i)=rPoint+t;
        end
        paint_kartesian_image(result,'m');
        hold on
    end
    equal=isequal(result,M_k);
end
