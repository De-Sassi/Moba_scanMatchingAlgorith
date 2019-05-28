close all;
test1=readtable('ScanFilesNew\Scan000000.dat');
test2=readtable('ScanFilesNew\Scan000020.dat');
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

lnght=length(X_k);
M_k=P_k;


mue_x=calc_mue(X_k);
%Delaunay triangulation for faster calculation of closest points
triang=delaunayn(X_k.');
while(true)

    %Mittelwert der Punktwolke
    mue_p=calc_mue(M_k);
    %get closest points
    %gives indices of the closest points in X_k for each point in M_k
    closest=dsearchn(X_k.',triang,M_k.');
    %Kreuzterm der Kostenfunktion
    W=zeros(2);
    for i=1:lnght
        test=closest(i);
        transp=X_k.';
        point=transp(test,:);
        point=point.';
        x=point-mue_x;
        p=M_k(:,i)-mue_p;
        W=W+x*p.';
    end
    %Singulärwertzerlegung
    [U,S,V]=svd(W);
    %Optimale Rotationsmatrix
    R=U*V.';
    %R=R.';
    %optimaler Translationsvektor
    t=mue_x-R*mue_p;
    %Neue Punktwolke
    for i=1:lnght
        point=M_k(:,i);
        rPoint=R*point;
        M_k(:,i)=rPoint+t;
    end
    clf;
    paint_kartesian_image(X_k,'k');
    hold on;
    paint_kartesian_image(M_k,'m');
    hold on;
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