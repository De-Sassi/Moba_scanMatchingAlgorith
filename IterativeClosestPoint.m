%**************************************************************************
%* Moba - Mobile Localization                                             *
%**************************************************************************
clc, clear all;
%Init
    %if first execution, read first data-cloud
    if exist('count_index', 'var') == 0
        count_index = 0;
        [filename_new] = get_next_filename(count_index);
        X=readtable(filename_new);                                %Fancy GUI needed ^^
    end
    %reads next data-cloud
    count_index = count_index+1;
    [filename_new] = get_next_filename(count_index);
    P=readtable(filename_new);
    
%     %translation vector
%     T=[1;1];
%     %rotation matrix
%     R=zeros(2);
    %minimize sum of the quadrats of the distance of the point pairs
    X=table2array(X);
    P=table2array(P);
%     matrix_height=X(90);                                                       %Wofür diese Variabeln?
%     matrix_width=max([X(1),X(180)]);                                           %Wofür diese Variabeln?

%First Round
    %Creats image of the data-cloud
    if count_index == 1
        %Execution if first loop-round
        [image]=create_image(X);
        [image]=paint_image(X, image);
        %converts polar data-cloud to a kartesic one as a matrix
        [X_k]=transvert_array(X);
        %Creates and keeps the origin-matrix as a copy of the first data-cloud
        [Z_k] = origin_matrix(X_k);
        %Objects position
        Pos=zeros(1500,2);

    end
    [P_k]=transvert_array(P);
    %calculates the optimal rotationmatrix and translationvector of the two
    %data-clouds
    [R, T] = optimal_R_and_T(X_k, P_k);
    [E] = qu_sum_distance(X_k, P_k,T,R);

%Transformation Loop
    %Adds new object-position to position-matrix
    Pos(count_index,:) = -T';                                                  %Noch nicht sicher ob korrect und zuverlässig
    %Data-Cloud Transformation Loop
    [M, d] = cloud_transformation (P_k, X_k, R, T, E);
    %Add new Datapoints to image
    [image]=paint_image(M);
    %Creates new data-files name
    [filename_new] = get_next_filename(count_index);
    %Increment count_index
    count_index = count_index+1;

%Results and Visualisations
    %Image Output
    imshow(image)



%Creates new image
function [image]=create_image(X)
if exist('image', 'var') == 0
    size=max(X);
    image=zeros(size*2);
end
end

%paint the data-cloud on a kartesic map
function [image]=paint_image(X, image)
size=max(X);
count=0;
for i=1:181
    alpha=i-1;
    d=X(i);
    if(d==0)
        %if there is no distance measured there is no measured point
        %so no point painting
        continue;
    end
    x_from_zero=cosd(alpha)*d;
    y_from_zero=sind(alpha)*d;
    %zero is in the middle of the image. so need to calculate the offset
    if(alpha>90)
        x_wOffset=round(x_from_zero+size);
    else
        x_wOffset=round(size-x_from_zero);
    end
     y_wOffset=round(size-y_from_zero);
%Für Verfeinerung der Image-Dartesstung.
%Evt Intensität steigern, zum optisch verstäken
%     if (image(x_wOffset,y_wOffset)==255)
%         .....
%     end
    image(x_wOffset,y_wOffset)=255;
    count=count+1;
end
end

%transverts polar to kartesic koordinates and creates a new matrix
function [kart_matrix] = transvert_array(X)
size=max(X);
count=0;
kart_matrix=zeros(181,2);
for i=1:181
    alpha=i-1;
    d=X(i);
    if(d==0)
        %if there is no distance measured there is no measured point
        %so no point painting
        continue;
    end
    %transformation polar to kartesic
    x_from_zero=cosd(alpha)*d;
    y_from_zero=sind(alpha)*d;
    %zero is in the middle of the image. so need to calculate the offset
    if(alpha>90)
        x_wOffset=round(x_from_zero+size);
    else
        x_wOffset=round(size-x_from_zero);
    end
    y_wOffset=round(size-y_from_zero);
    %fills kartesic matrix
    kart_matrix(i,1)=x_wOffset;
    kart_matrix(i,2)=y_wOffset;
    count=count+1;
end
end

%Origin Data-Cloud
function [Z_k] = origin_matrix(X_k)
Z_k=zeros(181,2);
Z_k=X_k;
end
%mean value of the datafiles
function [R, T] = optimal_R_and_T(X_k, P_k)
mue_x=sum(X_k)/length(X_k);
mue_p=sum(P_k)/length(P_k);
W=zeros(2);
for j=1:181
    W=W+(X_k(j)-mue_x)*transpose(P_k(j)-mue_p);
end
%Singular value decomposition of matrix W
[U,E,V] = svd(W);
%Optimal rotationsmatrix R
R = U * transpose(V);
%Optimaler Translationsvektor T
T = transpose(mue_x) - R * transpose(mue_p);
%Get filename of next datafile
end

%Quadratsumme der Abstände minimieren
function [E] = qu_sum_distance(X_k, P_k,T,R)
E=0;
    for j=1:181
%       M = X_k(j,:)'-R*P_k(j,:)'-T             
% Da bereits Betragsvektor, Negativ-Check üebrflüssig
%       if M(1,1)<=0
%             M(1,1)=-M(1,1);
%           end
%       if M(2,1)<=0
%           M(2,1)=-M(2,1);
%           end
        E=E+norm((X_k(j,:)'-R*P_k(j,:)'-T))^2;
    end
E = E/length(X_k);
end

%Transforms data-clouds to a new data_cloud via Rotationmatrix and translation vector
function [M] = new_data_cloud(P_k, R, T)
M=zeros(181,2);
    for j=1:181
        M(j,:) = P_k(j,:) * R + T';                 
    end
end
 
%Data-Cloud Transformation Loop:
%Repeats point-cloud transformation as long as average point-distance (d) 
%is bigger than minimal square-sum distance (E).
function [M, d] = cloud_transformation(P_k, X_k, R, T, E)
d=norm(T);
    while (d >= E)
        [M] = new_data_cloud(P_k, R, T);
        [R, T] = optimal_R_and_T(X_k, M);
        for j=1:4
            [M] = new_data_cloud(P_k, R, T);
            %Inproves data-cloud with random little changes
            M = M + 0.01*randn(181,2);
            [R, T] = optimal_R_and_T(X_k, M);
        end
        d=norm(T);
    end
end

%Creates the filename of the next upcoming datafile
function [filename_new] = get_next_filename(count_index)
   
   %incrementing index of the next datafile
   i = count_index;
   %String splits
   a = 'ScanFilesNew\Scan';
   b = '0';
   c = num2str(i);
   d = '.dat';

   %adding zeros depending on the index-length
     if (i>9999)
        b = '0';
     elseif (i>999)
        b = '00';
     elseif (i>99)
        b = '000';
     elseif (i>9)
        b = '0000';
     else
        b = '00000';
     end
    %Unite editted string splits
    filename_new = [a,b,c,d]  
end
