function [new_wf,new_wv]=WMVPCM(points,cluster_n,view_feat,b,alpha,beta,options)
 
%%INPUT VARIABLES
%points is a nxd matrix with n as the number of datapoints and d is the
%total number of dimensions
%cluster_n is the number of clusters
%view_feat is a 1xH vector consisting of the number of features in each h
%view, H is the total number of views
%label is the ground truth labels
%b is a seed number to be use for membership or cluster center
%initializations
 
% Change the following to set default options
default_options = [2;   % exponent for the partition matrix U
        100;    % max. number of iteration
        1e-4;   % min. amount of improvement for feature weights
        1e-4;   % min. amount of improvement for view weights
        1]; % info display during iteration 
 
if nargin == 2
    options = default_options;
else
    % If "options" is not fully specified, pad it with default values.
    if length(options) < 5;
        tmp = default_options;
        tmp(1:length(options)) = options;
        options = tmp;
    end
    % If some entries of "options" are nan's, replace them with defaults.
    nan_index = find(isnan(options)==1);
    options(nan_index) = default_options(nan_index);
    if options(1) <= 1;
        error('The exponent should be greater than 1!');
    end
end
 
m = options(1);         % Exponent for U
iter = options(2);      % Max. iteration
err_wf=options(3)         % Min. improvement for feature weights
err_wv = options(4);        % Min. improvement for view weights
display = options(5);       % Display info or not
 
rate=0;
thres=0.001; 
%Getting Ready:
points_n=size(points,1);
points_dim=size(points,2);
points_view=size(view_feat,2);
 
%Initialization of View weights
wv=ones(1,points_view)/points_view;
 
%Initialization of Feature Weights
%This can be adjusted depending how many views your dataset has
%For Example: for a two-view dataset, we define the following feature
%weight vetors in each view
wf1=ones(1,view_feat(1))/view_feat(1);
wf2=ones(1,view_feat(2))/view_feat(2);
wf=[wf1 wf2];
%you can add more depending on the number of views
 
%Initializing membership functions and cluster centers 
rand('seed',b)
%%%%%Initialization of mu
U = rand(cluster_n,points_n);
U_mat= rand(points_n,cluster_n);
u=U'; 
u_ori=u;
 
%%%%Initialize on cluster centers
rd=randi([1 points_n],1,cluster_n);
clust_cen=points(rd,:)+0.1;
orig_clust_cen=clust_cen;           



% calculate eta per view
 
    K=1;
    E=[];
    j=0;
    for h=1:points_view
        points_temp=points(:,j+1:j+view_feat(h));
        clust_cen_temp=clust_cen(:,j+1:j+view_feat(h));
       
        for k=1:cluster_n
        num=0;
        denom=0;
        for i = 1: size(points,1)
            num=num+u(i,k).^m.*(norm(points_temp(i)-clust_cen_temp(k))^2);
            denom=denom+u(i,k).^m;
        end
        eta(k)=K* (num/denom);
        
        end
     E=[E eta];
     j=j+view_feat(h);
    end
    
   

% 
% %%%Main Loop
for iter=1:5
time=[];
% while and(err_wf>thres,err_wv>thres)
%%% compute cluster center
tic;
rate=rate+1;
    new_clust_cen=[];
    for k=1:cluster_n
    temp4=zeros(1,points_dim);
    temp5=zeros(1,points_dim);
        for i=1:points_n
            temp4=temp4+u(i,k)^m*points(i,:);
            temp5=temp5+u(i,k)^m;
        end
        new_clust_cen=[new_clust_cen;temp4/temp5(k)];
    end
    new_clust_cen;

%     %update feature weight
    new_wf=[];
    j=0;
    W6=[];
    for h=1:points_view
        points_temp=points(:,j+1:j+view_feat(h));
        new_clust_cen_temp=new_clust_cen(:,j+1:j+view_feat(h));
        W5=[];
        for k=1:cluster_n
            W1=bsxfun(@minus,points_temp,new_clust_cen_temp(k,:)); 
            W2=W1.^2; 
            W3=bsxfun(@times,u(:,k).^m,W2); 
            W4=W3.*(wv(h).^alpha);
            W5=[W5;sum(W4,1)]; 
        end
%         W6=sum(W5,1).^(1./(beta-1));  %%for bbc-645
          W6=sum(W5,1).^(-1./(beta-1));
        W7=sum(W6,2);
        new_wf=[new_wf bsxfun(@rdivide,W6,W7)];
        j=j+view_feat(h);
    end
        new_wf;
%update view weight
    new_wv=[];
    j=0;
    V7=[];
    for h=1:points_view
        points_temp=points(:,j+1:j+view_feat(h));
        new_clust_cen_temp=new_clust_cen(:,j+1:j+view_feat(h));
        new_wf_temp=new_wf(:,j+1:j+view_feat(h));
        V5=[];
        for k=1:cluster_n
            V1=bsxfun(@minus,points_temp,new_clust_cen_temp(k,:)); 
            V2=V1.^2; 
            V3=bsxfun(@times,u(:,k).^m,V2);
            V4=bsxfun(@times,V3,new_wf_temp.^beta);
            V5=[V5;sum(V4,1)];
        end
        V6=sum(V5,1);
        V7=[V7 sum(V6,2)];
        j=j+view_feat(h);
    end
    V8=V7.^(-1./(alpha-1));
    V9=sum(V8,2);
    new_wv=[new_wv bsxfun(@rdivide,V8,V9)];

   
   new_u=[];
    u6=[];
        for k=1:cluster_n
        u5=[];
        j=0;
        y=0;
        for h=1:points_view
            points_temp=points(:,j+1:j+view_feat(h));
            new_clust_cen_temp=new_clust_cen(:,j+1:j+view_feat(h));
            new_wf_temp=new_wf(:,j+1:j+view_feat(h));     
           
        
            u1=bsxfun(@minus,points_temp,new_clust_cen_temp(k,:));
            u2=u1.^2;
            u3=bsxfun(@times,u2,new_wf_temp.^beta);
            u4=u3.*new_wv(h).^alpha;
            u5=[u5 (sum(u4,2))./E(1+cluster_n*y:cluster_n*(y+1))];  
            j=j+view_feat(h);
            y=y+1;
        end
        u6=[u6 1+sum(u5,2).^(1/(m-1))];
     end
     new_u=[new_u 1./u6];

    err_wf=norm(wf-new_wf);
    err_wv=norm(wv-new_wv);
    
    wf=new_wf;
    wv=new_wv;
    u=new_u;
    clust_cen=new_clust_cen;
    toc;
    time=[time toc];
end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clust=[];
T=u;
for i=1:points_n
    [num, idx]=max(T(i,:));
    clust=[clust;idx];
end
