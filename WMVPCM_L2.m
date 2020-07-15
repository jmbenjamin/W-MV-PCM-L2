% revised program for W-MV-PCM-L2

% %%%Cleaning the work space
%     clear all
%     close all
%     clc
%     
% % % %%%Loading the data
%  
% load mvc_ex5_2clusts2viewswn.mat
% cluster_n=2;
% view_feat=[3 3];
% points_n=size(points,1);
% points_dim=size(points,2);
% points_view=size(view_feat,2);
% label=label1;
% %%% Initialization of view weight
% wv=ones(1,points_view)/points_view;
% wf1=ones(1,view_feat(1))/view_feat(1);
% wf2=ones(1,view_feat(2))/view_feat(2);
% wf=[wf1 wf2];
% % % % 
% % % % % 
% % % % %%%Initial Inputs
% % m=2;
% % err_clust_cen=10;
% % rate=0;
% % err_wf=10;
% % err_wv=10;
% % thres=0.001;
% % max_iter=100;
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time=[];
% % 
% % %%%%Alpha and Beta
% alpha=[];
% alpha=[alpha cluster_n./view_feat];
% beta=[];
% beta=[beta var(points,1)];
% 
% 
% b=27;
% rand('seed',b)
% % %%%%%Initialization of mu
% U = rand(cluster_n,points_n);
% U_mat= rand(points_n,cluster_n);
% u=U'; 
% u_ori=u;
% 
% %%%%Initialize on cluster centers
% rd=randi([1 points_n],1,cluster_n);
% clust_cen=points(rd,:)+0.1;
% orig_clust_cen=clust_cen;                            
%  
% 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%calculate eta per view
%  
%     K=1;
%     E=[];
%     j=0;
%     for h=1:points_view
%         points_temp=points(:,j+1:j+view_feat(h));
%         clust_cen_temp=clust_cen(:,j+1:j+view_feat(h));
%        
%         for k=1:cluster_n
%         num=0;
%         denom=0;
%         for i = 1: size(points,1)
%             num=num+u(i,k).^m.*(norm(points_temp(i)-clust_cen_temp(k))^2);
%             denom=denom+u(i,k).^m;
%         end
%         eta(k)=K.* (num/denom);
%         
%         end
%      E=[E eta];
%      j=j+view_feat(h);
%     end
% 
% %%%%Main Loop
%   while and(err_wf>thres,err_wv>thres)
% 
%       tic;
%       rate=rate+1;
%       
%     %%%update cluster center
%     new_clust_cen=[];
%     for k=1:cluster_n
%     temp4=zeros(1,points_dim);
%     temp5=zeros(1,points_dim);
%         for i=1:points_n
%            temp4=temp4+u(i,k)^m*points(i,:);
%           temp5=temp5+u(i,k)^m;
%         end
%         new_clust_cen=[new_clust_cen;temp4./temp5(k)];
%     end
%     new_clust_cen;
%     
%     %update view weight
%     new_wv=[];
%     j=0;
%     V7=[];
%     V5=[];
%     V51=[];
%     for h=1:points_view
%         points_temp=points(:,j+1:j+view_feat(h));
%         new_clust_cen_temp=new_clust_cen(:,j+1:j+view_feat(h));
%         wf_temp=wf(:,j+1:j+view_feat(h));
%         V5=[];
%         
%         for k=1:cluster_n
%             V1=bsxfun(@minus,points_temp,new_clust_cen_temp(k,:)); 
%             V2=V1.^2; 
%             V3=bsxfun(@times,u(:,k).^m,V2);
%             V4=bsxfun(@times,V3,wf_temp.^2);
%             V5=[V5; sum(V4,1)];
%             V6=bsxfun(@plus,V5,alpha(h)); 
%         end
%             
%             V7=[V7 sum(V6,2)];
%            j=j+view_feat(h);
%     end
%         
%         V8=sum(V7,2);
%         V9=bsxfun(@rdivide,V7,V8);
%         new_wv=[new_wv V9(1, 1:size(V9,2))];
%         sum_new_wv=sum(new_wv,2);
%     
%     %update feature weight
%     new_wf=[];
%     j=0;
%     W6=[];
%     for h=1:points_view
%         points_temp=points(:,j+1:j+view_feat(h));
%         new_clust_cen_temp=new_clust_cen(:,j+1:j+view_feat(h));
%         beta_temp=beta(:,j+1:j+view_feat(h));
%         W5=[];
%         for k=1:cluster_n
%             W1=bsxfun(@minus,points_temp,new_clust_cen_temp(k,:)); 
%             W2=W1.^2; 
%             W3=bsxfun(@times,u(:,k).^m,W2); 
%             W4=W3*new_wv(h).^2;
%             W5=[W5;sum(W4,1)]; 
%         end
% 
%         W6=sum(W5,1)+beta_temp;
%         W7=sum(W6,2);
%         new_wf=[new_wf bsxfun(@rdivide,W6,W7)];
%         j=j+view_feat(h);
%     end
%        sum_wf=sum(new_wf,2);
%       
%       
%       new_u=[];
%       u6=[];
%         for k=1:cluster_n
%         u5=[];
%         j=0; 
%         y=0;
%         for h=1:points_view
%             points_temp=points(:,j+1:j+view_feat(h));
%             new_clust_cen_temp=new_clust_cen(:,j+1:j+view_feat(h));
%             new_wf_temp=new_wf(:,j+1:j+view_feat(h));     
%             beta_temp=beta(:,j+1:j+view_feat(h));
%         
%             u1=bsxfun(@minus,points_temp,new_clust_cen_temp(k,:));
%             u2=u1.^2;
%             u3=bsxfun(@times,u2,new_wf_temp.^2);
%             u4=u3.*new_wv(h).^2;
%             u5=[u5 ((sum(u4,2))./E(1+cluster_n*y:cluster_n*(y+1)))];  
%             j=j+view_feat(h);
%             y=y+1;
%         end
%         u6=[u6 1+(sum(u5,2).^(1/(m-1)))];
%      end
%      new_u=[new_u 1./u6];
%      u=new_u;
%      
%         err_wf=norm(wf-new_wf);
%         err_wv=norm(wv-new_wv);
%   
%     wf=new_wf;
%     wv=new_wv;
%     u=new_u;
%     toc;
%     time=[time toc];
%       end
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     clust=[];
% T=u;
% for i=1:points_n
%     [num, idx]=max(T(i,:));
%     clust=[clust;idx];
% end

function [new_wf,new_wv]=WMVPCM_L2(points,cluster_n,view_feat,b,options)
 
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
 
%Computing for Alpha and Beta
alpha=[];
alpha=[alpha cluster_n./view_feat];
beta=[];
beta=[beta var(points,1)];
 
%calculate eta per view
 
    K=1;         %default
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
        eta(k)=K.* (num/denom);
        
        end
     E=[E eta];
     j=j+view_feat(h);
    end
    
    time=[];
    %%%%Main Loop
%   while and(err_wf>thres,err_wv>thres)
 for iter=1:20
      tic;
      rate=rate+1;
      
    %%%update cluster center
    new_clust_cen=[];
    for k=1:cluster_n
    temp4=zeros(1,points_dim);
    temp5=zeros(1,points_dim);
        for i=1:points_n
           temp4=temp4+u(i,k)^m*points(i,:);
          temp5=temp5+u(i,k)^m;
        end
        new_clust_cen=[new_clust_cen;temp4./temp5(k)];
    end
    new_clust_cen;
    
    %update view weight
    %update view weight
    new_wv=[];
    j=0;
    V7=[];
    V5=[];
    V51=[];
    for h=1:points_view
        points_temp=points(:,j+1:j+view_feat(h));
        new_clust_cen_temp=new_clust_cen(:,j+1:j+view_feat(h));
        wf_temp=wf(:,j+1:j+view_feat(h));
        V5=[];
        
        for k=1:cluster_n
            V1=bsxfun(@minus,points_temp,new_clust_cen_temp(k,:)); 
            V2=V1.^2; 
            V3=bsxfun(@times,u(:,k).^m,V2);
            V4=bsxfun(@times,V3,wf_temp.^2);
            V5=[V5; sum(V4,1)];
            V6=bsxfun(@plus,V5,alpha(h)); 
        end
            
            V7=[V7 sum(V6,2)];
           j=j+view_feat(h);
    end
        
        V8=sum(V7,2);
        V9=bsxfun(@rdivide,V7,V8);
        new_wv=[new_wv V9(1, 1:size(V9,2))];
        sum_new_wv=sum(new_wv,2);
    
    
    %update feature weight
    new_wf=[];
    j=0;
    W6=[];
    for h=1:points_view
        points_temp=points(:,j+1:j+view_feat(h));
        new_clust_cen_temp=new_clust_cen(:,j+1:j+view_feat(h));
        beta_temp=beta(:,j+1:j+view_feat(h));
        W5=[];
        for k=1:cluster_n
            W1=bsxfun(@minus,points_temp,new_clust_cen_temp(k,:)); 
            W2=W1.^2; 
            W3=bsxfun(@times,u(:,k).^m,W2); 
            W4=W3*new_wv(h).^2;
            W5=[W5;sum(W4,1)]; 
        end
 
        W6=sum(W5,1)+beta_temp;
        W7=sum(W6,2);
        new_wf=[new_wf bsxfun(@rdivide,W6,W7)];
        j=j+view_feat(h);
    end
       sum_wf=sum(new_wf,2);
      
      
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
            beta_temp=beta(:,j+1:j+view_feat(h));
        
            u1=bsxfun(@minus,points_temp,new_clust_cen_temp(k,:));
            u2=u1.^2;
            u3=bsxfun(@times,u2,new_wf_temp.^2);
            u4=u3.*new_wv(h).^2;
            u5=[u5 ((sum(u4,2))./E(1+cluster_n*y:cluster_n*(y+1)))];  
            j=j+view_feat(h);
            y=y+1;
        end
        u6=[u6 1+(sum(u5,2).^(1/(m-1)))];
     end
     new_u=[new_u 1./u6];
     u=new_u;
     
     %Errors
     err_wf=norm(wf-new_wf);
     err_wv=norm(wv-new_wv);
 
    wf=new_wf;
    wv=new_wv;
    u=new_u;
    toc;
    time=[time toc];
   
  end
     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %New clust lables
    clust=[];
T=u;
for i=1:points_n
    [num, idx]=max(T(i,:));
    clust=[clust;idx];
end
end

 
