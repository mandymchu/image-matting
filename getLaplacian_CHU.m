function [L, op]=getLaplacian_CHU(IL,constsL,IR,constsR, epsilon,win_size)
  % 0531  try to use both samples of matched two images 
  % 0610  if cannot find matched pixel in image2, 
  % then winI = winI1,win_inds = win_inds1 not set winI2 to zeros
 
  load('matlab data\baby2_LtoR_vx')   %  load optical flow result to match pixels 
  
%  load('matlab data\baby2_L_smooth.mat');
% IL = IColor;   % original image
% load('matlab data\baby2_R_smooth.mat');
% IR = IColor;
% clear IColor pack 
% trimap = im2double(imread('baby2_L_tri.png'));
% constsL = ~(trimap(:,:,1) > 0 & trimap(:,:,1) < max(trimap(:)));
% trimap = im2double(imread('baby2_R_tri_ori.png'));
% constsR = ~(trimap(:,:,1) > 0 & trimap(:,:,1) < max(trimap(:)));

  if (~exist('epsilon','var'))
    epsilon=0.0000001;
  end  
  if (isempty(epsilon))
    epsilon=0.0000001;
  end
  if (~exist('win_size','var'))
    win_size=1;
  end     
  if (isempty(win_size))
    win_size=1;
  end     

%   win_size = 2;
   neb_size = (win_size*2+1)^2;   % neighbour size
  [h,w,c] = size(IR);
  img_size = w*h;
 
  indsM1 = reshape((1:img_size),h,w); % create indices for all pixels in image1, column first.
  indsM2 = indsM1+img_size;           % indices for all pixels in image2
    
if vx(1,1)>0
    I1 = IR;     I2 = IL;    % right image first, from right to left
    consts = imerode(constsR,ones(win_size*2+1)); 
%     consts = imerode(constsR,ones(3)); 
else
    I1 = IL;     I2 = IR;      % left image first
    consts = imerode(constsL,ones(win_size*2+1));    % erode scribbles
%     consts = imerode(constsL,ones(3)); 
end
   tlen1 = sum(sum(1-constsR(win_size+1:end-win_size,win_size+1:end-win_size)))*(neb_size^2); % calculate how many unknown pixels should be estimated.
   tlen2 = sum(sum(1-constsL(win_size+1:end-win_size,win_size+1:end-win_size)))*(neb_size^2); 
%     tlen1 = sum(sum(1-constsR(2:end-1,2:end-1)))*9;
%     tlen2 = sum(sum(1-constsL(2:end-1,2:end-1)))*9;
  tlen = tlen1+tlen2;
  row_inds = zeros(tlen ,1);
  col_inds = zeros(tlen,1);
  vals = zeros(tlen,1);  
  
  len = 0;  
  [~,w2] = size(vx);  % baby2 w=415, vx w2=413
  w2 = min(w,w2);
  
 for j = 1+win_size:w2-win_size  
     for i = win_size+1:h-win_size
      if (consts(i,j))    
        continue % ignore any pixels under scribble regions.
      end   
      neb_size = (win_size*2+1)^2;
%      win_inds1 = indsM1(i-winh:i+winh,j-winw:j+winw);
      win_inds1 = indsM1(i-win_size:i+win_size,j-win_size:j+win_size); % get neighbour pixel indices of pixel(i,j).
      win_inds1 = win_inds1(:); 
   %   winI1 = I1(i-winh:i+winh,j-winw:j+winw,:);
      winI1 = I1(i-win_size:i+win_size,j-win_size:j+win_size,:); % get neighbour pixel matrix from original image           
      winI1 = reshape(winI1,neb_size,c); % reshape neighbour pixel matrix into a neb_size*c matrix, each column represents a color channel.
          
      dx = round(vx(i,j));   % look for the matched pixel in I2
      %dy = round(vy(i,j));
%       m = i+dy;
      m = i;
      n = j+dx;
      if m-win_size>0 && n-win_size>0 && m+win_size<h+1 && n+win_size<w2+1
     %if m-winh>0 && n-winw>0 && m+winh<h+1 && n+winw<w2+1
         % win_inds2 = indsM2(m-winh:m+winh,n-winw:n+winw);
          win_inds2 = indsM2(m-win_size:m+win_size,n-win_size:n+win_size); % get neighbour pixel indices of pixel(i,j).
          win_inds2 = win_inds2(:);    
    %      winI2 = I2(m-winh:m+winh,n-winw:n+winw,:);
          winI2 = I2(m-win_size:m+win_size,n-win_size:n+win_size,:);  
          winI2 = reshape(winI2,neb_size,c);   
          winI = [winI1; winI2]; 
          win_inds = [win_inds1; win_inds2];
          neb_size = 2*neb_size;
      else 
          winI = winI1;
          win_inds = win_inds1;
      end          
      win_mu = mean(winI,1)';      
      win_var = inv(winI'*winI/neb_size-win_mu*win_mu' +epsilon/neb_size*eye(c)); % see equation 12.      
      winI = winI-repmat(win_mu',neb_size,1); % see equation 12.
      tvals = (1+winI*win_var*winI')/neb_size; % see equation 12.
      
      row_inds(1+len:neb_size^2+len) = reshape(repmat(win_inds,1,neb_size),...
                                             neb_size^2,1);
      col_inds(1+len:neb_size^2+len) = reshape(repmat(win_inds',neb_size,1),...
                                             neb_size^2,1);
      vals(1+len:neb_size^2+len) = tvals(:);           
      
      len = len+neb_size^2;
    end
  end  
    
  vals = vals(1:len); % cut unused space.
  row_inds = row_inds(1:len); % cut unused space.
  col_inds = col_inds(1:len); % cut unused space.
  L = sparse(row_inds,col_inds,vals,2*img_size,2*img_size); % Create the laplacian matrix, holding the sum of second part of equation 12. Note that the row,colomn pairs may have repeated position, so all elements of vals that have duplicate row and column'indices are added together, see the help.
  
  sumL = sum(L,2);
  L = spdiags(sumL(:),0,2*img_size,2*img_size)-L; % Sum of each row of Laplacian matrix is zero.
  op = vx(1,1);
 return



  
  

