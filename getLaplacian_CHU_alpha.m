function L = getLaplacian_CHU_alpha(I1, I2,alphaL,alphaR, consts,vx)
  % July  try to use adaptive window  vx is disparity map
  % 0717  get daptive window by alpha

  epsilon = 0.0000001;  
  [h,w,c] = size(I1);
  img_size = w*h; 
  indsM1 = reshape((1:img_size),h,w); % create indices for all pixels in image1, column first.
  indsM2 = indsM1+img_size;           % indices for all pixels in image2  
  tlen = sum(sum(1-consts(2:end-1,2:end-1)))*9^2;  % not accurately
  tlen = 2*tlen;
  row_inds = zeros(tlen ,1);  
  col_inds = zeros(tlen,1);
  vals = zeros(tlen,1);  
  
  len = 0;  
%   [~,w2] = size(vx);  % baby2 w=415, vx w2=413
%   w2 = min(w,w2);
  
   addpath adaptiveWindow
%    [~,pixelNum,edges,~] = getAdaptiveWindow_fV2(I1,I2,consts,vx); 

% ***********  add alpha to f***********
  f1 = I1;
  f2 = I2;
  if vx(1,1)>0  
        f1(:,:,4) = alphaR;
        f2(:,:,4) = alphaL;
  else
       f1(:,:,4) = alphaL;
       f2(:,:,4) = alphaR;
  end
  [~,pixelNum,edges,~] = getAdaptiveWindow_alpha(f1,f2,consts,vx);  
%   vx = newDisp;
  indx = 0;   % adaptive window's index   
   
 for j = 2:w-1
   for i = 2:h-1     %  same as adaptive window code
      if (consts(i,j))    
        continue % ignore any pixels under scribble regions.
      end  
      indx = indx+1;
      neb_size = pixelNum(indx); 
      top = edges(indx,1);
      right = edges(indx,2);
      bottom = edges(indx,3);
      left = edges(indx,4);

      win_inds1 = indsM1(top:bottom,left:right); % get neighbour pixel indices of pixel(i,j).
      win_inds1 = win_inds1(:);       
      winI1 = I1(top:bottom,left:right,:); % get neighbour pixel matrix from original image           
      winI1 = reshape(winI1,neb_size,c); % reshape neighbour pixel matrix into a neb_size*c matrix, each column represents a color channel.
      winI = winI1;
      win_inds = win_inds1;
      
%   if neb_size > 9
      dx = round(vx(i,j));   % look for the matched pixel in I2
%       dx = round(newD);
      shiftX = j+dx;
       if neb_size >= 25
          left  = shiftX-1;
          right = shiftX+1;
          top = i-1;
          bottom = i+1;
          neb_size2 = 9;
        else
         left = left+dx;
         right = right+dx;
         neb_size2 = neb_size;
       end    
      
       if left>0 && right<w+1  && top>0 && bottom <h+1 
          win_inds2 = indsM2(top:bottom,left:right);   % get neighbour pixel indices of pixel(i,j).
          win_inds2 = win_inds2(:);    
          winI2 = I2(top:bottom,left:right,:);  
          winI2 = reshape(winI2,neb_size2,c);   
          winI = [winI1; winI2]; 
          win_inds = [win_inds1; win_inds2];
          neb_size = neb_size+neb_size2;      
      end        
%   end
      
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
%    L = sparse(row_inds,col_inds,vals,img_size,img_size);
  sumL = sum(L,2);
  L = spdiags(sumL(:),0,2*img_size,2*img_size)-L; % Sum of each row of Laplacian matrix is zero.;
%   L = spdiags(sumL(:),0,img_size,img_size)-L; 

return

 % ********* adaptive window size ***************
%        oldD = vx(i, j);
%        newD = oldD;
%         % start with normalized window of size 3x3 centered at x, y
%         window = NewWindow(j, i);  % x=j, y=i
%          % initialize boolean flags for possible expansion
%         flags = [1 1 1 1];
%         pixelnum = 9;  % number of pixels in a window
% %         while abs(newD-oldD) > 1 && winSize < 30 % until disparity estimate converges or we reach max size   
%           while  pixelnum<30  && abs(newD-oldD)<3.5
%             % compute uncertainty
%             curUncert = uncertainty(stereoImg, window);
%               for k = 1:length(flags)                
%                  if flags(k)
%                     newWindow = NewWindow.expand(window, k);
%                     [newUncert, winReturned] = uncertainty(stereoImg, newWindow);                    
%                     if newUncert < curUncert
%                         uncerts(k) = newUncert;
%                         win(k) = winReturned;
%                     else
%                         uncerts(k) = Inf;
%                         win(k) = window;
%                         flags(k) = 0; % prohibit direction
%                     end
%                  end % if direction not prohibited
%              end % check all directions
%             m = min(uncerts);
%             if m == Inf
%                 break
%             end
%             direction = find(uncerts == m);
%             if(length(direction) > 1) % if find returns multiple mins
%                 break
%             end
%              % set new window
%             expandedWindow = win(direction);
%             currentUncert = uncerts(direction);
%             window = expandedWindow;
%              % compute the disparity increment    
%             increDisp = incrementDisp(stereoImg, expandedWindow, currentUncert);
%             stereoImg.DisparityMap(i,j) = stereoImg.DisparityMap(i,j) + increDisp;
%             newD = stereoImg.DisparityMap(i,j);            
%             pixelnum = (expandedWindow.edges(3) +1 - expandedWindow.edges(1)) ...
%                 *(expandedWindow.edges(2) - expandedWindow.edges(4)+1);            
%          end % while the disparity estimate hasn't converged
%           top = window.edges(1);
%           right = window.edges(2);
%           bottom = window.edges(3);
%           left = window.edges(4);   
%           neb_size = pixelnum;
          
%  ******** end of window size  *********
