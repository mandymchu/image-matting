function L = getLaplacian_CHU_V4(I1, I2, consts,vx)
  % V4: fix window size 5*5 in I1 and 3*3 in I2  

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
  r1 = 2;
  r2 = 1;
  
 for j = 1+r1:w-r1
   for i = 1+r1:h-r1    
      if (consts(i,j))    
        continue    % ignore any pixels under scribble regions.
      end  
      neb_size = 25;
      win_inds1 = indsM1(i-r1:i+r1,j-r1:j+r1); % get neighbour pixel indices of pixel(i,j).
      win_inds1 = win_inds1(:);      
      winI1 = I1(i-r1:i+r1,j-r1:j+r1,:);
      winI1 = reshape(winI1,neb_size,c); 
      winI = winI1;   
      win_inds = win_inds1;       
        
      dx = round(vx(i,j));   % look for the matched pixel in I2
      newj = j+dx;
      neb_size2 = 9;    
      if newj-r2 > 0 && newj+r2 <= w
         win_inds2 = indsM2(i-r2:i+r2,newj-r2:newj+r2);   % get neighbour pixel indices of pixel(i,j).
         win_inds2 = win_inds2(:);    
         winI2 = I2(i-r2:i+r2,newj-r2:newj+r2,:);  
         winI2 = reshape(winI2,neb_size2,c);   
         winI = [winI1; winI2]; 
         win_inds = [win_inds1; win_inds2];
         neb_size = neb_size+neb_size2;      
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
  L = spdiags(sumL(:),0,2*img_size,2*img_size)-L; % Sum of each row of Laplacian matrix is zero.;
 

return

