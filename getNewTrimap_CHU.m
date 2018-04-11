function [newtri,alpha1,confidence1] = getNewTrimap_CHU(trimap1,alpha1,alpha2,confidence1,vx) 

% 0722 

[h,w] = size(trimap1);
consts = ~(trimap1(:,:,1) > 0 & trimap1(:,:,1) < max(trimap1(:))); 
consts = double(consts);   % unkown pixels ==0
consts(consts~=0) = 3;         % kown pixels ==3
consts(consts==0) = 5;         % unkown pixels ==5

alpha = alpha1>=0.9 | alpha1<=0.05;  
alpha = double(alpha);
alpha(alpha==1) = 5;           % alpha>=0.9 or <=0.1 pixels ==5
alpha(alpha==0) = 7;

both = consts-alpha;    % 0 == pixels whose alpha >=0.9 or <=0.1 in unkown part
both(both~=0) = 1;
newtri = trimap1;
% new = alpha1;

for j = 1:w
    for i = 1:h
        if both(i,j) ~= 0  
            continue
        end
         dx = round(vx(i,j));   % look for the matched pixel in I2
         newj = j+dx;
         if newj-1 > 0 && newj+1 <= w
            t = alpha2(i,newj-1:newj+1);
            t = mean(t);
            if abs(alpha1(i,j)-t)<=0.05 && confidence1(i,j)>0.7                 
                  if alpha1(i,j)>=0.9 
                        alpha1(i,j) = 1;
                        newtri(i,j) = 1;
%                   else  
%                         alpha1(i,j) = 0;
%                         newtri(i,j) = 0;
                  end
                  confidence1(i,j) = min(1,confidence1(i,j)+0.2);
            end
         end
        
        
        
    end
end



end