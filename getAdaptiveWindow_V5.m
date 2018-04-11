function [edges,pixelnum] = getAdaptiveWindow_V5(I1,consts)
% 0723  get Adaptive Window V5 by comparing colors' standard deviation
% 0~255
% clear all
% I1 = double(imread('baby3_L.png'));
% trimap = im2double(imread('baby3_L_tri.png'));
% % trimap = rgb2gray(trimap);
% consts = ~(trimap(:,:,1) > 0 & trimap(:,:,1) < max(trimap(:)));

if max(I1(:)) <= 1
    I1 = I1*255;
end
[h,w,c] = size(I1);
r = 2;  % window's radius
pixelNum = zeros(h,w);
S1 = zeros(h,w);  % Standard deviation
S2 = zeros(h,w);
% deltaS = S1;
Std = zeros(1,4);
% [FX,FY,~] = gradient(I1);
% FX = abs(FX);
% FX = mean(FX(:));
% FY = abs(FY);
% FY = mean(FY(:));
% if FX > FY
%     class = 2;   % 3*1
% else 
%     class = 3;   % 1*3
% end
tlen = sum(sum(1-consts(2:end-1,2:end-1)));
pixelnum = zeros(tlen,1);
edges = zeros(tlen,4); 
inds = 0;

thr1 = 6;   % threshold
% for j = 200:201
%   for i = 9:11
for j = 2:w-1
  for i = 2:h-1    
      if (consts(i,j))    
        continue    % ignore any pixels under scribble regions.
      end   
      inds = inds+1;
       win1 = I1(i-1:i+1,j-1:j+1,:);
       win1 = reshape(win1,9,3);
%          s1 = var(win1);         
       s1 = std(win1);   % Standard deviation  0~255
       s1 = mean(s1);         
       if i-r<1 || i+r>h || j-r<1 || j+r>w
           window = buildWin(i,j,1); 
%            pixelNum(i,j) = 9;
           pixelnum(inds) = 9;
          edges(inds,:) = window.edges; 
           continue
       end
       if s1 > 11     %  choose small window 
           win1 = I1(i-1:i+1,j-1:j+1,:);
           [FX,FY,~] = gradient(win1);
           FX = mean(abs(FX(:)));         
           FY = mean(abs(FY(:)));
           if FX > FY
               category = 2;   % 3*1
           else
               category = 3;   % 1*3
           end
%          category = class;
         window = buildWin(i,j,category);     % start from 3*1 or 1*3
         win1 = I1( window.edges(1):window.edges(2),window.edges(3):window.edges(4),:);
         win1 = reshape(win1,3,3);
         s1 = std(win1);
         s1 = mean(s1);
       else  
         win2 = I1(i-r:i+r,j-r:j+r,:);  
         win2 = reshape(win2,(2*r+1)^2,3);
%          s2 = var(win2);
         s2 = std(win2);
         s2 = mean(s2); 
         if s2-s1 > thr1   % choose medium window
             category = 1;   
             window = buildWin(i,j,category);   % start from 3*3
     
         else
             pixelnum(inds) = 25;
             edges(inds,:) = [i-2,i+2,j-2,j+2];
             pixelNum(i,j) = 25;   % choose big window  5*5
             continue
         end         
       end
       
       oldS = s1;
       newS = oldS;
     if category == 1
        thr2 = 20;    % threshold
     else 
        thr2 = 9;
     end          
     flags = [1 1 1 1];
   
    while newS-oldS <= thr1 && window.pixelnum < thr2
         for k = 1:length(flags)                
           if flags(k)
              newWindow = buildWin.expand(window, k);
              top = newWindow.edges(1);
              bottom = newWindow.edges(2);
              left = newWindow.edges(3);
              right = newWindow.edges(4);
              neb_size = newWindow.pixelnum;
                      
             if top < 1 || left < 1 || bottom > h || right > w
                win(k) = window;
                Std(k) = Inf;
                flags(k) = 0;    % prohibit direction
             else
                win(k) = newWindow;       
                win1 = I1(top:bottom,left:right,:);
                win1 = reshape(win1,neb_size,c);
                t = std(win1);
                t = mean(t);   
                Std(k) = t;            
             end           
           end
         end
       newS = min(Std);
       direction = find(Std == newS);
       if(length(direction) > 1) % if find returns multiple mins
         break
       end
       if newS-oldS >= thr1
           break
       end
       window = win(direction);  % chosen expanded window
    end
%        pixelNum(i,j) = window.pixelnum;    
       edges(inds,:) = window.edges; 
       pixelnum(inds) = window.pixelnum;

         
  end
end
edges = edges(any(edges,2),:); 
pixelnum = pixelnum(1:inds);

% S1 = sqrt(S1);  % Standard deviation  0~255
% S2 = sqrt(S2);
% deltaS = S2-S1;
% con = S1~=0;
% s1 = S1(con);
% con = S2~=0;
% s2 = S2(con);
% 
% figure,h1 = histogram(s1,10);
% thr = h1.BinEdges(3);
% hold on
% h2 = histogram(s2,10);
% hold off
% con = deltaS~=0;
% delta = deltaS(con);
% figure,h3 = histogram(delta,10);
