% 0721  get Adaptive Window V4 by comparing colors' variance
function [edges,pixelnum,Vari] = getAdaptiveWindow_V4(I1,consts)
% clear all
% I1 = double(imread('baby3_L.png'));
% trimap = im2double(imread('baby3_L_tri.png'));
% trimap = rgb2gray(trimap);
% consts = ~(trimap(:,:,1) > 0 & trimap(:,:,1) < max(trimap(:)));

if max(I1(:)) <= 1
    I1 = I1*255;
end
[h,w,c] = size(I1);
tlen = sum(sum(1-consts(2:end-1,2:end-1)));
edges = zeros(tlen,6); 
pixelnum = zeros(tlen,1);
% pixelNum = zeros(h,w);
inds = 0;
Var = zeros(1,4);
Vari = zeros(h,w);
% 
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

% for j = 200:201
%   for i = 9:11
r = 3;  % window's radius
for j = 2:w-1
  for i = 2:h-1    
      if (consts(i,j))    
        continue    % ignore any pixels under scribble regions.
      end   
      inds = inds+1;
  % ******* decide to choose window size********** 
   if i-r<1 || i+r>h || j-r<1 || j+r>w
%           pixelNum(i,j) = 9;
          pixelnum(inds) = 9;
          edges(inds,:) = [i-1,i+1,j-1,j+1,i,j];
   else 
         win1 = I1(i-1:i+1,j-1:j+1,:);
         win1 = reshape(win1,9,3);
         var1 = var(win1);
         var1 = mean(var1);     
         win2 = I1(i-r:i+r,j-r:j+r,:);   % r=3 7*7
         win2 = reshape(win2,(2*r+1)^2,3);
         var2 = var(win2);
         var2 = mean(var2); 
         Vari(i,j) = var2-var1;
     if abs(var2-var1) <= 880
%          pixelNum(i,j) = 25;   % choose big window 
         edges(inds,:) = [i-2,i+2,j-2,j+2,i,j];
         pixelnum(inds) = 25;
         continue
     end
         if abs(var2-var1)>880 && abs(var2-var1) <=1000
           category = 1;
           window = buildWin(i,j,category);  % start from 3*3 
           oldVar = var1;
         end
       if  var2-var1>1000
           win1 = I1(i-1:i+1,j-1:j+1,:);
           [FX,FY,~] = gradient(win1);
           FX = mean(abs(FX(:)));         
           FY = mean(abs(FY(:)));
           if FX > FY
               category = 2;   % 3*1
           else
               category = 3;   % 1*3
           end
%           category = class;
          window = buildWin(i,j,category);    % start from 3*1 or 1*3
          win1 = I1( window.edges(1):window.edges(2),window.edges(3):window.edges(4),:);
          win1 = reshape(win1,3,3);
          oldVar = var(win1);
          oldVar = mean(oldVar);
%           win1 = I1(i-1:i+1,j,:);   % 3*1
%           win1 = reshape(win1,3,3);
%           var1 = var(win1);
%           var1 = mean(var1); 
%           win2 = I1(i,j-1:j+1,:);  % 1*3
%           win2 = reshape(win2,3,3);
%           var2 = var(win2);
%           var2 = mean(var2); 
%          if var1 < var2
%             category = 2;
%             window = buildWin(i,j,category);  %  start from 3*1
%             oldVar = var1;
%          else
%             category = 3;
%             window = buildWin(i,j,category);  %  start from 1*3
%             oldVar = var2;
%          end
       end
           
     newVar = oldVar;
     if category == 1
        thr1 = 20;    % threshold
     else 
        thr1 = 9;
     end
          
     flags = [1 1 1 1];
     thr2 = 5;
     while newVar-oldVar < thr2 && window.pixelnum < thr1
       %          oldVar = newVar;
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
                Var(k) = Inf;
                flags(k) = 0;    % prohibit direction
             else
                win(k) = newWindow;       
                win1 = I1(top:bottom,left:right,:);
                win1 = reshape(win1,neb_size,c);
                t = var(win1);
                t = mean(t);   
                Var(k) = t;            
             end           
           end
       end
       newVar = min(Var);
       direction = find(Var == newVar);
       if(length(direction) > 1) % if find returns multiple mins
         break
       end
       if newVar-oldVar >= thr2
           break
       end
       window = win(direction);  % chosen expanded window
      
     end
%        pixelNum(i,j) = window.pixelnum;  
       pixelnum(inds) = window.pixelnum;
       edges(inds,:) = [window.edges,i,j];     
   end
          
  end
end
%  save win_edges edges pixelnum         
     

     
     
     
     
     

               
                 
                 

% 
% top = edges(1);
% bottom = edges(2);
% left = edges(3);
% right = edges(4);


% end



