classdef buildWin
    % 0721 CHU  build a window
     properties
         i    %   the pixel's row index
         j    %  the pixel's column index
         edges   % array of edge coordinates  [top bottom left right]
         pixelnum   %  pixels number in this window
             
     end
     
     methods
         function obj = buildWin(i,j,category) 
             obj.i = i;
             obj.j = j;
             if category ==1  % 3*3 window
                 obj.edges = [i-1 i+1 j-1 j+1];
             else if category ==2  % 3*1 window
                     obj.edges = [i-1 i+1 j j];
                 end
                 if category == 3  % 1*3 window
                     obj.edges = [i i j-1 j+1];
                 end
             end
             h = obj.edges(2)-obj.edges(1)+1;
             w = obj.edges(4)-obj.edges(3)+1;
             obj.pixelnum = h*w;
         end
     end
         
     methods (Static)
         function newWin = expand(win, direction)
         %   DIRECTION: 1 - add row above           2 - add row below
               %        3 - add column on the left  4 - add column on the right
                if direction == 1 || direction == 3  
                     win.edges(direction) = win.edges(direction)-1;                                          
                else 
                     win.edges(direction) = win.edges(direction) + 1;
                end
                newWin = win;
                h = win.edges(2)-win.edges(1)+1;
                w = win.edges(4)-win.edges(3)+1;
                newWin.pixelnum = h*w;
         end
     end
          
end

