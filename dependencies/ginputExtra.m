function [neurite] = ginputExtra(ROI)
% INPUT
% n: Number of points to plot
% booText: Boolean (default false) command to display point number in
% the plot.

% Author: Lasse Nørfeldt (Norfeldt) 
% Date: 2012-04-09

xg = []; yg = [];
stop=0;
i=0;
while (stop==0)
    imshow(ROI,[0,mean(ROI(:))*2]);
    hold on
    if i>0
        plot(xg(1),yg(1),'ro');
        for (int=2:i)
            plot(xg([int-1:int]),yg([int-1:int]),'r');
        end
    end
    i=i+1;
    hold off;
    [xi, yi] = ginput(1);
  
    xg = [xg xi]; yg = [yg yi];
    if (length(xg)<i)
        break;
    end
end
close all;
for i=1:length(xg)
   neurite(i,1)=xg(i);
   neurite(i,2)=yg(i);
end

end


