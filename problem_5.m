problem(xs,n_y1,n_y2);
function [val1,val2,val3,val4] = compute_meanline(xs,n_y1)
    Xmean=mean(xs);
    Ymean=mean(n_y1);
    sum=0;
    sumX=0;
    sum2=0;
    for i =1:length(xs)
        X=Xmean-xs(i);
        sumX=sumX+X;
        Y=Ymean-n_y1(i);
        sum=sum+(X*Y);
        sum2=sum2+(X^2);
    end
    val1=sum;
    val2=sum2;
    val3=Xmean;
    val4=Ymean;

end

function [line] = compute_meanlineSVD(xs,n_y1)

pts = [xs;n_y1;ones(size(xs))]';

[U S V]=svd(pts);
line = transpose(V(:,2));

end

function[xran1,xran2,yran1,yran2]  = ransac(xs,n_y1)
    margin = 13;
    inliers=0;
    xran1=0;
    xran2=0;
    yran1=0;
    yran2=0;
    
    for i=1:20000
        temp=0;
        x1 = (randi(numel(xs)));
        x2 = (randi(numel(xs)));
        
        y1 = (randi(numel(n_y1)));
        y2 = (randi(numel(n_y1)));
       if((x1~=x2)&&(y1~=y2))
           for j=1:length(xs)

               v1 = [x1 y1 0];
               v2 = [x2 y2 0];
               pt = [xs(j) n_y1(j) 0];
               d = point_to_line(pt, v1, v2);
               if(d<margin)
                   temp=temp+1;
               end
               
           end
        if(temp>inliers)
          inliers=temp;
          xran1=x1;
          xran2=x2;
          yran1=y1;
          yran2=y2;
        end
        end

    end

end

function d = point_to_line(pt, v1, v2)
      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a);
end

function [] = problem(xs,n_y1,n_y2)
intercept=0;
y=0;
y1=0;
% [sum,sum2,Xmean,Ymean]=compute_meanline(xs,n_y1);
%     m=double(sum/sum2);
%     intercept=Ymean-m*Xmean;
%     hold on;
%     plot(xs,n_y1,'o');
%     y=(m*xs)+intercept;
%     plot(xs,y);
%     hold off;
% 
% [sum,sum2,Xmean,Ymean]=compute_meanline(xs,n_y2);
%     m=double(sum/sum2);
%     intercept=Ymean-m*Xmean;
%     hold on;
%     plot(xs,n_y1,'o');
%     y=(m*xs)+intercept;
%     plot(xs,y);
%     hold off;


[line1]=compute_meanlineSVD(xs,n_y1);
disp('line parameters of the first coordinates are')
a=line1(1)
b=line1(2)
c=line1(3)

[line2]=compute_meanlineSVD(xs,n_y2);
disp('line parameters of the second coordinates are')
a=line2(1)
b=line2(2)
c=line2(3)

disp('line parameters of the ransac algorithm is')
[xran1,xran2,yran1,yran2]=ransac(xs,n_y2)

figure
subplot(1,3,1); 
hold on
plot(xs,n_y1,'g o');
y=(-line1(1)*xs-line1(3))/line1(2);
plot(xs,y)
title('first coordinates')
hold off

subplot(1,3,2); 
hold on
plot(xs,n_y1,'ro');
y1=(-line2(1)*xs-line2(3))/line2(2);
plot(xs,y1,'b')
title('second coordinates')
hold off
 
subplot(1,3,3);
m = (yran2-yran1)/(xran2-xran1);
intercept = yran1-m*xran1;
y=m*xs+intercept;
hold on   
plot(xs,n_y1,'o');
plot(xs,y,'g');
title('ransac')
hold off
end

