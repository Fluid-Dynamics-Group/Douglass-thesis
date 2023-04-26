function plotting_adjacency_new1(Aij_sin,nosc,type)
%% plotting
yyy = {'I','2','3','4','5','6',...
    '7','8'};
xx = 1:nosc+1;yy = 1:nosc+1;
A_temp = zeros(nosc+1,nosc+1);
A_temp(1:nosc,1:nosc) =  Aij_sin(1:nosc,1:nosc);
for i = 1:nosc+1
    A_temp(i,i) = NaN;
end
figure;h1 = subplot(2,2,1);
e2 = pcolor(xx,yy,A_temp);
if type == 1
    colormap(h1,flipud(hot));
elseif type == 2
    %colormap(h1,redbluecmap);
    colormap(h1,1*phasemap);
    %colormap(h1,hsv);
    %colormap(h1,hsv);
end
% if nosc<6
%     set(h1,'XTick',yy(1)+0.5:1:yy(end),'Xticklabel',...
%         yyy(1:1:nosc),'FontSize',14);
%     set(h1,'YTick',yy(1)+0.5:1:yy(end),'Yticklabel',...
%         yyy(1:1:nosc),'FontSize',14);
% else
%     set(h1,'XTick',yy(1)+0.5:2:yy(end),'Xticklabel',...
%         yyy(1:2:nosc),'FontSize',14);
%     set(h1,'YTick',yy(1)+0.5:2:yy(end),'Yticklabel',...
%         yyy(1:2:nosc),'FontSize',14);
% end
set(h1,'Ydir','reverse');
if type == 1
    %caxis([0 1]);
    hcb = colorbar;
    %set(hcb,'Fontname','symbol');
    set(hcb,'YTicklabel','');
    %plotTickLatex2D('axis',hcb)
elseif type == 2
    caxis([-pi pi]);
    hcb=colorbar;
    set(hcb,'YTick',[-pi,-pi/2,0,pi/2,pi]);
    %set(hcb,'YTicklabel',{'-p','-p/2','0','p/2','p'});
    set(hcb,'YTicklabel','');
    set(hcb, 'Fontname', 'symbol');
    set(e2,'AlphaData',1)
    %set(e2,'facealpha',0.999);
    %shading interp;
end
if type == 1
    %title('|A_{mn}|');
    caxis([0 5]);
elseif type == 2
    %title('\angle A_{mn}');
end
axis off;
set(hcb,'Fontsize',16);
axis square;
print('-depsc',['A_',num2str(type),'.eps']);
end

function c = redbluecmap(m)
if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b];
end
