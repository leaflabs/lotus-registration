function plot_rot (h, i, param)
figure(h);

a = size(param.rotvec);
x = 1:a(1);
yrad = param.rotvec(:,i);

subplot(1,3,i);
yyaxis left;
plot(x,yrad),'k';
hold on;
plot(1,param.rotvec(1,i),'ro');
%plot(a(1),param.bestr(i),'gs');
hold off;

xlabel('iteration');
ylabel(['rotation around dim ' num2str(i) ' (radians)']);
xlim([1 a(1)]);
y = ylim;

scale = 180 / pi;
%ydeg = yrad * scale;
yyaxis right;
%h2 = plot(x,ydeg,'k');
ny = y*scale;
ylim(ny);
%set(h2, 'Visible' ,'off');
end
