function out = plotAllPSFs (handle, bead_center, fpath)
subplot(3,1,2);
% for i=1:a(1)
% x axis = z pos = bead_center(:,3)
% y axis = fwhm in z = bead_center(:,4)
% n=[7 6 1:4];
% plot(bead_center(n,3),bead_center(n,4),'-o');
plot(bead_center(:,3),bead_center(:,4),'o');
% %
xlabel('z, i.e. depth (um)');
ylabel('full-width half-max in z (um)');
% plot i's
% tweak = 1.00;
% %for i=n
text(0, 70,['N = ' num2str(length(bead_center)) ' beads'],'FontSize',10);
% %    text(tweak*double(bead_center(i,3)),tweak*double(bead_center(i,4)),sprintf('[%d,%d]',zrange(i,1),zrange(i,2)),'FontSize',10);
% %end
% %hold off;
subplot(3,1,3);
% % x axis = y pos = height = bead_center(:,1)
% % y axis = fwhm in y = bead_center(:,2)
% plot(bead_center(n,3),bead_center(n,2),'-o');
plot(bead_center(:,3),bead_center(:,2),'o');
% %
xlabel('z, i.e. depth (um)');
ylabel('full-width half-max in y (um)');
% %
% % plot i's
% %for i=n
% %    text(tweak*double(bead_center(i,3)),tweak*double(bead_center(i,2)),sprintf('[%d,%d]',zrange(i,1),zrange(i,2)),'FontSize',10);
% %end
subplot(3,1,1);
% % x axis = y pos = height = bead_center(:,1)
% % y axis = fwhm in y = bead_center(:,2)
% plot(bead_center(n,3),bead_center(n,2),'-o');
plot(bead_center(:,4),bead_center(:,2),'o');
% %
xlabel('full-width half-max in z (um)');
ylabel('full-width half-max in y (um)');
%axis equal;

% figure(handle);
% str=sprintf('%ssummary.png',fpath);
% f=getframe(gcf);
% [X, map] = frame2im(f);
% imwrite(X, str);
%disp(sprintf('# Output file = %s',str));
end