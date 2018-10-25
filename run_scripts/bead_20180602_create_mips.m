mypath = '/Users/justin/Dropbox (MIT)/mac air/6_2_18_beads_and_fish/vertical_registered/';
p = [mypath '*.mat'];
d = dir(p);
fprintf('%d mat files found for\n%s\n\n',numel(d),p);
N = length(d);
for i=1:N
    fprintf('file %d of %d\n',i,N);
    save_2d_max_projections_from_file (d(i).folder, d(i).name)
    close all;
end