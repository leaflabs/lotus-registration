% data from /Users/justin/Dropbox\ \(MIT\)/mac\ air/for\ corban/run_1_\(distribution_assesment\)_on_181010_at_1831_overnight_run_oct23-24 

clear all;
close all;

centroid = [318.7500 251.2500 202];

trans_2324 = [
[304.8073 238.5933 210.6438]
[308.7473 227.1119 203.1876]
[304.8881 238.3458 210.2563]
[304.7791 237.6391 209.2856]
[315.6180 219.5217 59.6870]
[303.4449 236.7027 209.6014]
[316.0959 244.2201 230.2684]
[305.0300 238.5568 211.1501]
[316.7314 215.8726 140.7544]
[305.1232 238.2979 212.1617]
[304.6654 237.4298 208.6964]
[305.0612 238.5551 210.7158]
[299.1379 215.7045 188.8153]
[305.0073 238.6611 211.1495]
[304.5786 236.2974 210.4057]
[304.7372 238.9673 210.1826]
[304.8557 238.1311 209.1218]
[305.0807 238.7404 210.7275]
[305.0707 238.7153 209.1503]
[304.5074 238.1353 209.9338]
[304.8053 237.8411 210.7301]
[304.6903 239.0257 209.9817]
[304.8271 237.8925 211.1116]
[307.7236 240.8448 211.5781]
[305.3067 238.1112 211.7338]
];

rot_2324 = [
[-1.5174 0.0374 -0.0717]
[-1.5196 0.0011 -0.0211]
[-1.5134 0.0343 -0.0700]
[-1.5063 0.0340 -0.0673]
[-1.4759 0.0271 0.1183]
[-1.5152 0.0277 -0.0732]
[-1.5478 -0.2423 0.0296]
[-1.5095 0.0345 -0.0704]
[-1.4879 0.0609 0.0844]
[-1.5197 0.0308 -0.0684]
[-1.5203 0.0333 -0.0709]
[-1.5055 0.0299 -0.0661]
[-1.2776 -0.0559 -0.0254]
[-1.5171 0.0400 -0.0679]
[-1.5265 0.0374 -0.0696]
[-1.5226 0.0339 -0.0698]
[-1.5175 0.0353 -0.0706]
[-1.5295 0.0376 -0.0632]
[-1.5212 0.0297 -0.0753]
[-1.5147 0.0378 -0.0711]
[-1.5161 0.0343 -0.0686]
[-1.5040 0.0319 -0.0735]
[-1.5231 0.0357 -0.0693]
[-1.5009 0.0157 -0.0677]
[-1.5277 0.0300 -0.0713]
];

trans_2425 = [
[304.8510 238.3692 210.4369]
[309.9032 224.1347 202.1209]
[304.9211 238.2823 210.6172]
[304.8408 237.6001 209.2643]
[308.9418 222.7231 77.7877]
[302.9674 237.0090 209.2801]
[304.7826 237.0401 210.8376]
[304.9279 238.2586 211.0395]
[258.7663 216.6802 49.0685]
];

rot_2425 = [
[-1.5199 0.0364 -0.0705]
[-1.5281 -0.0095 -0.0051]
[-1.5125 0.0347 -0.0702]
[-1.5064 0.0338 -0.0674]
[-1.5173 0.0085 0.0544]
[-1.5183 0.0236 -0.0733]
[-1.5184 0.0355 -0.0694]
[-1.5090 0.0345 -0.0687]
[-1.9679 0.0964 0.1894]
];

trans_bad_super_seed = [
[305.0124 236.8703 211.2418]
[304.9666 238.9645 211.9400]
[303.2932 236.8549 209.3815]
[304.7822 237.1655 210.8914]
[304.7537 238.1568 211.4822]
[304.9746 238.4809 210.2254]
[304.6590 237.6751 210.0719]
];

rot_bad_super_seed = [
[-1.5182 0.0394 -0.0690]
[-1.5111 0.0343 -0.0699]
[-1.5171 0.0267 -0.0732]
[-1.5209 0.0349 -0.0695]
[-1.5251 0.0359 -0.0672]
[-1.5225 0.0390 -0.0711]
[-1.5191 0.0353 -0.0724]
];

a = size(trans_2324);
b = size(trans_2425);

trans_ground_truth = [-5 -10 15];
rot_ground_truth = [-87 2 -4]*pi/180;

trans_2324_zerod = trans_2324-centroid-trans_ground_truth;
rot_2324_zerod = rot_2324-rot_ground_truth;

trans_2425_zerod = trans_2425-centroid-trans_ground_truth;
rot_2425_zerod = rot_2425-rot_ground_truth;

trans_bad_super_seed_zerod = trans_bad_super_seed-centroid-trans_ground_truth;
rot_bad_super_seed_zerod = rot_bad_super_seed-rot_ground_truth;

trans_comp = [trans_2425_zerod;trans_2324_zerod(b(1)+1:end,:)];
rot_comp = [rot_2425_zerod;rot_2324_zerod(b(1)+1:end,:)];

N = 1:a(1); 
bad = [2 5 9 13];
kinda_bad = [6 7 24];
all_bad = sort([bad kinda_bad]);
all_good = N(~ismember(N,all_bad));
good = N(~ismember(N,bad));

trans_comp_super_seed = [trans_2324_zerod(all_good,:);trans_bad_super_seed_zerod];
rot_comp_super_seed = [rot_2324_zerod(all_good,:);rot_bad_super_seed_zerod];

f = figure;
subplot(1,2,1);
plot(1*ones(1,a(1)),trans_comp(:,1),'ko');
hold on;
plot(2*ones(1,a(1)),trans_comp(:,2),'ko');
plot(3*ones(1,a(1)),trans_comp(:,3),'ko');
plot(1*ones(1,numel(bad)),trans_comp(bad,1),'ro');
plot(2*ones(1,numel(bad)),trans_comp(bad,2),'ro');
plot(3*ones(1,numel(bad)),trans_comp(bad,3),'ro');
hold off;
xlabel('translation dims');
ylabel('error [um]');

subplot(1,2,2);
plot(4*ones(1,a(1)),rot_comp(:,1)*180/pi,'ko');
hold on;
plot(5*ones(1,a(1)),rot_comp(:,2)*180/pi,'ko');
plot(6*ones(1,a(1)),rot_comp(:,3)*180/pi,'ko');
plot(4*ones(1,numel(bad)),rot_comp(bad,1)*180/pi,'ro');
plot(5*ones(1,numel(bad)),rot_comp(bad,2)*180/pi,'ro');
plot(6*ones(1,numel(bad)),rot_comp(bad,3)*180/pi,'ro');
hold off;
xlabel('rotation dims');
ylabel('error [degrees]');

print(f,'figure1.png','-dpng');

f = figure;
subplot(1,2,1);
plot(1*ones(1,numel(good)),trans_comp(good,1),'ko');
hold on;
plot(2*ones(1,numel(good)),trans_comp(good,2),'ko');
plot(3*ones(1,numel(good)),trans_comp(good,3),'ko');
plot(1*ones(1,numel(kinda_bad)),trans_comp(kinda_bad,1),'bo');
plot(2*ones(1,numel(kinda_bad)),trans_comp(kinda_bad,2),'bo');
plot(3*ones(1,numel(kinda_bad)),trans_comp(kinda_bad,3),'bo');
hold off;
xlabel('translation dims');
ylabel('error [um]');

subplot(1,2,2);
plot(4*ones(1,numel(good)),rot_comp(good,1)*180/pi,'ko');
hold on;
plot(5*ones(1,numel(good)),rot_comp(good,2)*180/pi,'ko');
plot(6*ones(1,numel(good)),rot_comp(good,3)*180/pi,'ko');
plot(4*ones(1,numel(kinda_bad)),rot_comp(kinda_bad,1)*180/pi,'bo');
plot(5*ones(1,numel(kinda_bad)),rot_comp(kinda_bad,2)*180/pi,'bo');
plot(6*ones(1,numel(kinda_bad)),rot_comp(kinda_bad,3)*180/pi,'bo');
hold off;
xlabel('rotation dims');
ylabel('error [degrees]');

print(f,'figure2.png','-dpng');

f = figure;
subplot(1,2,1);
plot(1*ones(1,a(1)),trans_comp_super_seed(:,1),'ko');
hold on;
plot(2*ones(1,a(1)),trans_comp_super_seed(:,2),'ko');
plot(3*ones(1,a(1)),trans_comp_super_seed(:,3),'ko');
hold off;
xlabel('translation dims');
ylabel('error [um]');

subplot(1,2,2);
plot(4*ones(1,a(1)),rot_comp_super_seed(:,1)*180/pi,'ko');
hold on;
plot(5*ones(1,a(1)),rot_comp_super_seed(:,2)*180/pi,'ko');
plot(6*ones(1,a(1)),rot_comp_super_seed(:,3)*180/pi,'ko');
hold off;
xlabel('rotation dims');
ylabel('error [degrees]');

print(f,'figure3.png','-dpng');
