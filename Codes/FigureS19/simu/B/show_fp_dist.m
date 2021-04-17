resultsfolder = 'results';
num_cycles = 300;
nbins_fp = 1000;
W = zeros(0,3);
R = zeros(0,3);

means = zeros(num_cycles,1);
vars = means;

for i = 1 : num_cycles
    pdt = [zeros(100,1),(1:100)'];
    s = load(strcat(resultsfolder,'/gen',num2str(i),'.mat'),'P','fp_manu','N_manu','L_manu');
    for j = 1 : 100
        pdt(j,1) = s.P{j}(end);
    end
    pdt = sortrows(pdt, 1);
    w_indx = pdt(end,2);
    r_indx = pdt(end-1,2);
    FL_w = [s.fp_manu{w_indx}, s.N_manu{w_indx} .* s.L_manu{w_indx}];
    FL_w(:,1) = round(nbins_fp * FL_w(:,1)) / nbins_fp;
    for j = 2 : length(FL_w)
        if isequal(FL_w(j,1),FL_w(j-1,1))
            FL_w(j,2) = FL_w(j,2) + FL_w(j-1,2);
            FL_w(j-1,2) = 0;
        end
    end
    %FL_w(or(FL_w(:,2) == 0,FL_w(:,1) == 0),:) = [];
    FL_w(FL_w(:,2) == 0,:) = [];
    FL_w(:,2) = FL_w(:,2) / sum(FL_w(:,2));
    %W = [W;[FL_w,ones(length(FL_w),1)*i + FL_w(:,2) * 0.5]];
    W = [W;[FL_w,ones(length(FL_w),1)*i]];
    means(i) = sum(FL_w(:,1) .* FL_w(:,2) / sum(FL_w(:,2)));
    vars(i) = sum(FL_w(:,1).^2 .* FL_w(:,2) / sum(FL_w(:,2))) ...
        - means(i)^2;
end
%%
f = figure;
tf1 = @(x) (log(x+.001)-log(.001)) .* 64 ./ (log(1.001)-log(.001) );
tf2 = @(x) x * 64;
subplot(2,2,[1,3])
scatter(W(:,3),W(:,1),14 + 3 * W(:,2),W(:,2),'filled')
%scatter3(W(:,3),W(:,1),W(:,2),14 + 3 * W(:,2),W(:,2),'filled')
set(gca,'Color','k')
subplot(2,2,2)
plot(means)
ylabel('mean')
subplot(2,2,4)
plot(vars)
ylabel('var')
colormap('hot')
savefig(f,'myfig')
