function b = slope_func(x, y)
y_low = mean(y) - 3*std(y);
y_high = mean(y) + 3*std(y);
I = y>y_low & y<y_high;
b = corr(x(I), y(I), 'Type', 'Pearson')*std(y, 1)/std(x, 1);