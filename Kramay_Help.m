ontourf(t,f,abs(s).^2, 40,'LineColor', 'none')
close all
contourf(t,f,abs(s).^2, 40,'LineColor', 'none')
powerSpec = abs(s).^2;
doc zscore
powerSpec = zscore(powerSpec);
contourf(t,f,powerSpec, 40,'LineColor', 'none')
contourf(t,f(1:100),powerSpec(1:100, :), 40,'LineColor', 'none')
powerSpec = abs(s).^2;
powerSpec2 = log(powerSpec./mean(powerSpec(:,1:9),2));
contourf(t,f(1:100),powerSpec2(1:100, :), 40,'LineColor', 'none')
close all
contourf(t,f(1:100),powerSpec2(1:100, :), 40,'LineColor', 'none')