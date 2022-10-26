format compact

t=[0 1 2 3 4 5]; % Vector time
y=[1 2 4 8 12 13]; % Your signal
% The derivative can be approximated :
diff(y)./diff(t)

y = t.^2

diff(y)./diff(t)