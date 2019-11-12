clear all
clc

neur=400;

%% Connection vector & matrix
w=load('w400(integer).dat');
%w=abs(rand(neur,1));
%w=w/(abs(sum(w)));
%w=w/max(w);

k=linspace(0.0,1.0,neur+1);
k=k(2:end);

disp('Start printing...');

for i=1:neur-1,
    disp([num2str(w(i),'%2.5f') '*x(t-' num2str(k(i),'%2.4f') '*T)+\']);
end

disp([num2str(w(neur),'%2.5f') '*x(t-' num2str(k(neur),'%2.4f') '*T)']);
