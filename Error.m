function [ce] = Error(ce)
e1=log(cell2mat(ce));
s1=size(e1);
x=1:1:s1(2);
figure
plot(x,e1)
hold on

title('Error AE for \Delta_3 in logarithmic scale')
xlabel('Iteration')
ylabel('Log (AE)')

end

