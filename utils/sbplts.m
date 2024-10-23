function [A,B] = sbplts(number)
% sbplts computes the distribution of the given number of subplots
%
% Date: 14/01/2021
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@mensa.cz

A = floor(sqrt(number));
while mod(number,A) ~= 0
    A = A-1;
end
B = number/A;

end

