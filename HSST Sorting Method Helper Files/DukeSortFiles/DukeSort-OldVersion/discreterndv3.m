function z = discreterndv3(p)
%
%   z = discreterndv3(p)
%
% Takes a m x n matrix and draws one discrete value for each column
% Output is row vector
%
%
%ALGORITHM
%=====================================================
%Grab a row for each column, based on how many values in that column are
%greater than a random variable being generated i.e. let's say we had some
%row with all 1s if our random variable was 0.5 that would give us an index
%of half the # of rows. If we only had a few non-zero points, for a given
%column we would only get indices that were in the range of those few
%non-zero points since they occupy the entirity of the non-zero part of the
%cumulative distribution.
%
%If we had 20 rows, of which only 2 were not-empty (for a given column)
%the output index for that column would either be 1 or 2, since all of the
%data is within those two samples

%It is not clear why a random integer sampling is not used
%NOTE: That would not be the same result, but it isn't clear
%why this result is desired

r            = rand(1,size(p,2)); % rand numbers are in the open interval (0,1)
p_cumulative = cumsum(p); 
r_n          = r.*p_cumulative(end,:);
z = sum(bsxfun(@gt,r_n,p_cumulative)) + 1;