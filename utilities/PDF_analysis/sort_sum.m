% $Id$
function [ sum_value ] = sort_sum( values )

% Function that sums all values in a given input vector of values.  This is
% done in a way to increase the accuracy of the sum.  First, all input values
% are sorted from smallest absolute value to largest absolute value, with each
% value keeping its sign.  The sum of the values is taken, from smallest
% magnitude to largest, which usually increases the accuracy of numerical
% addition.

sort_mag_values = abs_val_signed_sort( values );
sum_value = sum( sort_mag_values );
