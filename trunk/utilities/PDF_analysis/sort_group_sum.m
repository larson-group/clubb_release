% $Id$
function [ sum_value ] = sort_group_sum( values )

% Function that sums all values in a given input vector of values.  This is
% done in a way to increase the accuracy of the sum.  First, all input values
% are sorted from smallest absolute value to largest absolute value, keeping
% their respective signs.  Summing from smallest magnitude to largest usually
% increases the accuracy of numerical addition.  Accuracy is further increased
% by grouping together values that are approximately the same order of
% magnitude.  Partial sums are taken in order to find the overall sum.  The sum
% of the values in an order-of-magnitude group is found for each
% order-of-magnitude group.  This produces a new vector of values.  The process
% repeats itself until there is only one value from each order-of-magnitude
% group.  Then, the sum over those values is taken, from smallest magnitude
% to largest, to find the overall sum of the values.

% Sort values by absolute value, keeping the sign of the values.
sort_mag_values = abs_val_signed_sort( values );
% Find the number of values.
num_values = max(size(values));

% Pass through any and all zero values (which will be found at the beginning of
% the vector sorted by absolute value), recording the index of the first value
% with an absolute value greater than 0.  The partial sum of all zero values
% is 0, so the zero values can be ignored.  In a situation where all values in
% the input vector are 0, the sum is 0, and the function can be exited.
val_index = 1;
while ( sort_mag_values(val_index) == 0.0 )
   if ( val_index == num_values )
      sum_value = 0.0;
      return
   else
      val_index = val_index + 1;
   end
end

% Find the lowest and highest order of magnitude of all the remaining values.
lowest_order_mag = round( log10( abs( sort_mag_values(val_index) ) ) );
highest_order_mag = round( log10( abs( sort_mag_values(num_values) ) ) );

% Group values by approximate order of magnitude and find the partial sum for
% each order-of-magnitude group.
order_mag_index = 1;
for order_mag = lowest_order_mag:1:highest_order_mag

   index = 0;
   while ( round( log10( abs( sort_mag_values(val_index) ) ) ) == order_mag )
      index = index + 1;
      group(index) = sort_mag_values(val_index);
      if ( val_index == num_values )
         break
      else
         val_index = val_index + 1;
      end
   end

   if ( index >= 1 )
      order_mag_sum(order_mag_index) = sum( group );
      order_mag_num(order_mag_index) = max(size(group));
      clear group
      order_mag_index = order_mag_index + 1;
   end

end

% If any of the order-of-magnitude groups contain multiple values, the sum of
% the values in that group might have the same order of magnitude as the sum of
% the values from another order-of-magnitude group.  Call this function
% recursively until the input vector is reduced to one value (or no values) from
% each order-of-magnitude group.
% Otherwise, sum the remaining values, from smallest absolute value to largest
% absolute value, keeping the sign of each value, to produce the overall sum.
if any( order_mag_num >= 2 )
   sum_value = sort_group_sum( order_mag_sum );
else
   sum_value = sort_sum( order_mag_sum );
end
