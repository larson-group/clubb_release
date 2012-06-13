% $Id$
function [ mag_sort_values ] = abs_val_signed_sort( values )

% Function that sorts values from smallest absolute value to largest absolute
% value, keeping the sign of each of the values the same.
%
% Examples:
% 1) -3, 5, 10, -8, 0.1 is sorted as: 0.1, -3, 5, -8, 10.
% 2) 3, 5, 10, 8, 0.1 is sorted as: 0.1, 3, 5, 8, 10.
% 3) -3, -5, -10, -8, -0.1 is sorted as: -0.1, -3, -5, -8, -10.
%
% Adding values from smallest magnitude to largest magnitude usually increases
% the accuracy of numerical floating point addition.

num_values = max(size(values));
sort_values = sort(values);

if all( sort_values >= 0.0 )

   % Only positive values are found.  They are already sorted from
   % smallest absolute value to largest absolute value, which is the
   % same as from smallest value to largest value.
   mag_sort_values = sort_values;

elseif all ( sort_values <= 0.0 )

   % Only negative values are found.  They need to be sorted from
   % smallest absolute value to largest absolute value.
   mag_sort_values = sort(values,'descend');

else % Both positive and negative values are found.

   index = 1;
   while ( index <= num_values )
      if ( sort_values(index) >= 0.0 )
         break
      else
         neg(index) = sort_values(index);
         index = index + 1;
      end
   end
   for index_2 = 1:1:num_values-index+1
      pos(index_2) = sort_values(index_2+index-1);
   end
   neg_mag_sort = sort(neg,'descend');

   iter = 1;
   iter_1 = 1;
   iter_2 = 1;
   while ( iter <= num_values )
      if ( pos(iter_1) <= abs(neg_mag_sort(iter_2)) )
         mag_sort_values(iter) = pos(iter_1);
         if ( iter_1 ~= max(size(pos)) )
            iter_1 = iter_1 + 1;
         else % iter_1 = max(size(pos))
            mag_sort_values(iter+1:num_values) ...
               = neg_mag_sort(iter_2:max(size(neg_mag_sort)));
            iter = num_values;
         end
      else % pos(iter_1) > abs(neg_mag_sort(iter_2))
         mag_sort_values(iter) = neg_mag_sort(iter_2);
         if ( iter_2 ~= max(size(neg_mag_sort)) )
            iter_2 = iter_2 + 1;
         else % iter_2 = max(size(neg_mag_sort))
            mag_sort_values(iter+1:num_values) = pos(iter_1:max(size(pos)));
            iter = num_values;
         end
      end
      iter = iter + 1;
   end

end
