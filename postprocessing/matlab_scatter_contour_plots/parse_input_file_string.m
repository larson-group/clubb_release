% $Id$
function [ filename_array ] ...
= parse_input_file_string( input_file_string, num_input_files )

% Loop over all characters in the input file string to find spaces, which
% are the filename delimiters.
string_length = zeros( num_input_files, 1 );
str_len_idx = 0;
delim_idx_prev = 0;
for char_idx = 1:1:size( input_file_string, 2 )-1
   if ( isspace( input_file_string(1,char_idx) ) )
      str_len_idx = str_len_idx + 1;
      string_length(str_len_idx) = char_idx - delim_idx_prev - 1;
      delim_idx_prev = char_idx;
   end
end % char_idx = 1:1:size( input_file_string, 2 )-1
string_length(num_input_files) ...
= size( input_file_string, 2 ) - delim_idx_prev;

% Setup size for filename_array by longest string.
char_max_len = blanks( max( string_length ) );

% Initialize filename_array.
filename_array(1,:) = char_max_len;

% Enter filenames into filename_array.
remaining_string = input_file_string;
for clubb_idx = 1:1:num_input_files-1
   [ filename_array(clubb_idx,1:string_length(clubb_idx)) ...
     remain ] = strtok( remaining_string, ' ' );
   remaining_string = remain;
end % clubb_idx = 1:1:num_input_files
filename_array(num_input_files,1:string_length(num_input_files)) ...
= strtok( remaining_string, ' ' );
