package CaseReader;

# Read a case file
#   The arg can be a relative or full path, or
#   it can be a file located somewhere in @INC.
sub readCase
{
	shift(@_);
	my $file = shift(@_);

	our $err = undef;
	{   
        # Put case data into a separate namespace
		package CASE;

		# Process the contents of the case file
		my $rc = do($file);

		# Check for errors
		if ($@) 
		{
			$err = "ERROR: Failure compiling '$file' - $@";
		} 
		elsif (! defined($rc))
		{
			$err = "ERROR: Failure reading '$file' - $!";
		}
		elsif (! $rc) 
		{
			$err = "ERROR: Failure processing '$file'";
		}
	}

	return ($err);
}

1;
