#!/usr/bin/perl
#CHANGE THE FOLLOWING PATH TO YOUR HOC DIRECTORY
$the_path="/home/hoc_browser/hoc_v2.2_tuner/";

-e $the_path || die ("$the_path not set right");

$outdir="hoc_ls"; #relative path to output directory

mkdir($outdir,0755) if (! -e $outdir);

@the_F_dirs=("src/", "src/LH/", "src/HOC_parameterization/", "src/GCSS/", "src/BUGSrad/", "src/coamps_micro/" );

@the_inc_dirs=(	"src/BUGSrad/" );

@the_prol_dirs=( "src/coamps_micro/" );

$the_files="*.[Ff]*";
foreach $dir (@the_F_dirs){
	@split_dir=split /\//,$dir;
	$title=pop @split_dir;
	$file_name="$outdir/".$title.".ls";
	$ls_opt=$the_path.$dir.$the_files;
	print "\npreparing $file_name\n\t from $ls_opt\n";
	system "ls $ls_opt >$file_name " || die ("crash 1\n");
}

$the_files="*.h";
foreach $dir (@the_inc_dirs){
	@split_dir=split /\//,$dir;
	$title=pop @split_dir;
	$file_name="$outdir/".$title.".ls";
	$ls_opt=$the_path.$dir.$the_files;
	print "\npreparing $file_name\n\t from $ls_opt\n";
	system "ls $ls_opt >>$file_name " || die ("crash 2\n");
}

$the_files="*.prol";
foreach $dir (@the_prol_dirs){
	@split_dir=split /\//,$dir;
	$title=pop @split_dir;
	$file_name="$outdir/".$title.".ls";
	$ls_opt=$the_path.$dir.$the_files;
	print "\npreparing $file_name\n\t from $ls_opt\n";
	system "ls $ls_opt >>$file_name " || die ("crash 2\n");
}

print "DONE\nNow try typing f90tohtml hoc.f2h\n";
print "(but you may later want to comment out troublesome file names in the .ls files)\n"

