#!xPERLx -w

# creates the initial grid for solving using the laplacian thickness
# method

use strict;

use Getopt::Tabular;
use MNI::Spawn;
use MNI::Startup;
use MNI::FileUtilities qw(test_file check_output_dirs);

my ($stepValue);
my ($inputVolume, $grey, $white, $outputFile);
my ($tmpVolume);
my ($help, $usage);
my (@mincinfoOutput);

&Initialise;

# resample the input file if desired
if ($stepValue) {
    # find the current step sizes and dimlengths
    my @tmp;
    Spawn(['mincinfo', $inputVolume], stdout => \@mincinfoOutput);
    for (my $i = 5; $i <= 7; $i++) {
	@tmp = split /\s+/, $mincinfoOutput[$i];

	# get the dimension letter (i.e. x y z)
	my $dim = $tmp[1];
	$dim =~ s/(\w).+/$1/;

	# get the number of elements needed
	my $nelements = sprintf("%d", ($tmp[2] * $tmp[3]) / $stepValue);

	# now add the args
	AddDefaultArgs('mincresample', ["-${dim}nelements", $nelements,
		       "-${dim}step", $stepValue], 'pre');

    }
    Spawn(['mincresample', $inputVolume, $tmpVolume]);
}
else {
    $tmpVolume = $inputVolume;
}

# now generate the grid
Spawn(['create_laplacian_grid', $tmpVolume, $grey, $white, $outputFile]);

sub Initialise {

    # make sure necessary applications exist
    RegisterPrograms([qw(mincresample create_laplacian_grid mincinfo)]);

    &createInfoText;

    # defaults:
    
    my @leftOverArgs;
    my @argTbl = 
	(@DefaultArgs,
	 ["Programme specific options:", "section"],
	 ["-step", "float", 1, \$stepValue,
	  "Isostep size to use for the grid. Uses input sampling if unspecified."],
	 );

    GetOptions(\@argTbl, \@ARGV, \@leftOverArgs) or die "\n";

    $inputVolume = shift @leftOverArgs or die $usage;
    $grey = shift @leftOverArgs or die $usage;
    $white = shift @leftOverArgs or die $usage;
    $outputFile = shift @leftOverArgs or die $usage;
    
    # check if output exists and/or clobber was specified
    if (test_file('-e', $outputFile) && (! $Clobber)) {
	die "File $outputFile exists and -clobber was not specified\n";
    }

    # create the temporary volume
    $tmpVolume = "${TmpDir}/resampled_volume.mnc";
    check_output_dirs($TmpDir);

    AddDefaultArgs("mincresample", '-clobber');
    &self_announce if $Verbose;
}

sub createInfoText {
    $usage = <<USAGE;
Usage: $ProgramName [options] input.mnc grey_surface.obj ... 
           white_surface.obj output_grid.mnc
       $ProgramName -help for details
USAGE

$help = <<HELP;

No help text yet. Sorry.
HELP

Getopt::Tabular::SetHelp($help, $usage);

}


    
