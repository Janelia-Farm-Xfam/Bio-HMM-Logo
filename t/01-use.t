use strict;
use warnings;
use Test::More tests =>4;
use Test::Deep;
use Test::Warn;
use Test::Exception;
use File::Slurp;
use JSON;
use FindBin;

# load the module
BEGIN { use_ok 'Bio::HMM::Logo' }
# load the hmm and json file
my $hmmfile = $FindBin::Bin . '/data/dna.hmm';
my $expected = read_file($FindBin::Bin . '/data/test.json');
# create the json string from hmm
my $logo_json = undef;
$logo_json = Bio::HMM::Logo::hmmToLogoJson( $hmmfile );
# compare the json string to the load file
my $logo_obj = decode_json $logo_json;
my $expected_obj = decode_json $expected;

cmp_deeply($logo_obj, $expected_obj, 'Checking The logo json structure was created correctly. If this failed, check that you have the current version of Bio::HMM::Logo installed.');

# this is for testing score height calculations on an amino acid hmm
my $png_json = read_file( $FindBin::Bin . '/data/score_calc.json' );
my $png_input = decode_json $png_json;
my $aa_png = Bio::HMM::Logo::_build_png({
  data => $png_input,
  scaled => 1,
  debug => 'debug'
});
open my $build_file, '>', '/tmp/logo_build_png.png';
binmode $build_file;
print $build_file $aa_png;
close $build_file;

my $svg = Bio::HMM::Logo::hmmToLogoSVG( $hmmfile );
open my $svg_file, '>', '/tmp/logo_test.svg';
print $svg_file $svg;
close $svg_file;

# test corrupt files
$logo_json = undef;
my $corrupted = $FindBin::Bin . '/data/test.json';
dies_ok {$logo_json =  Bio::HMM::Logo::hmmToLogoJson( $corrupted ) } q(Expect the hmmToLogoJson method to die if hmm file is bad);
like $@, qr|unable to open HMM file at|, q(Error message should report inability to open HMM file);

my $png = undef;
$hmmfile = $FindBin::Bin . '/data/amino.hmm';
$png = Bio::HMM::Logo::hmmToLogoPNG( $hmmfile );
open my $png_file, '>', '/tmp/logo_test.png';
binmode $png_file;
print $png_file $png;
close $png_file;
