use strict;
use warnings;
use Test::More tests =>4;
use Test::Deep;
use Test::Warn;
use Test::Exception;
use File::Slurp;
use Data::Printer;
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

cmp_deeply($logo_obj, $expected_obj, 'The logo json structure was created correctly');

my $png = undef;
$png = Bio::HMM::Logo::hmmToLogoPNG( $hmmfile );

# this is for testing score height calculations on an amino acid hmm
my $png_json = read_file( $FindBin::Bin . '/data/score_calc.json' );
my $png_input = decode_json $png_json;
my $aa_png = Bio::HMM::Logo::_build_png( $png_input, 'aa', 1, 'debug');
open my $png_file, '>', '/tmp/logo_text.png';
binmode $png_file;
print $png_file $aa_png;
close $png_file;

# test corrupt files
$logo_json = undef;
my $corrupted = $FindBin::Bin . '/data/test.json';
dies_ok {$logo_json =  Bio::HMM::Logo::hmmToLogoJson( $corrupted ) } q(Expect the hmmToLogoJson method to die if hmm file is bad);
like $@, qr|unable to open HMM file at|, q(Error message should report inability to open HMM file);
