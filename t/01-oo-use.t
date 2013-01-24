use strict;
use warnings;
use Test::More tests => 6;
use Test::Deep;
use Test::Warn;
use File::Slurp;
use Data::Printer;
use JSON;

# load the module
BEGIN { use_ok 'Bio::HMM::Logo' }
# load the hmm and json file
my $hmmfile = './t/data/test.hmm';
my $expected = read_file('./t/data/test.json');
# create the json string from hmm

my $logo = Bio::HMM::Logo->new({ hmmfile => $hmmfile });
isa_ok($logo, 'Bio::HMM::Logo');

# check the file value is returned
my $filename =  $logo->hmm_file();

is($hmmfile, $filename, "expect to see the same filename that we entered on object creation.");

#check the json output
my $logo_json = undef;
warning_like {$logo_json = $logo->as_json} qr/Setting character height to default method/;
# compare the json string to the load file
my $logo_obj = decode_json $logo_json;
my $expected_obj = decode_json $expected;

cmp_deeply($logo_obj, $expected_obj, 'The logo json structure was created correctly');

my $png = undef;
warning_like {$png = $logo->as_png()} qr/Setting character height to default method/;
