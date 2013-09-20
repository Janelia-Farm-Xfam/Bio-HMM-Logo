use strict;
use warnings;
use Test::More tests => 10;
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

my $logo = Bio::HMM::Logo->new({ hmmfile => $hmmfile });
isa_ok($logo, 'Bio::HMM::Logo');

# check the file value is returned
my $filename =  $logo->hmm_file();

is($hmmfile, $filename, "expect to see the same filename that we entered on object creation.");

#check the json output
my $logo_json = undef;
lives_ok { $logo_json = $logo->as_json()} q/Expect the as_png method to work if file is available/;

# uncomment this if you want to replace the json test file.
# Should do a diff of the two to make sure you aren't changing
# something unexpected.

#open my $new_json, '>',  $FindBin::Bin . '/data/test_new.json' or die;
#print $new_json $logo_json;
#close $new_json;

# compare the json string to the load file
my $logo_obj = decode_json $logo_json;
my $expected_obj = decode_json $expected;

cmp_deeply($logo_obj, $expected_obj, 'The logo json structure was created correctly');

my $png = undef;
lives_ok{$png = $logo->as_png()} q/Expect the as_png method to work if file is available/;

# test missing files
my $missing_logo = Bio::HMM::Logo->new({ hmmfile => '/does/not/exist' });
$logo_json = undef;
dies_ok {$logo_json = $missing_logo->as_json()} q/Expect the as_json method to die if hmm file is missing/;
like $@, qr|/does/not/exist does not exist on disk!|, q(Error message should report HMM file is missing);

# test corrupt files
my $corrupt_logo = Bio::HMM::Logo->new({ hmmfile => $FindBin::Bin . '/data/test.json' });
$logo_json = undef;
dies_ok {$logo_json = $corrupt_logo->as_json()} q/Expect the as_json method to die if hmm file is bad/;
like $@, qr|unable to open HMM file at|, q(Error message should report inability to open HMM file);


