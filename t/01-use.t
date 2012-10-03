use strict;
use warnings;
use Test::More;
use File::Slurp;
use Data::Printer;

# load the module
BEGIN { use_ok 'Bio::HMM::Logo' }
# load the hmm and json file
my $hmmfile = './t/data/test.hmm';
# create the json string from hmm
my $json = Bio::HMM::Logo::hmmToLogoJson( $hmmfile );
# compare the json string to the load file
p $json;

done_testing();
