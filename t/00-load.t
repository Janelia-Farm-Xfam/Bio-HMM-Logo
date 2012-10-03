#!perl -T

use Test::More tests => 1;

BEGIN {
    use_ok( 'Bio::HMM::Logo' ) || print "Bail out!\n";
}

diag( "Testing Bio::HMM::Logo $Bio::HMM::Logo::VERSION, Perl $], $^X" );
