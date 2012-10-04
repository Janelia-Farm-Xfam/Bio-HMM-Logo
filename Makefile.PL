use strict;
use warnings;
use Inline::MakeMaker;
use Cwd;

my $orig_location = getcwd;

# need to get included src compiled here:
# cd to src dir
my $hmmer_src = __FILE__;

warn $hmmer_src, "\n";

$hmmer_src =~ s/Makefile.PL/src/;

chdir $hmmer_src;

system("./configure");
system("make");

chdir $orig_location;



WriteMakefile(
    NAME                => 'Bio::HMM::Logo',
    AUTHOR              => q{Jody Clements <clementsj@janelia.hhmi.org>},
    VERSION_FROM        => 'lib/Bio/HMM/Logo.pm',
    ABSTRACT_FROM       => 'lib/Bio/HMM/Logo.pm',
    ($ExtUtils::MakeMaker::VERSION >= 6.3002
      ? ('LICENSE'=> 'perl')
      : ()),
    PL_FILES            => {},
    PREREQ_PM => {
        'File::Slurp' => 0,
        'Test::More' => 0,
        'Test::Deep' => 0,
        'Test::Warn' => 0,
    },
    CONFIGURE_REQUIRES  =>  {
      'Inline::MakeMaker'     => 0.45,
      'ExtUtils::MakeMaker'   => 6.52,
    },
    dist                => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean               => { FILES => 'Bio-HMM-Logo-* _Inline Logo.inl' },
);
