use strict;
use warnings;
use Inline::MakeMaker;

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
    },
    CONFIGURE_REQUIRES  =>  {
      'Inline::MakeMaker'     => 0.45,
      'ExtUtils::MakeMaker'   => 6.52,
    },
    dist                => { COMPRESS => 'gzip -9f', SUFFIX => 'gz', },
    clean               => { FILES => 'Bio-HMM-Logo-* _Inline Logo.inl' },
);
