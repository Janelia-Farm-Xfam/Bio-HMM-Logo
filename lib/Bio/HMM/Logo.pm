package Bio::HMM::Logo;

use strict;
use warnings;
use JSON;
use File::Spec;


=head1 NAME

Bio::HMM::Logo - The great new Bio::HMM::Logo!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

my $src_file   = undef;
my $typemaps   = undef;
my $src_dir    = undef;

BEGIN {
  $src_file = __FILE__;
  $src_file =~ s/\.pm/\.c/;

  my $file = __FILE__;
  ($src_dir) = $file =~ /^(.*)\/blib/;
  $src_dir = File::Spec->catfile($src_dir, 'src');
  $typemaps = __FILE__;
  $typemaps =~ s/\.pm/\.typemap/;
}

use Inline
  C        => "$src_file",
  VERSION  => '0.01',
  ENABLE   => 'AUTOWRAP',
  INC      => "-I$src_dir -I$src_dir/easel",
  LIBS     => "-L$src_dir/easel -L$src_dir -lhmmer -leasel",
  TYPEMAPS => $typemaps,
  NAME     => 'Bio::HMM::Logo';

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use Bio::HMM::Logo;

    my $foo = Bio::HMM::Logo->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=cut

=head2 hmmToLogoJson

=cut

sub hmmToLogoJson {
  my ( $hmmfile, $method ) = @_;

  if ( !$method || $method !~ /^(emission|posscore|score)$/ ) {
    warn "Setting character height to default method [emission]\n";
    $method = 'emission';
  }

  unless ( -e $hmmfile ) {
    die "$hmmfile does not exist on disk!\n";
  }
  my $hmm = inline_read_hmm($hmmfile);
  my $abc = inline_get_abc($hmm);

  my $alph     = inline_get_alphabet_string($abc);
  my @alph_arr = split( //, $alph );

  my $max_height;
  my $height_arr_ref;
  if ( $method eq "emission" ) {
    $max_height     = inline_hmmlogo_maxHeight($abc);
    $height_arr_ref = inline_get_emission_heights($hmm);
  }
  elsif ( $method eq "posscore" ) {
    $max_height     = inline_hmmlogo_maxHeight($abc);
    $height_arr_ref = inline_get_posscore_heights($hmm);
  }
  elsif ( $method eq "score" ) {
    $height_arr_ref = inline_get_score_heights($hmm);
  }

  foreach my $row (@$height_arr_ref) {
    my %char_heights;
    for my $i ( 0 .. $#{$row} ) {
      $char_heights{ $alph_arr[$i] } = 0 + sprintf( "%.3f", ${$row}[$i] );
    }

    #sort by height
    my @sorted_keys =
      sort { $char_heights{$a} <=> $char_heights{$b} } keys %char_heights;

    for my $i ( 0 .. $#{$row} ) {
      my $key = $sorted_keys[$i];
      ${$row}[$i] = "$key:$char_heights{$key}";
    }
  }

  my $insertP = inline_get_insertP($hmm);
  foreach my $v (@$insertP) {
    $v = 0 + sprintf( "%.0f", 100 * $v );
  }

  my $insert_len = inline_get_insertLengths($hmm);
  foreach my $v (@$insert_len) {
    $v = 0 + sprintf( "%.0f", $v );
  }

  my $deleteP = inline_get_deleteP($hmm);
  foreach my $v (@$deleteP) {
    $v = 0 + sprintf( "%.0f", 100 * ( 1 - $v ) );
  }


  my $mm = inline_get_MM_array($hmm);


  my $height_data_hashref = {
    max_height      => $max_height,
    height_arr      => $height_arr_ref,
    insert_probs    => $insertP,
    insert_lengths  => $insert_len,
    occupancy_probs => $deleteP,
    mmline          => $mm,
  };

  #This destory was causing issues (no return string) when it occured just
  #before the return.  I think that perl would actually look after us and
  #garbage collect anyway, but this seems to solve the problem. I can not
  #explain why it was destroying the JSON object, possibly a memory violation.

 # inline_destroy_abc($abc);
 # inline_destroy_hmm($hmm);

  my $json             = JSON->new->allow_nonref;
  my $height_data_json = $json->encode($height_data_hashref);

  return $height_data_json;
}

=head2 dl_load_flags
=head2 inline_destroy_abc
=head2 inline_destroy_hmm
=head2 inline_get_MM_array
=head2 inline_get_abc
=head2 inline_get_alphabet_string
=head2 inline_get_deleteP
=head2 inline_get_emission_heights
=head2 inline_get_insertLengths
=head2 inline_get_insertP
=head2 inline_get_posscore_heights
=head2 inline_get_score_heights
=head2 inline_hmmlogo_maxHeight
=head2 inline_read_hmm

=head1 AUTHORS

Jody Clements, C<< <clementsj at janelia.hhmi.org> >>
Rob Finn, C<< <finnr at janelia.hhmi.org> >>
Travis Wheeler, C<< <wheelert at janelia.hhmi.org> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-hmm-logo at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-HMM-Logo>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::HMM::Logo


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-HMM-Logo>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-HMM-Logo>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/Bio-HMM-Logo>

=item * Search CPAN

L<http://search.cpan.org/dist/Bio-HMM-Logo/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2012 Jody Clements.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of Bio::HMM::Logo
