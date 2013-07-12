package Bio::HMM::Logo;

use strict;
use warnings;
use JSON;
use File::Spec;
use Imager ':handy';


=head1 NAME

Bio::HMM::Logo - The great new Bio::HMM::Logo!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.02';

my $src_file   = undef;
my $typemaps   = undef;
my $hmmer_src_dir    = undef;

BEGIN {
  $src_file = __FILE__;
  $src_file =~ s/\.pm/\.c/;

  my $file = __FILE__;
  ($hmmer_src_dir) = $file =~ /^(.*)\/blib/;
  $hmmer_src_dir = File::Spec->catfile($hmmer_src_dir, 'src');
  $typemaps = __FILE__;
  $typemaps =~ s/\.pm/\.typemap/;
}

use Inline
  C        => "$src_file",
  VERSION  => '0.02',
  ENABLE   => 'AUTOWRAP',
  INC      => "-I$hmmer_src_dir/src -I$hmmer_src_dir/easel",
  LIBS     => "-L$hmmer_src_dir/easel -L$hmmer_src_dir/src -lhmmer -leasel",
  TYPEMAPS => $typemaps,
  NAME     => 'Bio::HMM::Logo';

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use Bio::HMM::Logo;

    my $foo = Bio::HMM::Logo->new();
    ...

=head1 SUBROUTINES/METHODS

=cut

=head2 hmmToLogo

=cut

sub hmmToLogo {
  my ( $hmmfile, $method ) = @_;

  if ( !$method || $method !~ /^(emission|posscore|score)$/ ) {
    $method = 'emission';
  }

  unless ( -e $hmmfile ) {
    die "$hmmfile does not exist on disk!\n";
  }
  my $hmm = inline_read_hmm($hmmfile);
  my $abc = inline_get_abc($hmm);

  my $alph     = inline_get_alphabet_string($abc);
  my @alph_arr = split( //, $alph );

  my $max_height_theoretical = 0;
  my $max_height_observed = 0;
  my $min_height_observed = 0;
  my $height_arr_ref;
  if ( $method eq "emission" ) {
    $max_height_theoretical  = inline_hmmlogo_maxHeight($abc);
    $height_arr_ref          = inline_get_emission_heights($hmm);
  }
  elsif ( $method eq "posscore" ) {
    $max_height_theoretical  = inline_hmmlogo_maxHeight($abc);
    $height_arr_ref          = inline_get_posscore_heights($hmm);
  }
  elsif ( $method eq "score" ) {
    $height_arr_ref = inline_get_score_heights($hmm);
    $max_height_theoretical  = inline_hmmlogo_maxHeight($abc);
  }

  foreach my $row (@$height_arr_ref) {
    my %char_heights;
    my $height_sum = 0;
    my $neg_height_sum = 0;
    for my $i ( 0 .. $#{$row} ) {
      $char_heights{ $alph_arr[$i] } = 0 + sprintf( "%.3f", ${$row}[$i] );
      if ($char_heights{ $alph_arr[$i] } > 0) {
        $height_sum += $char_heights{ $alph_arr[$i] } ;
      }
      else {
        $neg_height_sum += $char_heights{ $alph_arr[$i] } ;
      }
    }
    $max_height_observed = $height_sum  if ($height_sum > $max_height_observed);
    $min_height_observed = $neg_height_sum  if ($neg_height_sum < $min_height_observed);

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
    $v = 0 + sprintf( "%.0f", 100 * $v );
  }


  my $mm = inline_get_MM_array($hmm);


  my $height_data_hashref = {
    max_height_theory => $max_height_theoretical,
    max_height_obs    => $max_height_observed,
    min_height_obs    => $min_height_observed,
    height_arr        => $height_arr_ref,
    insert_probs      => $insertP,
    insert_lengths    => $insert_len,
    delete_probs      => $deleteP,
    mmline            => $mm,
  };

  #This destory was causing issues (no return string) when it occured just
  #before the return.  I think that perl would actually look after us and
  #garbage collect anyway, but this seems to solve the problem. I can not
  #explain why it was destroying the JSON object, possibly a memory violation.

 # inline_destroy_abc($abc);
 # inline_destroy_hmm($hmm);

 return $height_data_hashref;

}

=head2 hmmToLogoJson

=cut

sub hmmToLogoJson {
  my ($hmmfile, $method) = @_;

  my $height_data_hashref = hmmToLogo($hmmfile, $method);

  my $json             = JSON->new->allow_nonref;
  my $height_data_json = $json->encode($height_data_hashref);

  return $height_data_json;

}

=head2 hmmToLogoPNG

=cut

sub hmmToLogoPNG {
  my ($hmmfile, $method, $alphabet, $scaled) = @_;

  my $dna_colors = {
    'A'=> '#cbf751',
    'C'=> '#5ec0cc',
    'G'=> '#ffdf59',
    'T'=> '#b51f16',
    'U'=> '#b51f16'
  };

  my $aa_colors = {
    'A'=> '#FF9966',
    'C'=> '#009999',
    'D'=> '#FF0000',
    'E'=> '#CC0033',
    'F'=> '#00FF00',
    'G'=> '#f2f20c',
    'H'=> '#660033',
    'I'=> '#CC9933',
    'K'=> '#663300',
    'L'=> '#FF9933',
    'M'=> '#CC99CC',
    'N'=> '#336666',
    'P'=> '#0099FF',
    'Q'=> '#6666CC',
    'R'=> '#990000',
    'S'=> '#0000FF',
    'T'=> '#00FFFF',
    'V'=> '#FFCC33',
    'W'=> '#66CC66',
    'Y'=> '#006600'
  };

  my $colors = $dna_colors;
  if ($alphabet && $alphabet eq 'aa') {
    $colors = $aa_colors;
  }

  my $path = __FILE__;
  $path =~ s|[^/]*$||;

  my $regfont  = $path .'Logo/Fonts/SourceCodePro-Semibold.ttf';
  my $boldfont = $path .'Logo/Fonts/SourceCodePro-Bold.ttf';

  my $font = Imager::Font->new(file => $regfont) or die "$!\n";
  my $bold_font = Imager::Font->new(file => $boldfont) or die "$!\n";

  my $height_data_hashref = hmmToLogo($hmmfile, $method);

  #create PNG

  # determine image width and height
  my $height       = 300;
  my $column_width = 32;
  my $left_gutter  = 40;
  my $column_count = scalar @{$height_data_hashref->{height_arr}};
  my $width        = $column_count * $column_width;

  my $max_height = $height_data_hashref->{max_height_theory};
  if (defined $scaled) {
    $max_height = $height_data_hashref->{max_height_obs};
  }

  $width += $left_gutter;

  # create the image
  my $image = Imager->new(
    xsize => $width,
    ysize => $height,
  );

  $image->box(
    filled => 1,
    color  => 'white'
  );


  for (my $i = 0; $i < $column_count; $i++) {
    # draw the divider lines
    $image->line(
      color => '#EEEEEE',
      x1 => $left_gutter + ($i * $column_width),
      x2 => $left_gutter + ($i * $column_width),
      y1 => 0,
      y2 => $height,
      aa => 1,
      endp => 1
    );
    # draw the ticks
    $image->line(
      color => '#999999',
      x1 => $left_gutter + ($i * $column_width),
      x2 => $left_gutter + ($i * $column_width),
      y1 => 0,
      y2 => 5,
      aa => 1,
      endp => 1
    );
    $image->line(
      color => '#999999',
      x1 => $left_gutter + ($i * $column_width),
      x2 => $left_gutter + ($i * $column_width),
      y1 => $height - 30,
      y2 => $height - 25,
      aa => 1,
      endp => 1
    );
    $image->line(
      color => '#999999',
      x1 => $left_gutter + ($i * $column_width),
      x2 => $left_gutter + ($i * $column_width),
      y1 => $height - 15,
      y2 => $height - 10,
      aa => 1,
      endp => 1
    );

    # draw the column number
    $image->align_string(
      x => $left_gutter + ($i * $column_width) + ($column_width / 2),
      y => 2,
      font => $font,
      string => $i + 1,
      color => '#999999',
      halign => 'center',
      valign => 'top',
      size => 10,
      aa => 1
    );
    # fill in the insert odds
    my $insert_odds = $height_data_hashref->{insert_probs}[$i] / 100;
    my $insert_fill = '#ffffff';
    my $insert_text = '#666666';
    if ($insert_odds > 0.1 ) {
      $insert_fill = '#d7301f';
      $insert_text = '#ffffff';
    }
    elsif ( $insert_odds > 0.05 ) {
      $insert_fill = '#fc8d59';
    }
    elsif ( $insert_odds > 0.03 ) {
      $insert_fill = '#fdcc8a';
    }

    $image->box(
      color => $insert_fill,
      xmin => $left_gutter + ($i * $column_width) + 1,
      ymin => $height - 30,
      xmax => ($left_gutter + ($i * $column_width) + $column_width) - 1,
      ymax => $height - 15,
      filled => 1
    );

    $image->align_string(
      x => $left_gutter + ($i * $column_width) + ($column_width / 2),
      y => $height - 27,
      font => $font,
      string => $insert_odds,
      color => $insert_text,
      halign => 'center',
      valign => 'top',
      size => 10,
      aa => 1
    );
    # fill in the insert length
    my $insert_len = $height_data_hashref->{insert_lengths}[$i];
    my $length_fill = '#ffffff';
    my $length_text = '#666666';

    if ($insert_len > 9 ) {
      $length_fill = '#2171b5';
      $length_text = '#ffffff';
    }
    elsif ( $insert_len > 7 ) {
      $length_fill = '#6baed6';
    }
    elsif ( $insert_len > 4 ) {
      $length_fill = '#bdd7e7';
    }

    $image->box(
      color => $length_fill,
      xmin => $left_gutter + ($i * $column_width) + 1,
      ymin => $height - 15,
      xmax => ($left_gutter + ($i * $column_width) + $column_width) - 1,
      ymax => $height,
      filled => 1
    );
    $image->align_string(
      x => $left_gutter + ($i * $column_width) + ($column_width / 2),
      y => $height - 12,
      font => $font,
      string => $height_data_hashref->{insert_lengths}[$i],
      color => $length_text,
      halign => 'center',
      valign => 'top',
      size => 10,
      aa => 1
    );

    # draw the logo letters
    if ($height_data_hashref->{mmline}[$i] == 1) {# the column is masked
      $image->box(
        color => '#cccccc',
        xmin => $left_gutter + ($i * $column_width) + 1,
        ymin => 1,
        xmax => ($left_gutter + ($i * $column_width) + $column_width) - 1,
        ymax => $height - 30,
        filled => 1
      );
    }
    else { # column is not masked.
      my $column = $height_data_hashref->{height_arr}[$i];
      my $previous_height = 0;
      for my $letter (@$column) {
        my @values = split ':', $letter, 2;
        if ($values[1] > 0.01) { # the letter is significant enough to draw
          my $letter_color = $colors->{$values[0]};
          my $letter_height = (1 * $values[1]) / $max_height;
          my $glyph_height = $letter_height * ($height - 30);

          # there seems to be a reproducible difference between the font height
          # requested and the height that is rendered. This attempts to correct
          # that difference.
          my $fudge_factor = 1.52;

          my $bbox = $font->bounding_box(
            string => $values[0],
            size   => $glyph_height * $fudge_factor,
            sizew  => 60
          );

          $image->string(
            font => $bold_font,
            string => $values[0],
            x => $left_gutter + ($i * $column_width) + ($column_width / 2) - 16,
            y => ($height - 30) - $previous_height ,
            size => $glyph_height * $fudge_factor,
            sizew => 55,
            color => $letter_color,
            aa => 1
          );

          $previous_height += $bbox->text_height;

        }
      }
    }


  }

  # draw the axes
  # y-axes
  $image->line(
    color => '#999999',
    x1 => $left_gutter - 5,
    x2 => $width,
    y1 => 0,
    y2 => 0,
    aa => 1,
    endp => 1
  );
  $image->align_string(
    font => $font,
    string => sprintf('%.2f', $max_height),
    x => $left_gutter - 5,
    y => 0,
    size => 10,
    halign => 'right',
    valign => 'top',
    color => '#666666',
    aa => 1
  );
  $image->line(
    color => '#999999',
    x1 => $left_gutter,
    x2 => $width,
    y1 => $height - 15,
    y2 => $height - 15,
    aa => 1,
    endp => 1
  );
  $image->line(
    color => '#999999',
    x1 => $left_gutter - 5, # extend a little for the 0 tick mark
    x2 => $width,
    y1 => $height - 30,
    y2 => $height - 30,
    aa => 1,
    endp => 1
  );
  $image->align_string(
    font => $font,
    string => '0',
    x => $left_gutter - 5,
    y => $height - 30,
    size => 10,
    halign => 'right',
    valign => 'center',
    color => '#666666',
    aa => 1
  );
  $image->line(
    color => '#999999',
    x1 => $left_gutter - 5, # extend a little for the midpoint tick mark
    x2 => $left_gutter,
    y1 => ($height - 30) / 2,
    y2 => ($height - 30) / 2,
    aa => 1,
    endp => 1
  );
  $image->align_string(
    font => $font,
    string => sprintf('%.2f', $max_height / 2),
    x => $left_gutter - 5,
    y => ($height - 30) / 2,
    size => 10,
    halign => 'right',
    valign => 'center',
    color => '#666666',
    aa => 1
  );

  # x-axis
  $image->line(
    color => '#999999',
    x1 => $left_gutter,
    x2 => $left_gutter,
    y1 => 0,
    y2 => $height,
    aa => 1,
    endp => 1
  );
  my $png = undef;
  $image->write(data => \$png, type => 'png') or die $image->errstr;
  return $png;
}


######
# OO interface
#

=head2 new

=cut

sub new {
  my ($class, $args) = @_;
  my $self = {};
  bless($self, $class);
  if ($args->{hmmfile}) {
    $self->hmm_file($args->{hmmfile});
  }
  return $self;
}

=head2 hmm_file

=cut

sub hmm_file {
  my ($self, $file) = @_;
  if ($file) {
    $self->{_file} = $file;
  }
  return $self->{_file};
}

=head2 raw

=cut

sub raw {
  my ($self, $method) = @_;
  return hmmToLogo($self->hmm_file, $method);
}

=head2 as_json

=cut

sub as_json {
  my ($self, $method) = @_;
  my $height_data_hashref = hmmToLogo($self->hmm_file, $method);
  my $json             = JSON->new->allow_nonref;
  my $height_data_json = $json->encode($height_data_hashref);
  return $height_data_json;
}

=head2 as_png

=cut

sub as_png {
  my ($self, $method, $alphabet, $scaled) = @_;
  return hmmToLogoPNG($self->hmm_file, $method, $alphabet, $scaled);
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
