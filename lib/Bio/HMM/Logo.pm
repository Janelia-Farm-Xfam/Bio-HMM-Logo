package Bio::HMM::Logo;

use strict;
use warnings;
use JSON;
use File::Spec;
use Imager ':handy';
use SVG;


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

# as defined in esl_alphabet.h
##define eslUNKNOWN     0
##define eslRNA         1
##define eslDNA         2
##define eslAMINO       3
my @alphabet = qw( unk rna dna aa );

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

  if ( !$method || $method !~ /^(entropy_all|entropy_above|score)$/ ) {
    $method = 'entropy_all';
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
  my $height_arr_ref = undef;
  my $prob_arr_ref = undef;
  if ( $method eq "entropy_all" ) {
    $max_height_theoretical  = inline_hmmlogo_maxHeight($abc);
    my $arr_pair_ref         = inline_get_relative_entropy_all($hmm);
    ($height_arr_ref,$prob_arr_ref) = @$arr_pair_ref;
  }
  elsif ( $method eq "entropy_above" ) {
    $max_height_theoretical  = inline_hmmlogo_maxHeight($abc);
    my $arr_pair_ref         = inline_get_relative_entropy_above_bg($hmm);
    ($height_arr_ref,$prob_arr_ref) = @$arr_pair_ref;
  }
  elsif ( $method eq "score" ) {
    $height_arr_ref = inline_get_score_heights($hmm);
    $max_height_theoretical  = -1; # this field is not meaningful for the "score" method
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

  foreach my $p_row (@$prob_arr_ref) {
    my %char_probs = ();
    for my $i ( 0 .. $#{$p_row} ) {
      $char_probs{ $alph_arr[$i] } = sprintf( "%.3f", ${$p_row}[$i] );
    }
    my @sorted_keys =
      sort { $char_probs{$a} <=> $char_probs{$b} } keys %char_probs;

    for my $i ( 0 .. $#{$p_row} ) {
      my $key = $sorted_keys[$i];
      ${$p_row}[$i] = "$key:$char_probs{$key}";
    }
  }

  my $insertP = inline_get_insertP($hmm);
  foreach my $v (@$insertP) {
    $v = sprintf( "%.2f", $v );
  }

  my $insert_len = inline_get_insertLengths($hmm);
  foreach my $v (@$insert_len) {
    $v = 0 + sprintf( "%.0f", $v );
  }

  my $deleteP = inline_get_deleteP($hmm);
  foreach my $v (@$deleteP) {
    $v = sprintf( "%.2f", $v );
  }


  my $abc_type = inline_get_abc_type($hmm);
  my $mm = inline_get_MM_array($hmm);


  my $height_data_hashref = {
    alphabet          => $alphabet[$abc_type],
    max_height_theory => $max_height_theoretical,
    max_height_obs    => $max_height_observed,
    min_height_obs    => "$min_height_observed",
    height_arr        => $height_arr_ref,
    insert_probs      => $insertP,
    insert_lengths    => $insert_len,
    delete_probs      => $deleteP,
    mmline            => $mm,
    height_calc       => $method,
  };

  if ($prob_arr_ref) {
    $height_data_hashref->{probs_arr} = $prob_arr_ref;
  }

  #This destroy was causing issues (no return string) when it occurred just
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
  my ($hmmfile, $method, $scaled) = @_;
  my $height_data_hashref = hmmToLogo($hmmfile, $method);
  return _build_png($height_data_hashref, $scaled);
}

=head2 hmmToLogoSVG

=cut

sub hmmToLogoSVG {
  my ($hmmfile, $method, $scaled) = @_;
  my $height_data_hashref = hmmToLogo($hmmfile, $method);
  return _build_svg($height_data_hashref, $scaled);
}


=head2 _dna_colors

=cut

sub _dna_colors {
  return {
    'A'=> '#cbf751',
    'C'=> '#5ec0cc',
    'G'=> '#ffdf59',
    'T'=> '#b51f16',
    'U'=> '#b51f16'
  };
}

=head2 _aa_colors

=cut

sub _aa_colors {
  return {
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
}


=head2 _colors_by_alphabet

=cut

sub _colors_by_alphabet {
  my $alphabet = shift;
  if ($alphabet && $alphabet eq 'aa') {
    return _aa_colors;
  }
  else {
    return _dna_colors;
  }
}

=head2 _image_height

=cut

sub _image_height {
  return 300;
}

=head2 _build_png

=cut

sub _build_png {
  my ($height_data_hashref, $scaled, $debug) = @_;
  my $alphabet = (exists $height_data_hashref->{alphabet}) ?
                   $height_data_hashref->{alphabet} : 'dna';


  my $colors = _colors_by_alphabet($alphabet);

  my $path = __FILE__;
  $path =~ s|[^/]*$||;

  my $regfont  = $path .'Logo/Fonts/SourceCodePro-Semibold.ttf';
  my $boldfont = $path .'Logo/Fonts/SourceCodePro-Bold.ttf';

  my $font = Imager::Font->new(file => $regfont, color => '#eeeeee') or die "$!\n";
  my $bold_font = Imager::Font->new(file => $boldfont, color => '#eeeeee') or die "$!\n";

  #create PNG

  # determine image width and height
  my $height       = _image_height;
  my $column_width = 32;
  my $left_gutter  = 40;
  my $column_count = scalar @{$height_data_hashref->{height_arr}};
  my $width        = $column_count * $column_width;



  my $max_height = $height_data_hashref->{max_height_obs};
  if (defined $scaled) {
    if (exists $height_data_hashref->{height_calc}
      && $height_data_hashref->{height_calc} eq 'emission') {
      $max_height = $height_data_hashref->{max_height_theory};
    }
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
    # draw delete probability section ticks
    $image->line(
      color => '#999999',
      x1 => $left_gutter + ($i * $column_width),
      x2 => $left_gutter + ($i * $column_width),
      y1 => $height - 45,
      y2 => $height - 40,
      aa => 1,
      endp => 1
    );
    # draw insert probability section ticks
    $image->line(
      color => '#999999',
      x1 => $left_gutter + ($i * $column_width),
      x2 => $left_gutter + ($i * $column_width),
      y1 => $height - 30,
      y2 => $height - 25,
      aa => 1,
      endp => 1
    );
    # draw insert length section ticks
    $image->line(
      color => '#999999',
      x1 => $left_gutter + ($i * $column_width),
      x2 => $left_gutter + ($i * $column_width),
      y1 => $height - 15,
      y2 => $height - 10,
      aa => 1,
      endp => 1
    );




    # fill in the delete odds
    my $delete_odds = $height_data_hashref->{delete_probs}[$i];
    my $delete_fill = '#ffffff';
    my $delete_text = '#666666';

    if ( $delete_odds < 0.95 ) {
      $delete_fill = '#bdd7e7';
    }
    elsif ( $delete_odds < 0.85 ) {
      $delete_fill = '#6baed6';
    }
    elsif ($delete_odds < 0.75 ) {
      $delete_fill = '#2171b5';
      $delete_text = '#ffffff';
    }

    $image->box(
      color => $delete_fill,
      xmin => $left_gutter + ($i * $column_width) + 1,
      ymin => $height - 45,
      xmax => ($left_gutter + ($i * $column_width) + $column_width) - 1,
      ymax => $height - 30,
      filled => 1
    );
    $image->align_string(
      x => $left_gutter + ($i * $column_width) + ($column_width / 2),
      y => $height - 42,
      font => $font,
      string => $delete_odds,
      color => $delete_text,
      halign => 'center',
      valign => 'top',
      size => 10,
      aa => 1
    );



    # fill in the insert odds
    my $insert_odds = $height_data_hashref->{insert_probs}[$i];
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
      $length_fill = '#d7301f';
      $length_text = '#ffffff';
    }
    elsif ( $insert_len > 7 ) {
      $length_fill = '#fc8d59';
    }
    elsif ( $insert_len > 4 ) {
      $length_fill = '#fdcc8a';
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
        ymax => $height - 45,
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
          my $glyph_height = $letter_height * ($height - 45);

          # there seems to be a reproducible difference between the font height
          # requested and the height that is rendered. This attempts to correct
          # that difference.
          my $fudge_factor = 1.52;

          # the Q in the font we use has a large descender that makes it taller
          # than it should be, so we have to offset that by making it proportionally
          # smaller here.
          if ($values[0] eq 'Q') {
            $fudge_factor = 1.18;
          }
          if ($values[0] =~ /C|G|S|O/) {
            $fudge_factor = 1.46;
          }
          if ($values[0] =~ /J|U/) {
            $fudge_factor = 1.48;
          }

          my $bbox = $font->bounding_box(
            string => $values[0],
            size   => $glyph_height * $fudge_factor,
            sizew  => 60
          );

          if ($debug) {
            my $xmin = $left_gutter + ($i * $column_width);
            my $ymax = ($height - 45) - $previous_height;
            $image->box(
              color => $colors->{'Q'},
              xmin => $xmin,
              xmax => $xmin + $bbox->display_width,
              ymin => $ymax - $bbox->text_height,
              ymax => $ymax,
            );
          }

          my $x = $left_gutter + ($i * $column_width);
          my $y = ($height - 45) - $previous_height - $bbox->text_height;

          $image->string(
            font => $bold_font,
            string => $values[0],
            x => $x,
            align => 0,
            y => $y,
            size => $glyph_height * $fudge_factor,
            sizew => 55,
            color => $letter_color,
            aa => 1,
          );

          $previous_height += $bbox->text_height;

        }
      }
    }

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
  # draw the line above the insert length section
  $image->line(
    color => '#999999',
    x1 => $left_gutter,
    x2 => $width,
    y1 => $height - 15,
    y2 => $height - 15,
    aa => 1,
    endp => 1
  );
  # draw the line above the insert probability section
  $image->line(
    color => '#999999',
    x1 => $left_gutter,
    x2 => $width,
    y1 => $height - 30,
    y2 => $height - 30,
    aa => 1,
    endp => 1
  );
  # draw the line above the delete probability section
  $image->line(
    color => '#999999',
    x1 => $left_gutter - 5, # extend a little for the 0 tick mark
    x2 => $width,
    y1 => $height - 45,
    y2 => $height - 45,
    aa => 1,
    endp => 1
  );
  $image->align_string(
    font => $font,
    string => '0',
    x => $left_gutter - 5,
    y => $height - 45,
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
    y1 => ($height - 45) / 2,
    y2 => ($height - 45) / 2,
    aa => 1,
    endp => 1
  );
  $image->align_string(
    font => $font,
    string => sprintf('%.2f', $max_height / 2),
    x => $left_gutter - 5,
    y => ($height - 45) / 2,
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


=head2 _build_svg

=cut

sub _build_svg {
  my ($height_data_hashref, $scaled, $debug) = @_;
  my $alphabet = (exists $height_data_hashref->{alphabet}) ?
                   $height_data_hashref->{alphabet} : 'dna';

  my $colors = _colors_by_alphabet($alphabet);

  # determine image width and height
  my $height       = _image_height;
  my $column_width = 32;
  my $left_gutter  = 40;
  my $column_count = scalar @{$height_data_hashref->{height_arr}};
  my $width        = $column_count * $column_width;



  my $max_height = $height_data_hashref->{max_height_obs};
  if (defined $scaled) {
    if (exists $height_data_hashref->{height_calc}
      && $height_data_hashref->{height_calc} eq 'emission') {
      $max_height = $height_data_hashref->{max_height_theory};
    }
  }

  $width += $left_gutter;
  # create svg version and render it here.
  my $svg = SVG->new(width => $width, height => $height);

  # clear the background and set as white.
  $svg->rect(
    x      => 0,
    y      => 0,
    width  => $width,
    height => $height,
    style  => {
      fill => '#fff',
    }
  );

  for (my $i = 0; $i < $column_count; $i++) {
    # draw the divider lines
    $svg->line(
      x1 => $left_gutter + ($i * $column_width),
      x2 => $left_gutter + ($i * $column_width),
      y1 => 0,
      y2 => $height,
      style => {
        stroke => '#ddd',
      }
    );

    # draw the ticks
    $svg->line(
      x1 => $left_gutter + ($i * $column_width),
      x2 => $left_gutter + ($i * $column_width),
      y1 => 0,
      y2 => 5,
      style => {
        stroke => '#999',
      }
    );
    # draw delete probability section ticks
    $svg->line(
      x1 => $left_gutter + ($i * $column_width),
      x2 => $left_gutter + ($i * $column_width),
      y1 => $height - 45,
      y2 => $height - 40,
      style => {
        stroke => '#999',
      }
    );
    # draw insert probability section ticks
    $svg->line(
      x1 => $left_gutter + ($i * $column_width),
      x2 => $left_gutter + ($i * $column_width),
      y1 => $height - 30,
      y2 => $height - 25,
      style => {
        stroke => '#999',
      }
    );
    # draw insert length section ticks
    $svg->line(
      x1 => $left_gutter + ($i * $column_width),
      x2 => $left_gutter + ($i * $column_width),
      y1 => $height - 15,
      y2 => $height - 10,
      style => {
        stroke => '#999',
      }
    );

    # fill in the delete odds
    my $delete_odds = $height_data_hashref->{delete_probs}[$i];
    my $delete_fill = '#ffffff';
    my $delete_text = '#666666';

    if ( $delete_odds < 0.95 ) {
      $delete_fill = '#bdd7e7';
    }
    elsif ( $delete_odds < 0.85 ) {
      $delete_fill = '#6baed6';
    }
    elsif ($delete_odds < 0.75 ) {
      $delete_fill = '#2171b5';
      $delete_text = '#ffffff';
    }

    $svg->rect(
      x => $left_gutter + ($i * $column_width) + 1,
      y => $height - 45,
      width => $column_width,
      height => 15,
      style  => {
        fill => $delete_fill,
      }
    );
    $svg->text(
      x => $left_gutter + ($i * $column_width) + ($column_width / 2),
      y => $height - 34,
      style => {
        'font-family' => 'Arial',
        'font-size'   => '10px',
        'fill'        => $delete_text,
        'text-anchor' => 'middle',
      }
    )->cdata($delete_odds);


    # fill in the insert odds
    my $insert_odds = $height_data_hashref->{insert_probs}[$i];
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

    $svg->rect(
      x => $left_gutter + ($i * $column_width) + 1,
      y => $height - 30,
      width => $column_width,
      height => 15,
      style  => {
        fill => $insert_fill,
      }
    );
    $svg->text(
      x => $left_gutter + ($i * $column_width) + ($column_width / 2),
      y => $height - 19,
      style => {
        'font-family' => 'Arial',
        'font-size'   => '10px',
        'fill'        => $insert_text,
        'text-anchor' => 'middle',
      }
    )->cdata($insert_odds);

    # fill in the insert length
    my $insert_len = $height_data_hashref->{insert_lengths}[$i];
    my $length_fill = '#ffffff';
    my $length_text = '#666666';

    if ($insert_len > 9 ) {
      $length_fill = '#d7301f';
      $length_text = '#ffffff';
    }
    elsif ( $insert_len > 7 ) {
      $length_fill = '#fc8d59';
    }
    elsif ( $insert_len > 4 ) {
      $length_fill = '#fdcc8a';
    }

    $svg->rect(
      x => $left_gutter + ($i * $column_width) + 1,
      y => $height - 15,
      width => $column_width,
      height => 15,
      style  => {
        fill => $length_fill,
      }
    );
    $svg->text(
      x => $left_gutter + ($i * $column_width) + ($column_width / 2),
      y => $height - 4,
      style => {
        'font-family' => 'Arial',
        'font-size'   => '10px',
        'fill'        => $length_text,
        'text-anchor' => 'middle',
      }
    )->cdata($height_data_hashref->{insert_lengths}[$i]);

    # draw the letters
    if ($height_data_hashref->{mmline}[$i] == 1) {# the column is masked
      $svg->rect(
        x => $left_gutter + ($i * $column_width) + 1,
        y => 1,
        width => $column_width,
        height => $height - 45,
        filled => 1,
        style  => {
          fill => '#ccc',
        }
      );
    }
    else {
      # TODO: Letter drawing goes here
      my $column = $height_data_hashref->{height_arr}[$i];
      my $previous_height = 0;
      my @coordinates = ();
      for my $letter (@$column) {
        my @values = split ':', $letter, 2;
        if ($values[1] > 0.01) { # the letter is significant enough to draw
          my $letter_color = $colors->{$values[0]};
          my $letter_height = (1 * $values[1]) / $max_height;
          my $glyph_height = $letter_height * ($height - 45);

          my $x = $left_gutter + ($i * $column_width) + ($column_width / 2);
          my $y = ($height - 45) - $previous_height;

          my $rect_y = $y;

          my $h_ratio = $glyph_height / 60;
          my $w_ratio = $column_width / 70;

          my $font_size = '85px';

          if($values[0] =~ /G|C/) {
            $font_size = '80px';

            $y -= ($glyph_height * 2) / 100;
          }

          push @coordinates, {
            x       => $x,
            y       => $y,
            f_size  => $font_size,
            w_ratio => $w_ratio,
            h_ratio => $h_ratio,
            color   => $letter_color,
            char    => $values[0],
          };

          if ($debug) {
            $svg->rect(
              x => $left_gutter + ($i * $column_width),
              y => $rect_y - $glyph_height,
              width => $column_width,
              height => $glyph_height,
              style => {
                stroke => $letter_color,
                'fill-opacity' => 0,
              }
            );
          }

          $previous_height += $glyph_height;

        }
      }

      # Draw the letters in reverse order so that the larger letters
      # don't overlap the smaller ones below. This makes it easier to
      # read the smaller letters.

      for my $letter (reverse @coordinates) {
        $svg->text(
          x => 0,
          y => 0,
          style => {
            'font-family' => 'Arial',
            'font-size'   => $letter->{f_size},
            'font-weight' => 'bold',
            'fill'        => $letter->{color},
            'text-anchor' => 'middle',
          },
          transform => "matrix($letter->{w_ratio}, 0, 0, $letter->{h_ratio}, $letter->{x}, $letter->{y})",

        )->cdata($letter->{char});
      }

    }

    #draw the column number
    $svg->text(
      x => $left_gutter + ($i * $column_width) + ($column_width / 2),
      y => 10,
      style => {
        'font-family' => 'Arial',
        'font-size'   => '10px',
        'fill'        => '#999',
        'text-anchor' => 'middle',
      }
    )->cdata($i + 1);

  }

  # draw the axes
  # y-axes
  $svg->line(
    x1 => $left_gutter - 5,
    x2 => $width,
    y1 => 0,
    y2 => 0,
    style => {
      stroke => '#999',
    }
  );
  $svg->text(
    x => $left_gutter - 5,
    y => 8,
    style => {
      'font-family' => 'Arial',
      'font-size'   => '10px',
      'fill'        => '#666',
      'text-anchor' => 'end',
    }
  )->cdata(sprintf('%.2f', $max_height));
  # draw the line above the insert length section
  $svg->line(
    x1 => $left_gutter,
    x2 => $width,
    y1 => $height - 15,
    y2 => $height - 15,
    style => {
      stroke => '#999',
    }
  );
  # draw the line above the insert probability section
  $svg->line(
    x1 => $left_gutter,
    x2 => $width,
    y1 => $height - 30,
    y2 => $height - 30,
    style => {
      stroke => '#999',
    }
  );
  # draw the line above the delete probability section
  $svg->line(
    x1 => $left_gutter - 5, # extend a little for the 0 tick mark
    x2 => $width,
    y1 => $height - 45,
    y2 => $height - 45,
    style => {
      stroke => '#999',
    }
  );
  $svg->text(
    x => $left_gutter - 5,
    y => $height - 45,
    style => {
      'font-family' => 'Arial',
      'font-size'   => '10px',
      'fill'        => '#666',
      'text-anchor' => 'end',
    }
  )->cdata('0');
  $svg->line(
    x1 => $left_gutter - 5, # extend a little for the midpoint tick mark
    x2 => $left_gutter,
    y1 => ($height - 45) / 2,
    y2 => ($height - 45) / 2,
    style => {
      stroke => '#999',
    }
  );
  $svg->text(
    x => $left_gutter - 5,
    y => (($height - 45) / 2) + 3,
    style => {
      'font-family' => 'Arial',
      'font-size'   => '10px',
      'fill'        => '#666',
      'text-anchor' => 'end',
    }
  )->cdata(sprintf('%.2f', $max_height / 2));

  # x-axis
  $svg->line(
    x1 => $left_gutter,
    x2 => $left_gutter,
    y1 => 0,
    y2 => $height,
    style => {
      stroke => '#999',
    }
  );

  return $svg->xmlify;
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
  my $json                = JSON->new->allow_nonref;
  my $height_data_json    = $json->encode($height_data_hashref);
  return $height_data_json;
}

=head2 as_png

=cut

sub as_png {
  my ($self, $method, $scaled) = @_;
  return hmmToLogoPNG($self->hmm_file, $method, $scaled);
}

=head2 as_svg

=cut

sub as_svg {
  my ($self, $method, $scaled) = @_;
  return hmmToLogoSVG($self->hmm_file, $method, $scaled);
}

=head2 dl_load_flags
=head2 inline_destroy_abc
=head2 inline_destroy_hmm
=head2 inline_get_MM_array
=head2 inline_get_abc
=head2 inline_get_abc_type
=head2 inline_get_alphabet_string
=head2 inline_get_deleteP
=head2 inline_get_relative_entropy_all
=head2 inline_get_insertLengths
=head2 inline_get_insertP
=head2 inline_get_relative_entropy_above_bg
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
