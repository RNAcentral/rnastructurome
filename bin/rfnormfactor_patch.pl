#!/usr/bin/env perl
# Author: Victoria Begley (begley@ebi.ac.uk)
# License: MIT
# Writes a copy of rf-normfactor (RNA Framework 2.9.8) with rf-norm's guards added; exits 3 without output if any edit no longer matches.
# Usage: rfnormfactor_patch.pl <rf-normfactor> <patched copy>
use strict;
use warnings;

my ($in, $out) = @ARGV;
die "usage: rfnormfactor_patch.pl <rf-normfactor> <patched copy>\n" unless (defined $out);
open(my $fh, "<", $in) or die "$in: $!\n";
my $src = do { local $/; <$fh> };
close($fh);

my @edits = (
    # RT-stop: skip a transcript whose log-count mean is 0 (no stops), as rf-norm does, instead of dividing by it
    [ "Ding mean guard",
      "    \$cMean = mean((map { log(\$ctrlCounts[\$_] + \$pseudoCount) } \@alli), log(max(0, \$ctrlReads) + \$pseudoCount)); \n",
      "    \$cMean = mean((map { log(\$ctrlCounts[\$_] + \$pseudoCount) } \@alli), log(max(0, \$ctrlReads) + \$pseudoCount)); \n\n                    next TRANSCRIPT if (!\$tMean || !\$cMean);\n" ],
    # denatured rates at each remaining position, so they stay aligned with \@react after the untreated filters
    [ "denatured rates aligned",
      "\@denatRatio = map { \$denatCounts[\$_] / \$denatCov[\$_] } \@i;",
      "\@denatRatio = map { my (\$p) = \$_->[0] =~ /:(\\d+)\$/; \$denatCounts[\$p] / \$denatCov[\$p] } \@react;" ],
    # the denatured filter also trims the untreated rates
    [ "untreated rates trimmed",
      "    \@denatRatio = \@denatRatio[\@i];\n                            \@react = \@react[\@i];",
      "    \@denatRatio = \@denatRatio[\@i];\n                            \@ctrlRatio = \@ctrlRatio[\@i];\n                            \@react = \@react[\@i];" ],
    # a denatured rate of 0 gives reactivity 0, as in rf-norm, instead of a division by zero
    [ "denatured zero guard",
      "max(0, (\$react[\$_]->[1] - \$ctrlRatio[\$_]) / \$denatRatio[\$_])",
      "\$denatRatio[\$_] ? max(0, (\$react[\$_]->[1] - \$ctrlRatio[\$_]) / \$denatRatio[\$_]) : 0" ],
);

for my $e (@edits) {
    my ($name, $from, $to) = @{$e};
    my $n = () = $src =~ /\Q$from\E/g;
    if ($n != 1) { print STDERR "rfnormfactor_patch.pl: '$name' matched $n times, expected 1; not patching\n"; exit 3; }
    $src =~ s/\Q$from\E/$to/;
}

open(my $oh, ">", $out) or die "$out: $!\n";
print $oh $src;
close($oh);
