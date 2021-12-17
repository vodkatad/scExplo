use strict;

sub revcomp {
	my $seq = shift;
	my $rev = reverse($seq);
	$rev =~ tr/ACTG/TGAC/;
	return $rev;
}

while (<STDIN>) {
	chomp;
	my @line = split("-", $_);
	my $rc = &revcomp($line[0]);
	if ($rc eq $line[1]) {
		print("YAY!\t$line[0]\t$line[1]\t$rc\n");
	} else {
		print("NUOOH!\t$line[0]\t$line[1]\t$rc\n");
	}
}
