#!/usr/bin/perl -w
# Author: Sanzhen Liu
# Date: 3/28/2021

package dotplotpm;

sub seqnameCheck {
	my ($fasta, $list) = @_;
	my %seqname;

	# seqnames in FASTA
	open(FAS, "<", $fasta) || die;
	while (<FAS>) {
		chomp;
		if (/^>(\S+)/) {
			$seqname{$1}++;
		}
	}
	close FAS;

	# seqnames in list
	my @list = split(/,/, $list);
	
	# compare each item in the list to seqnames in FASTA
	for (@list) {
		if (! exists($seqname{$_})) {
			print STDERR "$_ does not exist in $fasta\n";
			return 0;
		}
	}
	
	# return 1 if no errors
	return 1;
}

1;
