#!/bin/bash
less <&0| \
	perl -pe '/^>/?s/^>/\n>/:s/\s*$// if$.>1' | \
	perl -nse 'push @a, $_; @a = @a[@a-2..$#a]; 
	if ($. % 2 == 0){
		chomp $a[0]; 
		chomp $a[1]; 
		$r=qx/newcpgreport -sequence=asis:$a[1] $o -stdout -auto 2>\/dev\/null/;
		$r=~s/^(.*?\n)+?FT/FT/g;
		$r=~s/FT\s+CpG island\s+(\d+)\.\.(\d+)\n/\n\1\t\2\t/g;
		$r=~s/FT\s+\/size=(\d+)\n/\1\t/g;
		$r=~s/FT\s+\/Sum C\+G=(\d+)\n/\1\t/g;
		$r=~s/FT\s+\/Percent CG=(\d+\.\d+)\n/\1\t/g;
		$r=~s/FT\s+\/ObsExp=(\d+\.\d+)\n/\1/g;
		$r=~s/FT\s+numislands\s+\d+\s*\n\/\///g;
		$r=~s/^\n//g;
		if (defined $r and length $r){
			print substr($a[0],1)."\t".join("\n".substr($a[0],1)."\t",split("\n",$r))."\n"}}' -- -o=$1
