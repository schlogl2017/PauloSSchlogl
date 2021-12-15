#!/bin/bash
less <&0| \
	perl -pe '/^>/?s/^>/\n>/:s/\s*$// if$.>1' | \
	perl -nse 'push @a, $_; @a = @a[@a-2..$#a]; 
	if ($. % 2 == 0){
		chomp $a[0]; 
		chomp $a[1]; 
		$r=qx/freak -seqall=asis:$a[1] $o -stdout -auto 2>\/dev\/null/;
		$r=~s/\s+\n/\n/g; 
		$r=~s/FREAK.*?\n//g;
		$r=~s/^\s+//g;
		$r=~s/\n\s+/\n/g;
		$r=~s/\h+/\t/g;
		if (defined $r and length $r){
			print substr($a[0],1)."\t".join("\n".substr($a[0],1)."\t",split("\n",$r))."\n"}}' -- -o=$1
