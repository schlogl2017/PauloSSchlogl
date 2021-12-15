#!/bin/bash

less <&0| \
	perl -pe '/^>/?s/^>/\n>/:s/\s*$// if$.>1' | \
	perl -nse 'push @a, $_; @a = @a[@a-2..$#a]; 
	if ($. % 2 == 0){
		chomp $a[0]; 
		chomp $a[1]; 
		$r=qx/compseq -sequence=asis:$a[1] $o -stdout -auto /; 
		$r=~s/#.*\n//g;$r=~s/\s+\n//g;
		$r=~s/^Word size\s+\d\nTotal count\s+\d+//g;
		$r=~s/Other.*?\n//g;
		$r=~s/\n\s+/\n/g;
		$r=~s/\h+/\t/g;
		if (defined $r and length $r){
			print substr($a[0],1)."\t".join("\n".substr($a[0],1)."\t",split("\n",$r))."\n"}}' -- -o=$1
