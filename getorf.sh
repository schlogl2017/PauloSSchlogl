#!/bin/bash
less <&0| \
	perl -pe '/^>/?s/^>/\n>/:s/\s*$// if$.>1' | \
	perl -nse 'push @a, $_; @a = @a[@a-2..$#a]; 
	if ($. % 2 == 0){
		chomp $a[0];
		chomp $a[1]; 
		$r=qx/getorf -sequence=asis:$a[1] $o -stdout -auto 2>\/dev\/null/;
		$r =~ s/>asis/$a[0]/g;
		print $r}' -- -o=$2 | \
		perl -pe '/^>/?s/^>/\n>/:s/\s*$// if$.>1' | \
		perl -nse 'push @a, $_; @a = @a[@a-2..$#a]; 
		if ($. % 2 == 0){
			chomp $a[0];
			$a[0]=~/>(.*?) \[(\d+) - (\d+)\]\s*(.*)/g; 
			$s=(($4=~y===c)=="0")?"+":"-";
			if($f eq "f"){
				print ">".$1."_".$2."_".$3."_".$s."\n".$a[1]}
			elsif($f eq "t"){
				print $1."\t".$2."\t".$3."\t".$4."\t".$s."\t".$a[1]}}' -- -f=$1
