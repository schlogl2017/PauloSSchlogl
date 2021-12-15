#!/bin/bash
less <&0| \
	perl -ane '$r{$F[0].":".$F[1]}=$F[2];
		unless($F[0]~~@s){
			push @s,$F[0];}
		unless($F[1]~~@m){
			push @m,$F[1];}
	END{
	print "Contigs\t".join("\t",@s)."\n";
	for($i=0;$i<@m;$i++){
		print $m[$i];
		for($j=0;$j<@s;$j++){
			(not defined $r{$s[$j].":".$m[$i]})?print "\t".0:print"\t".$r{$s[$j].":".$m[$i]};}
		print "\n";}}'
