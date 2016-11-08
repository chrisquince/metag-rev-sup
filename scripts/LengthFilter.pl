#!/usr/bin/perl

my $fastaFile = $ARGV[0];
my $minLength = $ARGV[1];
my @Seq = ();
my @id       = ();

my $count = 0;

open(FILE, $fastaFile) or die "Can't open $fastaFile\n";

my $seq = "";

while($line = <FILE>){ 
    chomp($line);
    
    if($line =~ />(.*)/){
	
	$id[$count] = $1;
	
	if($seq ne ""){
	    $Seq[$count - 1] = $seq;

	    $seq = "";
	}

	$count++;
    }
    else{
	$seq .= $line;
    }
}

$Seq[$count - 1] = $seq;
$total = $count;



for($i = 0; $i < $total; $i++){
    $seqlength = length($Seq[$i]);
    if($seqlength > $minLength){
    	print ">$id[$i]\n";

    	print "$Seq[$i]\n";
    }
}
