# Title:  Make a Exon TPM, Counts, GeneLengths supermatrix for all samples
# Author: Jingxun Chen
# Date: code compiled on 20220805
# Related publication: Andrew McKay, Emma K. Costa, and Jingxun Chen, eLife, 2022
#

use strict;
use warnings;

## -------------------------------------------------------------------------------------------------------------------
my $designfile = 'GetCounts/Input/ExperimentDesign_batch_Z01-F001_F002_combined_update.csv';
my $tpmfile = 'GetCounts/Output/';
my $outct = 'Counts_Atlas_Plate01-02_240319.csv';
my $outtpm = 'TPM_Atlas_Plate01-02_240319.csv';
my $outlen = 'GeneLengths_Atlas_Plate01-02_240319.csv';

#note - might be outputting elsewhere - for me it was outputting in the R_scripts folder
## -------------------------------------------------------------------------------------------------------------------

# Open design file to see which files to read
open DN, $designfile or die $!;
my @design = <DN>;
shift @design;

my @libOrder;
my %sampleNames; # understandable names for samples
foreach (@design){
	
	my @line = split ',', $_;
	map {$_=~s/\n//g} @line;
	
	print "$line[1]\n";
	push @libOrder, $line[1];
	if ($line[0] =~/^\d+/g){$line[0] = "X$line[0]"} # Append an X if name starts with number because R won't allow it
	$sampleNames{$line[1]} = $line[0]; # understandable names for samples
	print "Sample name: $line[0]\n";
}
my $lastlib = $libOrder[-1];
print "Last library: $lastlib\n";

my %SuperMatrix;
my %Idlist;

foreach my $lib (@libOrder){
	
	#print "$tpmfile/$lib\_counts.csv\t";
	open FH, "$tpmfile/$lib\_counts.csv" or die "$tpmfile/$lib\_counts.csv ", $!;
	my @file = <FH>;
	close (FH);	
	shift @file;
	print scalar @file, "\n";
	
	foreach (@file){

		my @line = split ",", $_;
		map {$_=~s/\n//g} @line;	
		die "Lengths differ at $line[0] $lib" if ((exists $Idlist{$line[0]}) && $Idlist{$line[0]} != $line[1]);
		$Idlist{$line[0]} = $line[1]; # global id list common to all samples and their lengths
		$SuperMatrix{$line[0]}{$lib} = [$line[2], $line[3]]; # Counts and TPM
	} 
}

print "Total gene: ", scalar keys %Idlist, "\n";


# print to a file
open CT, ">$outct" or die $!;
open TPM, ">$outtpm" or die $!;
open LEN, ">$outlen" or die $!;

# Print header with proper sample names
foreach (@libOrder){
	if ($_ eq $lastlib){ # If last library print without comma
		print CT "\"$sampleNames{$_}\""; print TPM "\"$sampleNames{$_}\""; print LEN "\"$sampleNames{$_}\"";
	}
	else {
		print CT "\"$sampleNames{$_}\","; print TPM "\"$sampleNames{$_}\","; print LEN "\"$sampleNames{$_}\",";
	}
}
print "Made it to line 80\n";
print CT "\n"; print TPM "\n"; print LEN "\n";

foreach my $gene (keys %Idlist){
	
	print CT "\"$gene\","; print LEN "\"$gene\",";	print TPM "\"$gene\",";
	
	foreach my $lib (@libOrder){
		
		if ($lib eq $lastlib){ # If last library print without comma
			if (exists $SuperMatrix{$gene}{$lib}){
				print LEN "$Idlist{$gene}";
				print CT "$SuperMatrix{$gene}{$lib}[0]";
				print TPM "$SuperMatrix{$gene}{$lib}[1]";

			}
			else {
				print LEN "$Idlist{$gene},";
				print CT "0,";
				print TPM "0,";			
			}
		}
		else {
			if (exists $SuperMatrix{$gene}{$lib}){
				print LEN "$Idlist{$gene},";
				print CT "$SuperMatrix{$gene}{$lib}[0],";
				print TPM "$SuperMatrix{$gene}{$lib}[1],";

			}
			else {
				print LEN "$Idlist{$gene},";
				print CT "0,";
				print TPM "0,";
			}
		}
	}
print CT "\n"; print TPM "\n"; print LEN "\n";
}
close (CT); close (TPM); close (LEN);




