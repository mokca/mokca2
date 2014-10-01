#!/usr/bin/perl

use Getopt::Long;
use FileHandle;

STDERR->autoflush(1);

# Ticker

my $tick_c = 0;
my $col_wrap = 78;
my $tick_step = 1000;
my $tick_char = '.';

# Explain purpose

print STDERR "::: Split Uniprot XML - Split Uniprot XML files for parallel processing\n";

# Setup universal paths

use path_setup;

# Local uniprot files

my $sprot_fn = $mokcatanic_data_root . "uniprot_sprot.xml.gz";
my $trembl_fn = $mokcatanic_local_root . "/uniprot_trembl.xml.gz";
my $split_base_fn = $mokcatanic_local_root . "/uniprot_split_";
my $split_base_c = 36;

# Create XML headers in all split files

for (my $ii=0;$ii<$split_base_c;$ii++) {
    my $split_fn = $split_base_fn . $ii . ".xml";
    open my $split_fh, ">", $split_fn || die "\n+++ Cannot open $split_fn for writing: $?\n";
    print $split_fh "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" . 
    "<uniprot xmlns=\"http://uniprot.org/uniprot\"\n" . 
    " xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n" .
    " xsi:schemaLocation=\"http://uniprot.org/uniprot http://www.uniprot.org/support/docs/uniprot.xsd\">\n";
    close $split_fh || die "\n+++ Cannot close $split_fn ($split_fh): $?\n";
}

# Proccess sprot then trembl files

open my $uniprot_fh, "-|", "gzip -dc < $sprot_fn" || die "\n+++ Cannot open $sprot_fn to read data: $?\n";
split_xml($uniprot_fh, $split_base_fn, $split_base_c);
close $uniprot_fh || die "\n+++ Cannot close $sprot_fn to stop reading data\n";

open $uniprot_fh, "-|", "gzip -dc < $trembl_fn" || die "\n+++ Cannot open $sprot_fn to read data: $?\n";
split_xml($uniprot_fh);
close $uniprot_fh || die "\n+++ Cannot close $trembl_fn to stop reading data: $?\n";

# Finish xml in all split files

for (my $ii=0;$ii<$split_base_c;$ii++) {
    my $split_fn = $split_base_fn . $ii . ".xml";
    open my $split_fh, ">>", $split_fn || die "\n+++ Cannot open $split_fn for writing: $?\n";
    print $split_fh "</uniprot>\n";
    close $split_fh || die "\n+++ Cannot close $split_fh: $?\n";
}

sub split_xml {
    ( my $current_fh )= @_;

    my $entry_buffer = "";
    my $in_org_def = 0;
    my $is_human = 0;
    my $muted = 0;
    my $entry_c = 0;
    my $h_entry_c = 0;
    my $slice_size = 100;
    
    # Generate first split filename
    
    my $split_fc = 0;
    my $split_fn = $split_base_fn . $split_fc . ".xml";
    open my $split_fh, ">>", $split_fn || die "\n+++ Cannot open $split_fn for appending: $?\n";
    
    # For each line, accumulate in buffer until end of entry detected
    
    while ($line = <$current_fh>) {
        
        # Do not accumulate references or evidence in buffer, some entries have
        # a pathologically large number of these features
        
        if ($line =~ /<reference/ || $line=~ /<evidence type="ECO:0000313"/) {
            $muted = 1;
        } elsif (! $muted) {
            $entry_buffer .= $line;
        }
        
        if (($line =~ /<\/reference/ || $line=~ /<\/evidence>/) && $muted == 1) {
            $muted = 0;
        }
        
        if ($line =~ /<organism>/ || $line =~ /<organism evidence/) {
            $in_org_def = 1;
        }
        
        if ($line =~ /<\/organism>/) {
            $in_org_def = 0;
        }
        
        # Flag human entries as being human, only these are added to split files
        
        if ($in_org_def == 1 && $line =~ /<name type="scientific">Homo sapiens<\/name>/) {
            $is_human = 1;
        }
        
        # On seeing the end of an entry, add to a split file if human
        
        if ($line eq "</entry>\n") {
            if ($is_human) {
                print $split_fh $entry_buffer;
                $h_entry_c++;
                
                # Pointless ticker
                
                if ($tick_c % $tick_step == 0) {
                    print STDERR $tick_char;
                    if ($tick_c % ($tick_step * $col_wrap) == 0 && $tick_c) {
                        print STDERR "\n";
                    }
                }
                $tick_c++;
                
                if (($h_entry_c % $slice_size) == 0) {
                    
                    # After a slice, close current split file & open next
                    
                    close $split_fh || die "\n+++ Cannot close $split_fh: $?\n";
                    
                    $split_fc++;
                    $split_fc %= $split_base_c;
                    $split_fn = $split_base_fn . $split_fc . ".xml";
                    
                    open $split_fh, ">>", $split_fn || die "\n+++ Cannot open $split_fn for appending: $?\n";
                }
            }
            $entry_buffer = "";
            $is_human = 0;
            $muted = 0;
            $entry_c++;
        }
    }
    
    # Tidy up
    
    close $split_fh || die "\n+++ Cannot close $split_fh: $?\n";
    
    print STDERR "\n--- File finished: $entry_c entries, $h_entry_c human\n";
    $tick_c = 0;
}