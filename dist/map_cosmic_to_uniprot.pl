#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use FileHandle;
use DBI();
use List::MoreUtils qw(uniq);

# Autoflush STDERR

STDERR->autoflush(1);

# Ticker

my $tick_c = 1;
my $col_wrap = 78;
my $tick_step = 1000;
my $tick_char = '.';

# Explain purpose

print STDERR "::: Map COSMIC to Uniprot - use mysteriously arcane rules to map genes to proteins\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

my $cosmic_fasta_fn = $mokcatanic_data_root . "All_COSMIC_Genes.fasta.gz";
my $biomart_fn = $mokcatanic_data_root . "biomart.tsv";

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_populate_map_table = $dbh_titanic->prepare(q{INSERT INTO cosmic_vs_uniprot VALUES(?,NULL,0.0)});
my $dbq_uniprot_id = $dbh_titanic->prepare(q{SELECT uniprot_id FROM uniprot_accession WHERE accession = ?});
my $dbq_make_mapping = $dbh_titanic->prepare(q{UPDATE cosmic_vs_uniprot SET uniprot_id = ?, identity = 100.0 WHERE cosmic_id = ?}); # The only maps make in this script should be 100% identical
my $dbq_cosmic_id = $dbh_titanic->prepare(q{SELECT cosmic_id FROM cosmic_hgnc WHERE cosmic_gene_name = ?});
my $dbq_unmatched_entries = $dbh_titanic->prepare(q{SELECT cosmic_id FROM cosmic_vs_uniprot WHERE uniprot_id IS NULL});
my $dbq_cosmic_sequence = $dbh_titanic->prepare(q{SELECT sequence, sequence_md5 FROM sequence WHERE source='cosmic' AND codes_for='protein' AND source_id = ?});
my $dbq_uniprot_md5_match_c = $dbh_titanic->prepare(q{SELECT COUNT(*) FROM sequence WHERE source='uniprot' AND codes_for='protein' AND isoform=1 AND sequence_md5=?});
my $dbq_uniprot_sequence_by_md5 = $dbh_titanic->prepare(q{SELECT source_id, sequence FROM sequence WHERE source='uniprot' AND codes_for='protein' AND isoform=1 AND sequence_md5=?});

# Go

# Read all of the biomart mappings into an array so we can search through them

print STDERR "--- Loading BioMart data\n";
my (@biomart, $BIOMART);
open($BIOMART, "<$biomart_fn") || die "\n+++ Cannot open $biomart_fn to read data: $?\n";
while (<$BIOMART>) {
    chomp;
    push @biomart, [ split /\t/];
    
    # Pointless ticker
    
    if ($tick_c % $tick_step == 0) {
        print STDERR $tick_char;
        if ($tick_c % ($tick_step * $col_wrap) == 0 && $tick_c) {
            print STDERR "\n";
        }
    }
    $tick_c++;
}
close($BIOMART) || die "\n+++ Cannot close $biomart_fn: $?\n";

# First pass: run through COSMIC file containing FASTAs and - conveniently - ENSTS to match
# BioMart ENSTs to COSMIC ENSTs and hope they give us a unique SPROT or TREMBL accession

$tick_c = 1;
print STDERR "\n--- Matching COSMIC names using BioMart\n";

my %unmapped_cosmic_ids;
my ($cosmic_c, $missing_c, $multiple_c, $matched_c, $cosmic_hgnc_name, $cosmic_enst, $fasta_fh, $cosmic_id);
$cosmic_c = 0;
$missing_c = 0;
$multiple_c = 0;
$matched_c = 0;
open $fasta_fh, "gunzip -c $cosmic_fasta_fn |" || die "\n+++ Cannot open $cosmic_fasta_fn to read data: $?\n";
while (my $line = <$fasta_fh>) {
    chomp($line);
    if ($line =~ /^>(\w+)\s(\w+)/) { # We have a FASTA sequence header
        $cosmic_hgnc_name = $1;
        $cosmic_enst = $2;
        $cosmic_c++;
        
        # Grab the cosmic ID
        
        $cosmic_id = 0;
        $dbq_cosmic_id->execute($cosmic_hgnc_name);
        ( $cosmic_id ) = $dbq_cosmic_id->fetchrow_array;
        $dbq_cosmic_id->finish;
        if (! $cosmic_id) {
            die "\n+++ Cannot get numeric COSMIC ID for $cosmic_hgnc_name\n";
        }
        
        # Put the cosmic ID into cosmic_vs_uniprot
        
        $dbq_populate_map_table->execute($cosmic_id);
        
        # Iterate through the BioMart data and pull out all SPROT and TREMBL accessions that are mapped
        
        my ( @trembl_accs, @sprot_accs );
        
        my ($i, $biomart_enst, $biomart_hgnc_name, $biomart_trembl, $biomart_sprot);
        for $i ( 0 .. $#biomart ) {
            $biomart_enst = $biomart[$i][1];
            $biomart_hgnc_name = $biomart[$i][4];
            $biomart_trembl = $biomart[$i][2];
            $biomart_sprot = $biomart[$i][3];
            if (!$biomart_hgnc_name) {
                $biomart_hgnc_name = "NULL";
            }

            # Removed "$biomart_hgnc_name eq $cosmic_hgnc_name || " from comparison below
            # Now only comparing ENSTs
            
            if ($biomart_enst eq $cosmic_enst) {
                if ($biomart_trembl) {
                    push(@trembl_accs, $biomart_trembl);
                }
                if ($biomart_sprot) {
                    push(@sprot_accs, $biomart_sprot);
                }
            }
        }
        
        # Look for unique mappings
        
        # Add multiply-mapped or unmapped cosmic genes to %$unmapped_cosmic_ids for the next round of
        # mapping
        
        my @uniq_sprot_accs = uniq @sprot_accs;
        my @uniq_trembl_accs = uniq @trembl_accs;
        
        if (@uniq_sprot_accs == 0) {            # If we have no SPROT hits, consider TREMBL
            if (@uniq_trembl_accs == 0) {       # No TREMBL HITS - no hits at all
                $missing_c++;
                $unmapped_cosmic_ids{$cosmic_hgnc_name} = $cosmic_id;
            } elsif (@uniq_trembl_accs > 1) {   # Multiple TREMBL hits are a bad thing
                $multiple_c++;
                $unmapped_cosmic_ids{$cosmic_hgnc_name} = $cosmic_id;
            } else {                            # Unique TREMBL hit will be used in absence of SPROT
                $dbq_uniprot_id->execute($uniq_trembl_accs[0]);
                ( my $uniprot_id ) = $dbq_uniprot_id->fetchrow_array;
                $dbq_uniprot_id->finish;
                if (!$uniprot_id) {
                    print STDERR "\n--- Missing TREMBL/uniprot id for " . $uniq_trembl_accs[0] . "\n";
                    $tick_c = 1;
                    $missing_c++;
                } else {
                    $dbq_make_mapping->execute($uniprot_id, $cosmic_id);
                    $matched_c++;
                }
            }
        } elsif (@uniq_sprot_accs > 1) {        # Multpile SPROT hits are a bad thing
            $multiple_c++;
            $unmapped_cosmic_ids{$cosmic_hgnc_name} = $cosmic_id;
        } else {                                # Unique SPROT hit will be used
            $dbq_uniprot_id->execute($uniq_sprot_accs[0]);
            ( my $uniprot_id ) = $dbq_uniprot_id->fetchrow_array;
            $dbq_uniprot_id->finish;
            if (!$uniprot_id) {
                print STDERR "\n--- Missing SPROT/uniprot id for " . $uniq_sprot_accs[0] . "\n";
                $tick_c = 1;
                
                # Our unique SPROT hit failed.  Do we have a unique TREMBL hit?
                
                if (@uniq_trembl_accs == 1) {
                    $dbq_uniprot_id->execute($uniq_trembl_accs[0]);
                    ( my $uniprot_id ) = $dbq_uniprot_id->fetchrow_array;
                    $dbq_uniprot_id->finish;
                    if (!$uniprot_id) {
                        print STDERR "\n--- Missing TREMBL/uniprot id for " . $uniq_trembl_accs[0] . " (discovered in fall-back)\n";
                        $tick_c = 1;
                    } else {
                        $dbq_make_mapping->execute($uniprot_id, $cosmic_id);
                        $matched_c++;
                    }
                } else {
                    $missing_c++;
                }
            } else {
                $dbq_make_mapping->execute($uniprot_id, $cosmic_id);
                $matched_c++;
            }
        }
    }
    
    # Pointless ticker
    
    if ($tick_c % $tick_step == 0) {
        print STDERR $tick_char;
        if ($tick_c % ($tick_step * $col_wrap) == 0 && $tick_c) {
            print STDERR "\n";
        }
    }
    $tick_c++;
}
close($fasta_fh) || die "\n+++ Cannot close $cosmic_fasta_fn after reading data: $?\n";

$tick_c = 1;
print STDERR "\n--- $cosmic_c cosmic entries; $missing_c with no match, $multiple_c with multiple matches, $matched_c matched.\n";

# Second pass: now use the sequences

my $sequence_exact_c = 0;
$dbq_unmatched_entries->execute;        # For each unmatched cosmic entry
while (( $cosmic_id ) = $dbq_unmatched_entries->fetchrow_array) {
    my ( $cosmic_sequence, $cosmic_seq_md5, $uniprot_sequence, $uniprot_id, $match_c );
    $dbq_cosmic_sequence->execute($cosmic_id);      # Fetch cosmic sequence and md5 hash
    ( $cosmic_sequence, $cosmic_seq_md5 ) = $dbq_cosmic_sequence->fetchrow_array;
    $dbq_cosmic_sequence->finish;
    
    $dbq_uniprot_md5_match_c->execute($cosmic_seq_md5);
    ( $match_c ) = $dbq_uniprot_md5_match_c->fetchrow_array;
    $dbq_uniprot_md5_match_c->finish;
    
    if ($match_c > 1) {
        print STDERR "\n--- Multiple mappings for $cosmic_seq_md5\n";
        $tick_c = 1;
        
        # XXX With multiple mappings we need some sort of strategy.  SPROT first, etc. ???
    } elsif ($match_c == 1) {
        $dbq_uniprot_sequence_by_md5->execute($cosmic_seq_md5);
        ( $uniprot_id, $uniprot_sequence ) = $dbq_uniprot_sequence_by_md5->fetchrow_array;
        $dbq_uniprot_sequence_by_md5->finish;
        if ($uniprot_sequence eq $cosmic_sequence) {
            $dbq_make_mapping->execute($uniprot_id, $cosmic_id);
            $sequence_exact_c++;
        }
    }
}

$tick_c = 1;
print STDERR "\n--- Further $sequence_exact_c sequences matched by sequence identity.\n";

# Tick finish

print STDERR "\n";

$dbh_titanic->disconnect;
