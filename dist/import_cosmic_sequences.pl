#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use FileHandle;
use DBI();
use Digest::MD5 qw(md5 md5_hex);

# Autoflush STDERR

STDERR->autoflush(1);

# Ticker

my $tick_c = 0;
my $col_wrap = 78;
my $tick_step = 10;
my $tick_char = '.';

# Explain purpose

print STDERR "::: Import Cosmic Sequences - import gene and protein sequences from Cosmic\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

my $cosmic_fasta_fn = $mokcatanic_data_root . "All_COSMIC_Genes.fasta.gz";

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_cosmic_id = $dbh_titanic->prepare(q{SELECT cosmic_id FROM cosmic_hgnc WHERE cosmic_gene_name = ?});
my $dbq_insert_gene = $dbh_titanic->prepare(q{INSERT INTO sequence VALUES(NULL,'cosmic',?,1,'gene',?,?)});
my $dbq_insert_protein = $dbh_titanic->prepare(q{INSERT INTO sequence VALUES(NULL,'cosmic',?,1,'protein',?,?)});

# Go

# Translation matrix

my %aacode = (
TTT => "F", TTC => "F", TTA => "L", TTG => "L",
TCT => "S", TCC => "S", TCA => "S", TCG => "S",
TAT => "Y", TAC => "Y", TAA => "", TAG => "",
TGT => "C", TGC => "C", TGA => "", TGG => "W",
CTT => "L", CTC => "L", CTA => "L", CTG => "L",
CCT => "P", CCC => "P", CCA => "P", CCG => "P",
CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
CGT => "R", CGC => "R", CGA => "R", CGG => "R",
ATT => "I", ATC => "I", ATA => "I", ATG => "M",
ACT => "T", ACC => "T", ACA => "T", ACG => "T",
AAT => "N", AAC => "N", AAA => "K", AAG => "K",
AGT => "S", AGC => "S", AGA => "R", AGG => "R",
GTT => "V", GTC => "V", GTA => "V", GTG => "V",
GCT => "A", GCC => "A", GCA => "A", GCG => "A",
GAT => "D", GAC => "D", GAA => "E", GAG => "E",
GGT => "G", GGC => "G", GGA => "G", GGG => "G",
);

# For all sequences in the COSMIC fasta gene dump

my ( $cosmic_gene_seq, $cosmic_protein_seq, $cosmic_hgnc_name, $cosmic_enst, $cosmic_id, $sequence_md5 );
$cosmic_gene_seq = "";
open my $fasta_fh, "gunzip -c $cosmic_fasta_fn |" || die "\n+++ Cannot open $cosmic_fasta_fn to read data: $?\n";
while (my $line = <$fasta_fh>) {
    chomp($line);
    if ($line =~ /^>([\w\-.]+)\s([\w\-.]+)/) {  # FASTA information line
        
        if ($cosmic_gene_seq) {                 # Start of one sequence implies the end of the previous one

            # Translate the previous sequence
            my @codons = unpack '(A3)*', $cosmic_gene_seq;
            my @aminoAcids = map { exists $aacode{$_} ? $aacode{$_} : "X" } @codons;
            $cosmic_protein_seq = join '', @aminoAcids;
        
            #Fetch the COSMIC id
            $dbq_cosmic_id->execute($cosmic_hgnc_name);
            ( $cosmic_id ) = $dbq_cosmic_id->fetchrow_array;
            $dbq_cosmic_id->finish;
            if (! $cosmic_id) {
                die "\n+++ Cannot fetch ID for $cosmic_hgnc_name\n";
            }
            
            # Insert gene and protein into database
            $dbq_insert_gene->execute($cosmic_id, $cosmic_gene_seq, md5_hex($cosmic_gene_seq));
            $dbq_insert_protein->execute($cosmic_id, $cosmic_protein_seq, md5_hex($cosmic_protein_seq));
        }
        
        # Get new FASTA information
        $cosmic_hgnc_name = $1;
        $cosmic_enst = $2;
        $cosmic_gene_seq = "";
        
        # Pointless ticker
        if ($tick_c % $tick_step == 0) {
            print STDERR $tick_char;
            if ($tick_c % ($tick_step * $col_wrap) == 0 && $tick_c) {
                print STDERR "\n";
            }
        }
        $tick_c++;
        
    } else {                                    # Sequence or continuation thereof
        $cosmic_gene_seq .= uc $line;
    }
}

# Tick finish

print STDERR "\n";

$dbh_titanic->disconnect;
