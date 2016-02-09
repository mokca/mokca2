#!/usr/bin/perl -w

use Getopt::Long;
use FileHandle;
use DBI();
use Digest::MD5 qw(md5 md5_hex);

use database_setup;
use path_setup;

# Autoflush STDERR

STDERR->autoflush(1);

# Ticker

my $tick_c = 0;
my $col_wrap = 78;
my $tick_step = 10000;
my $tick_char = '.';

# Explain purpose

print STDERR "::: Import Uniprot sequences - Import sequences for all uniprot entries, including isoforms\n";

# Setup universal paths, database details

my $sprot_seq_fn = $mokcatanic_data_root . "uniprot_sprot.fasta.gz";
my $trembl_seq_fn = $mokcatanic_local_root . "uniprot_trembl.fasta.gz";
my $uniprot_iso_fn = $mokcatanic_data_root . "uniprot_sprot_varsplic.fasta.gz";

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_accs = $dbh_titanic->prepare(q{SELECT DISTINCT accession FROM uniprot});
my $dbq_uniprot_id = $dbh_titanic->prepare(q{SELECT uniprot_id FROM uniprot WHERE accession = ?});
my $dbq_insert_seq = $dbh_titanic->prepare(q{INSERT INTO sequence VALUES (NULL, 'uniprot', ?, ?, 'protein', ?, ?)});

# Go

my %accessions;
my ( $uniprot_accession, $uniprot_isoform, $uniprot_id, $sequence_md5, $sequence_c );
my $sequence;

$sequence = "";
$uniprot_accession = "";
$sequence_c = 0;

open my $uniprot_fh, "gunzip -c $sprot_seq_fn |" || die "\n+++ Cannot open $sprot_seq_fn to read data: $?\n";
while ($line = <$uniprot_fh>) {
    chomp $line;
    if ($line =~ /^>..\|(.{6,10})\|/) {

        if ($uniprot_accession) {
            
            $dbq_uniprot_id->execute($uniprot_accession);
            ( $uniprot_id ) = $dbq_uniprot_id->fetchrow_array;
            $dbq_uniprot_id->finish;
            if (!$uniprot_id) {
                
                # Sequence files are for all ogranisms, so most sequences won't be in our uniprot
                # table of human proteins
                
            } else {
                $sequence_md5 = md5_hex($sequence);
                $dbq_insert_seq->execute($uniprot_id, $uniprot_isoform, $sequence, $sequence_md5);
                $sequence_c++;
            }
            
            $sequence = "";
        }
        
        $uniprot_accession = $1;
        $uniprot_isoform = 1;
        ++$accessions{$uniprot_accession};
        
        # Pointless ticker
        
        if ($tick_c % $tick_step == 0) {
            print STDERR $tick_char;
            if ($tick_c % ($tick_step * $col_wrap) == 0 && $tick_c) {
                print STDERR "\n";
            }
        }
        $tick_c++;
    } else {
        $sequence .= $line;
    }
}
close $uniprot_fh || die "\n+++ Cannot close $sprot_seq_fn to stop reading data: $?\n";

print STDERR "\n--- $sequence_c sequences imported from $sprot_seq_fn\n";
$tick_c = 0;

$sequence = "";
$uniprot_accession = "";
$uniprot_isoform = 0;
$sequence_c = 0;

open $uniprot_fh, "gunzip -c $uniprot_iso_fn |" || die "\n+++ Cannot open $uniprot_iso_fn to read data: $?\n";
while ($line = <$uniprot_fh>) {
    chomp $line;
    if ($line =~ /^>..\|(.{6,10})-(\d+)\|/) {
        
        if ($uniprot_accession) {
            
            $dbq_uniprot_id->execute($uniprot_accession);
            ( $uniprot_id ) = $dbq_uniprot_id->fetchrow_array;
            $dbq_uniprot_id->finish;
            if (!$uniprot_id) {
                
                # Sequence files are for all ogranisms, so most sequences won't be in our uniprot
                # table of human proteins
                
            } else {
                $sequence_md5 = md5_hex($sequence);
                $dbq_insert_seq->execute($uniprot_id, $uniprot_isoform, $sequence, $sequence_md5);
                $sequence_c++;
            }
            
            $sequence = "";
        }
        
        $uniprot_accession = $1;
        $uniprot_isoform = $2;
        ++$accessions{$uniprot_accession};
        
        # Pointless ticker
        
        if ($tick_c % $tick_step == 0) {
            print STDERR $tick_char;
            if ($tick_c % ($tick_step * $col_wrap) == 0 && $tick_c) {
                print STDERR "\n";
            }
        }
        $tick_c++;
    } else {
        $sequence .= $line;
    }
}
close $uniprot_fh || die "\n+++ Cannot close $uniprot_iso_fn to stop reading data: $?\n";

print STDERR "\n--- $sequence_c sequences imported from $uniprot_iso_fn\n";
$tick_c = 0;

$sequence = "";
$uniprot_accession = "";
$sequence_c = 0;

open $uniprot_fh, "gunzip -c $trembl_seq_fn |" || die "\n+++ Cannot open $trembl_seq_fn to read data: $?\n";
while ($line = <$uniprot_fh>) {
    chomp $line;
    if ($line =~ /^>..\|(.{6,10})\|/) {
        
        if ($uniprot_accession) {
            
            $dbq_uniprot_id->execute($uniprot_accession);
            ( $uniprot_id ) = $dbq_uniprot_id->fetchrow_array;
            $dbq_uniprot_id->finish;
            if (!$uniprot_id) {
                
                # Sequence files are for all ogranisms, so most sequences won't be in our uniprot
                # table of human proteins
                
            } else {
                $sequence_md5 = md5_hex($sequence);
                $dbq_insert_seq->execute($uniprot_id, $uniprot_isoform, $sequence, $sequence_md5);
                $sequence_c++;
            }
            
            $sequence = "";
        }
        
        $uniprot_accession = $1;
        $uniprot_isoform = 1;
        ++$accessions{$uniprot_accession};
    
        # Pointless ticker
        
        if ($tick_c % $tick_step == 0) {
            print STDERR $tick_char;
            if ($tick_c % ($tick_step * $col_wrap) == 0 && $tick_c) {
                print STDERR "\n";
            }
        }
        $tick_c++;
    } else {
        $sequence .= $line;
    }
}
close $uniprot_fh || die "\n+++ Cannot close $trembl_seq_fn to stop reading data: $?\n";

print STDERR "\n--- $sequence_c sequences imported from $trembl_seq_fn\n";
$tick_c = 0;

$dbq_accs->execute();
while ((my $acc) = $dbq_accs->fetchrow_array) {
    if (!$accessions{$acc}) {
        print "--- Missing sequence for $acc\n";
    }
}
