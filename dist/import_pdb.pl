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

print STDERR "::: Import PDB - Import PDB chain sequences\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});
my $dbq_insert_pdb = $dbh_titanic->prepare(q{INSERT INTO pdb VALUES(NULL, ?, ?)});
my $dbq_pdb_id = $dbh_titanic->prepare(q{SELECT pdb_id FROM pdb WHERE pdb_code = ? AND pdb_chain = ?});
my $dbq_insert_sequence = $dbh_titanic->prepare(q{INSERT INTO sequence VALUES(NULL, 'pdb', ?, NULL, 'protein', ?, ?)});

# Prepare queries

# Go

my $pdb_seqres_fn = $mokcatanic_local_root . "pdb_seqres.txt.gz";

# pdb_seqres.txt is a typical FASTA format sequence file with sequences on single lines.

my ($pdb_code, $pdb_chain, $pdb_mol, $pdb_length, $pdb_name, $pdb_sequence);
open(PDB, "gunzip -c $pdb_seqres_fn |") || die "\n+++ Cannot open $pdb_seqres_fn to read sequence: $!\n";
$pdb_code = "NO_CODE";
while (my $line = <PDB>) {
    
    chomp($line);
    
    # Parse FASTA entry header
    
    if ($line =~ />(....)_(.) mol:(\w+) length:(\d+)  (.*)/) {
        $pdb_code = $1;
        $pdb_chain = $2;
        $pdb_mol = $3;
        $pdb_length = $4;
        $pdb_name = $5;
    } elsif ($pdb_code ne "NO_CODE") {
        $pdb_sequence = $line;
        
        # The file has DNA & RNA sequences for bound nucleotides at the end.
        # We only want to get protein sequences
        
        if ($pdb_mol eq "protein") {
            $dbq_insert_pdb->execute($pdb_code, $pdb_chain);
            my $pdb_id;
            $dbq_pdb_id->execute($pdb_code, $pdb_chain);
            ( $pdb_id ) = $dbq_pdb_id->fetchrow_array;
            $dbq_pdb_id->finish;
            
            my $sequence_md5 = md5_hex($pdb_sequence);
            $dbq_insert_sequence->execute($pdb_id, $pdb_sequence, $sequence_md5);
        }
    } else {
        die "Unexpected input: $line\n";
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

# Tick finish

print STDERR "\n";

$dbh_titanic->disconnect;
