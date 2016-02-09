#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use FileHandle;
use DBI();
use XML::LibXML;
use Parallel::ForkManager;

# Autoflush STDERR

STDERR->autoflush(1);

# Ticker

my $tick_c = 0;
my $col_wrap = 78;
my $tick_step = 10;
my $tick_char = '.';

# Explain purpose

print STDERR "::: BLAST COSMIC to UNIPROT - Match remaining sequences by BLASTing against UNIPROT human sequences\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

my $blast_db_path = $mokcatanic_local_root . "uniprot_human.seq";

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_unmapped = $dbh_titanic->prepare(q{SELECT cosmic_id FROM cosmic_vs_uniprot WHERE uniprot_id IS NULL});

# Go

my $fm = new Parallel::ForkManager(10);

my $unmapped_cosmic_id;
$dbq_unmapped->execute();
while (($unmapped_cosmic_id) = $dbq_unmapped->fetchrow_array) {
    
    $fm->start and next;
    
    blast_in_a_basket($unmapped_cosmic_id);
    
    $fm->finish;
    
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

sub blast_in_a_basket {
    
    # Parse arguments
    
    my $cosmic_id;
    ( $cosmic_id ) = @_;
    
    # Connect to databases
    
    my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});
    
    # Prepare queries

    my $dbq_sequence = $dbh_titanic->prepare(q{SELECT sequence FROM sequence WHERE source_id = ? AND source = 'cosmic' AND codes_for = 'protein'});
    my $dbq_accession = $dbh_titanic->prepare(q{SELECT uniprot_id FROM uniprot_accession WHERE accession = ?});
    my $dbq_uniprot_sequence = $dbh_titanic->prepare(q{SELECT sequence FROM sequence WHERE source_id = ? AND source = 'uniprot' AND codes_for = 'protein' AND isoform = 1});
    my $dbq_make_mapping = $dbh_titanic->prepare(q{UPDATE cosmic_vs_uniprot SET uniprot_id = ?, identity = ? WHERE cosmic_id = ?});
    my $dbq_commit_alignment = $dbh_titanic->prepare(q{INSERT INTO alignment_mapping VALUES('cosmic', 'uniprot', ?, ?, ?, ?, ?, ?, ?)});
    
    # Fetch sequence, BLAST sequence
    
    my ( $sequence, $query_fn, $result_fn, $uniprot_fn, $align_fn, $query_fh, $result_fh, $uniprot_fh, $align_fh );
    $query_fn = "/tmp/$cosmic_id.seq";
    $result_fn = "/tmp/$cosmic_id.xml";
    
    $dbq_sequence->execute($cosmic_id);
    ( $sequence ) = $dbq_sequence->fetchrow_array;
    $dbq_sequence->finish;
    
    open ($query_fh, ">$query_fn") || die "\n+++ Cannot open $query_fn to write sequence: $?\n";
    print $query_fh ">$cosmic_id\n$sequence\n";
    close ($query_fh) || die "\n+++ Cannot close $query_fn: $?\n";
    
    my $blast_command = "psiblast -num_iterations 1 -evalue 0.001 -outfmt 5 -db $blast_db_path -query $query_fn -out $result_fn";
    system($blast_command);
    
    # Parse XML blast results
    
    my ( $parser, $tree, $root, $BlastOutput_iterations, $Iteration, $Iteration_hits, $Hit, $Hit_def );
    $parser = XML::LibXML->new(load_ext_dtd => 0);
    $tree = $parser->parse_file($result_fn);
    $root = $tree->getDocumentElement;
    $BlastOutput_iterations= $root->findnodes('BlastOutput_iterations')->get_node(1);
    $Iteration = $BlastOutput_iterations->findnodes('Iteration')->get_node(1);
    $Iteration_hits = $Iteration->findnodes('Iteration_hits')->get_node(1);
    $Hit = $Iteration_hits->findnodes('Hit')->get_node(1);
    if ($Hit) { # It's possible to not get *any* hits - the data is awful
        my ( $uniprot_id, $uniprot_sequence, $uniprot_accession);
        $uniprot_accession = $Hit->findvalue('Hit_def');
        $dbq_accession->execute($uniprot_accession);
        ($uniprot_id) = $dbq_accession->fetchrow_array;
        $dbq_accession->finish;
        if (! $uniprot_id ) {
            die "\n+++ Really?  No id for $uniprot_accession\n";
        }
        $dbq_uniprot_sequence->execute($uniprot_id);
        ( $uniprot_sequence ) = $dbq_uniprot_sequence->fetchrow_array;
        $dbq_uniprot_sequence->finish;
        $uniprot_fn = "/tmp/$uniprot_accession.seq";
        $align_fn = "/tmp/$cosmic_id.aln";
        open($uniprot_fh, ">$uniprot_fn") || die "\n+++ Cannot open $uniprot_fn to write sequence: $?\n";
        print $uniprot_fh ">$uniprot_accession\n$uniprot_sequence\n";
        close($uniprot_fh) || die "\n+++ Cannot close $uniprot_fn: $?\n";
        
        # Do a further alignment with needle rather than trying to piece it together from (possible)
        # multiple HSPs
        
        if (system("needle", $query_fn, $uniprot_fn, "-auto", "-aformat", "fasta", "-gapopen", "10.0", "-gapextend", "0.5", "-datafile", "EBLOSUM62", "-outfile", $align_fn) != 0) {
            die "\n+++ Failed to run needle on $cosmic_id vs $uniprot_accession\n";
        }
        my ($line, $cosmic_alignment, $uniprot_alignment, $alignment_c);
        $alignment_c = 0;
        $cosmic_alignment = "";
        $uniprot_alignment = "";
        open($align_fh, "<$align_fn") || die "\n+++ Cannot open $align_fn to read sequence alignment: $?\n";
        while ($line = <$align_fh>) { # Read alignment
            chomp($line);
            if ($line =~  />\w+/) {
                $alignment_c++;
            } else {
                if ($alignment_c == 1) {
                    $cosmic_alignment .= $line;
                } else {
                    $uniprot_alignment .= $line;
                }
            }
        }
        close($align_fh) || die "\n+++ Cannot close $align_fn: $?\n";
        if (length($cosmic_alignment) != length($uniprot_alignment)) {
            die "\n+++ Alignment length mismatch $cosmic_alignment $uniprot_alignment\n";
        }
        my ($cosmic_res, $uniprot_res, $cosmic_i, $uniprot_i, $cosmic_gap, $uniprot_gap, $gap_location, $identity_c, $identity);
        $cosmic_i = 0;
        $uniprot_i = 0;
        $identity_c = 0;
        for (my $ii=0;$ii<length($cosmic_alignment);$ii++) {
            $cosmic_res = substr($cosmic_alignment, $ii, 1);
            $uniprot_res = substr($uniprot_alignment, $ii, 1);
            if ($cosmic_res eq $uniprot_res) {
                $identity_c++;
            }
            $gap_location = "none";
            if ($cosmic_res eq '-') {
                $cosmic_gap = 1;
                $gap_location = "query";
            } else {
                $cosmic_gap = 0;
                $cosmic_i++;
            }
            if ($uniprot_res eq '-') {
                $uniprot_gap = 1;
                $gap_location = "hit";
            } else {
                $uniprot_gap = 0;
                $uniprot_i++;
            }
            $dbq_commit_alignment->execute($cosmic_id, $uniprot_id, $cosmic_i, $uniprot_i, $cosmic_res, $uniprot_res, $gap_location);
            $dbq_commit_alignment->finish;
        }
        
        # As the data is so shocking, we need an artibtrary %ge identity cutoff
        
        $identity = 100.0 * $identity_c / length($sequence);
        if ($identity > 75.0) {
            $dbq_make_mapping->execute($uniprot_id, $identity, $cosmic_id);
            $dbq_make_mapping->finish;
        }
    }
    
    # Tidy up
    
    $dbh_titanic->disconnect;
}