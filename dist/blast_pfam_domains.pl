#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use FileHandle;
use DBI();
use Parallel::ForkManager;
use XML::LibXML;

# Autoflush STDERR

STDERR->autoflush(1);

# Ticker

my $tick_c = 0;
my $col_wrap = 78;
my $tick_step = 10;
my $tick_char = '.';

# Explain purpose

print STDERR "::: BLAST mutated pFam domains - \n";

# Setup universal paths, database details

use path_setup;
use database_setup;

my $pfilt_fn = $mokcatanic_local_root . "uniref_pdb_pfilt.txt";

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_domains = $dbh_titanic->prepare(q{SELECT uniprot_domain_id, uniprot_id, uniprot_start, uniprot_end FROM uniprot_domain WHERE domain_type = 'pfam_A' AND mutated = 1});

# Go

my $ARB_IDENT_CUTOFF = 30.0;

my $FORKTHREADS = 10;
my $pm = new Parallel::ForkManager($FORKTHREADS);

my ($uniprot_domain_id, $uniprot_id, $uniprot_start, $uniprot_end);
$dbq_domains->execute;
while (($uniprot_domain_id, $uniprot_id, $uniprot_start, $uniprot_end) = $dbq_domains->fetchrow_array) {
    
    my $thread_pid = $pm->start and next;
    
    blast_in_a_basket($uniprot_domain_id, $uniprot_id, $uniprot_start, $uniprot_end);
    
    if ($tick_c % $tick_step == 0) {
        print STDERR $tick_char;
        if ($tick_c % ($tick_step * $col_wrap) == 0 && $tick_c) {
            print STDERR "\n";
        }
    }
    $tick_c++;
    
    $pm->finish;
}

# Tick finish

print STDERR "\n";

$dbh_titanic->disconnect;

sub blast_in_a_basket {
    
    # Connect to databases
    
    my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});
    
    # Prepare queries
    
    my $dbq_map_exists = $dbh_titanic->prepare(q{SELECT COUNT(*) FROM alignment_mapping WHERE query_db = 'domain' AND query_id = ?});
    my $dbq_domain_aligned = $dbh_titanic->prepare(q{SELECT COUNT(query_pos) FROM alignment_mapping WHERE query_db = 'domain' AND hit_db = 'pdb' AND query_id = ?});
    my $dbq_uniprot_seq = $dbh_titanic->prepare(q{SELECT sequence FROM sequence WHERE source = 'uniprot' AND source_id = ?});
    
    # Initialise variables, output
    
    my ($uniprot_domain_id, $uniprot_id, $uniprot_start, $uniprot_end);
    ( $uniprot_domain_id, $uniprot_id, $uniprot_start, $uniprot_end ) = @_;
    my $o_buff = "*---\n";
    $o_buff .= "--- Domain $uniprot_domain_id in $uniprot_id from $uniprot_start to $uniprot_end\n";
    
    # Check to see if this has already been done
    
    my $alignment_c = 0;
    $dbq_domain_aligned->execute($uniprot_domain_id);
    ( $alignment_c ) = $dbq_domain_aligned->fetchrow_array;
    $dbq_domain_aligned->finish;
    
    if ($alignment_c > 0) {
        $o_buff .= "SELECT COUNT(query_pos) FROM alignment_mapping WHERE query_db = 'domain' 
        AND hit_db = 'pdb' AND query_id = $uniprot_domain_id\n";
        $o_buff .= "--- Already aligned ($alignment_c)\n";
        
        print $o_buff;
        return;
    }
    
    # Get uniprot sequence & thence domain sequence for this domain
    
    my ($uniprot_seq, $domain_seq);
    $dbq_uniprot_seq->execute($uniprot_id);
    ( $uniprot_seq ) = $dbq_uniprot_seq->fetchrow_array;
    $dbq_uniprot_seq->finish;
    
    $domain_seq = substr($uniprot_seq, $uniprot_start - 1, $uniprot_end - $uniprot_start + 1);
    
    # Run PSIBLAST
    
    my $query_fn = $tmp_root . $uniprot_domain_id . "_query.fasta";
    my $result_fn = $tmp_root . $uniprot_domain_id . "_result.xml";
    
    open (QUERY, "> $query_fn") || die "\n++++ Cannot open $query_fn to write sequence: $!\n";
    print QUERY ">$uniprot_domain_id/$uniprot_start-$uniprot_end\n$domain_seq\n";
    close (QUERY) || die "\n+++ Cannot close $query_fn: $!\n";
    
    my $psiblast_sys = "psiblast -num_iterations 1 -evalue 0.001 -outfmt 5 -db $pfilt_fn -query $query_fn -out $result_fn";
    
    system($psiblast_sys) == 0 || die "\n+++ $psiblast_sys failed: $!\n";
    
    # Parse PSIBLAST XML output
    
    my $parser = XML::LibXML->new(load_ext_dtd => 0);
    my $tree = $parser->parse_file($result_fn);
    my $root = $tree->getDocumentElement;
    
    my $query_len = $root->findvalue('BlastOutput_query-len');
    $o_buff .= "---    Query length: $query_len\n";
    foreach my $BlastOutput_iterations ($root->findnodes('BlastOutput_iterations')) {
        my $iteration_max = 0;
        
        # Count the number of iterations performed
        
        foreach my $iteration ($BlastOutput_iterations->findnodes('Iteration')) {
            my $iteration_c = $iteration->findvalue('Iteration_iter-num');
            $o_buff .= "---    Iteration $iteration_c\n";
            if ($iteration_c > $iteration_max) {
                $iteration_max = $iteration_c;
            }
        }
        $o_buff .=  "---    Max iteration: $iteration_max\n";
        
        # We're only interested in results from the last iteration, so skip to that one
        
        foreach my $iteration ($BlastOutput_iterations->findnodes('Iteration')) {
            my $iteration_c = $iteration->findvalue('Iteration_iter-num');
            if ($iteration_c == $iteration_max) {
                $o_buff .= "---    Iteration $iteration_c\n";
                foreach my $Iteration_hits ($iteration->findnodes('Iteration_hits')) {
                    my $hit_c = 0;
                    my $PDB_hit_c = 0;
                    
                    # Look at all hits in the last iteration, picking hits on PDB sequences
                    
                    foreach my $Hit ($Iteration_hits->findnodes('Hit')) {
                        $hit_c++;
                        my $hit_def = $Hit->findvalue('Hit_def');
                        if ($hit_def =~ /UniRef90/) {
                            
                            # Ignore hits on UniRef sequences
                        
                        } else {
                            
                            # This is a hit on a PDB sequence
                            
                            $o_buff .= "---       PDB Hit: $hit_def\n";
                            $PDB_hit_c++;
                            
                            $hit_def =~ /(....)_(.)/;
                            my $pdb_code = $1;
                            my $pdb_chain = $2;
                            
                            $o_buff .= "---       PDB Hit: $pdb_code $pdb_chain\n";
                            
                            my $hit_hsps = $Hit->findnodes('Hit_hsps')->get_node(0);
                            my $hsp_c = 0;
                            foreach my $Hsp ($hit_hsps->findnodes('Hsp')) {
                                my $query_from = $Hsp->findvalue('Hsp_query-from');
                                my $query_to = $Hsp->findvalue('Hsp_query-to');
                                my $hit_from = $Hsp->findvalue('Hsp_hit-from');
                                my $hit_to = $Hsp->findvalue('Hsp_hit-to');
                                my $e_value = $Hsp->findvalue('Hsp_evalue');
                                my $hsp_query_seq = $Hsp->findvalue('Hsp_qseq');
                                my $hsp_hit_seq = $Hsp->findvalue('Hsp_hseq');
                                my $hsp_align_len = $Hsp->findvalue('Hsp_align-len');
                                my $hsp_identity = $Hsp->findvalue('Hsp_identity');
                                my $hsp_positive = $Hsp->findvalue('Hsp_positive');
                                
                                my $query_aligned_len = $query_to - $query_from + 1;
                                my $percent_ident = 100.0 * $hsp_identity / $query_aligned_len;
                                my $percent_positive = 100.0 * $hsp_positive / $query_aligned_len;
                                
                                # Quick sanity check
                                
                                if (length($hsp_query_seq) != $hsp_align_len || length($hsp_hit_seq) != $hsp_align_len) {
                                    $o_buff .= "+++       XXX Length disparity: hsp length $hsp_align_len, query " . length($hsp_query_seq) . " vs. hit " . length($hsp_hit_seq) . "\n";
                                } elsif ($percent_ident < $ARB_IDENT_CUTOFF) {
                                    $o_buff .= "---       Percentage identity below cutoff: $percent_ident\n";
                                } else {
                                 
                                    # Passed sanity checks
                                }
                                $hsp_c++;
                            }
                            
                            if ($hsp_c != 1) {
                                $o_buff .= " +++        XXX HSP Count $hsp_c\n";
                            }
                        }
                    }
                }
            }
        }
    }
    
    # Tidy up
    
    $dbh_titanic->disconnect;
    
    $o_buff .= "*---\n";
    print $o_buff;
}

