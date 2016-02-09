#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use FileHandle;
use DBI();

# Autoflush STDERR

STDERR->autoflush(1);

# Ticker

my $tick_c = 0;
my $col_wrap = 78;
my $tick_step = 10;
my $tick_char = '.';

# Explain purpose

print STDERR "::: Enumerate Domains - create a list of pfamA and pfamB domains and inter-domain regions\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});
my $dbh_pfam = DBI->connect("DBI:mysql:database=pfam_27_0;host=apsu","mokca","iwt4lp", {'RaiseError' => 1});

# Prepare queries

my $dbq_uniprot_id = $dbh_titanic->prepare(q{SELECT uniprot_id FROM uniprot});
my $dbq_cosmic_id = $dbh_titanic->prepare(q{SELECT cosmic_id, identity FROM cosmic_vs_uniprot WHERE uniprot_id = ?});
my $dbq_uniprot_acc = $dbh_titanic->prepare(q{SELECT accession FROM uniprot WHERE uniprot_id = ?});
my $dbq_uniprot_seq = $dbh_titanic->prepare(q{SELECT sequence FROM sequence WHERE source = 'uniprot' AND source_id = ? AND isoform = 1 AND codes_for = 'protein'});
my $dbq_pfamA = $dbh_pfam->prepare(q{SELECT pfamA_acc, pfamA_id, seq_start, seq_end FROM pfamseq, pfamA, pfamA_reg_full_significant
    WHERE pfamseq_acc = ?
    AND in_full = 1 AND pfamseq.auto_pfamseq = pfamA_reg_full_significant.auto_pfamseq AND pfamA_reg_full_significant.auto_pfamA = pfamA.auto_pfamA
    ORDER BY seq_start ASC});
my $dbq_pfamB = $dbh_pfam->prepare(q{SELECT DISTINCT pfamB.pfamB_acc, pfamB_id, seq_start, seq_end
    FROM   pfamB_reg, pfamB, pfamseq
    WHERE  pfamseq_acc = ?
    AND    pfamB_reg.auto_pfamseq = pfamseq.auto_pfamseq AND pfamB_reg.auto_pfamB = pfamB.auto_pfamB
    ORDER BY seq_start ASC});

# Go

my $uniprot_id;
$dbq_uniprot_id->execute;
while (($uniprot_id) = $dbq_uniprot_id->fetchrow_array) { # For each uniprot ID...
 
    # ... fetch primary accession code
    
    # ( XXX And ignore the fact that there are alternative accession codes, some of which have
    #       the same pfam domains, some of which have different domains, and some of which have
    #        none. XXX )
    
    my $uniprot_acc;
    $dbq_uniprot_acc->execute($uniprot_id);
    ($uniprot_acc) = $dbq_uniprot_acc->fetchrow_array;
    $dbq_uniprot_acc->finish;
    
    # ... fetch cosmic ID
    
    my ($cosmic_id, $identity);
    $dbq_cosmic_id->execute($uniprot_id);
    ( $cosmic_id, $identity ) = $dbq_cosmic_id->fetchrow_array;
    $dbq_cosmic_id->finish;
    
    # ... fetch sequence for isoform 1
    
    my ($uniprot_seq, $uniprot_len);
    $dbq_uniprot_seq->execute($uniprot_id);
    ( $uniprot_seq ) = $dbq_uniprot_seq->fetchrow_array;
    $dbq_uniprot_seq->finish;
    
    $uniprot_len = length($uniprot_seq);
    if (!$uniprot_len) {
        print STDERR "\n+++ No sequence found for $uniprot_acc\n";
        $tick_c = 0;
    }
    
    # Create map of uniprot sequence, and remove pfamA+B domains to leave inter-domain regions
    
    my @uniprot_map;
    $uniprot_map[0] = 0;
    for (my $ii=1;$ii<=$uniprot_len;$ii++) {
        $uniprot_map[$ii] = 1;
    }
    
    # Fetch pfamA domains associated with the accession code
    
    my ( $pfam_acc, $pfam_id, $seq_start, $seq_end );
    $dbq_pfamA->execute($uniprot_acc);
    while (( $pfam_acc, $pfam_id, $seq_start, $seq_end ) = $dbq_pfamA->fetchrow_array) {
        #print "::: Domain $seq_start -> $seq_end\n";

        for (my $ii=$seq_start;$ii<=$seq_end;$ii++) {
            $uniprot_map[$ii] = 0;
        }
        
        map_and_insert_domain("pfam_A", $pfam_acc, $pfam_id, $cosmic_id, $identity, $uniprot_id, $uniprot_acc, $seq_start, $seq_end);
    }
    $dbq_pfamA->finish;
    
    # Fetch pfamB domains associated with the accession code
    
    $dbq_pfamB->execute($uniprot_acc);
    while (( $pfam_acc, $pfam_id, $seq_start, $seq_end ) = $dbq_pfamB->fetchrow_array) {
        #print "::: Domain $seq_start -> $seq_end\n";
        
        # XXX At the moment we use just pfamA domains for domain/inter-domain residue assignments XXX
        
        #for (my $ii=$seq_start;$ii<=$seq_end;$ii++) {
        #    $uniprot_map[$ii] = 0;
        #}
        
        map_and_insert_domain("pfam_B", $pfam_acc, $pfam_id, $cosmic_id, $identity, $uniprot_id, $uniprot_acc, $seq_start, $seq_end);

    }
    $dbq_pfamB->finish;
    
    # Enumerate the inter-domain regions
    
    my $in_inter = 0;
    my $inter_start = 0;
    my $inter_end = 0;
    for (my $ii=0;$ii<=$uniprot_len;$ii++) {
        if ($uniprot_map[$ii] == 0 && $in_inter) {
            $inter_end = $ii-1;
            #print "---    inter domain $inter_start -> $inter_end\n";
            map_and_insert_domain("inter", "inter", "inter", $cosmic_id, $identity, $uniprot_id, $uniprot_acc, $inter_start, $inter_end);
            $in_inter = 0;
        } elsif ($uniprot_map[$ii] == 1 && !$in_inter) {
            $inter_start = $ii;
            $in_inter = 1;
        }
    }
    if ($in_inter) {
        $inter_end = $uniprot_len;
        #print "---    inter domain $inter_start -> $inter_end\n";
        
         map_and_insert_domain("inter", "inter", "inter", $cosmic_id, $identity, $uniprot_id, $uniprot_acc, $inter_start, $inter_end);
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

sub map_and_insert_domain {
    
    # Database connection
    
    my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

    # Prepare queries
    
    my $dbq_insert_domain = $dbh_titanic->prepare(q{INSERT INTO uniprot_domain VALUES (NULL,?,?,?,?,?,?,?,?,?)});
    my $dbq_domain_id = $dbh_titanic->prepare(q{SELECT uniprot_domain_id FROM uniprot_domain WHERE domain_type = ? AND domain_acc = ? AND domain_id = ? AND uniprot_id = ? AND uniprot_start = ? AND uniprot_end = ?});
    my $dbq_update_mutation = $dbh_titanic->prepare(q{UPDATE uniprot_domain SET mutated = ? WHERE uniprot_domain_id = ?});
    my $dbq_insert_agg = $dbh_titanic->prepare(q{INSERT INTO domain_map VALUES (?,?)});
    my $dbq_map = $dbh_titanic->prepare(q{SELECT query_pos, gap FROM alignment_mapping WHERE query_db = 'cosmic' AND hit_db = 'uniprot' AND query_id = ? AND hit_id = ? AND hit_pos = ?});
    my $dbq_map_exists = $dbh_titanic->prepare(q{SELECT COUNT(hit_pos) FROM alignment_mapping WHERE query_db = 'cosmic' AND hit_db = 'uniprot' AND query_id = ?});
    my $dbq_aggs = $dbh_titanic->prepare(q{SELECT aggregate_id FROM mut_aggregate WHERE cosmic_id = ? AND map_aa_start >= ? AND map_aa_start <= ? AND map_aa_stop >= ? AND map_aa_stop <= ?});
    
    # Get subroutine parameters
    
    my ($dom_type, $dom_acc, $dom_id, $cosmic_id, $identity, $uniprot_id, $uniprot_acc, $seq_start, $seq_end);
    ($dom_type, $dom_acc, $dom_id, $cosmic_id, $identity, $uniprot_id, $uniprot_acc, $seq_start, $seq_end)= @_;
    
    #print "$dom_type, $dom_acc, $dom_id, $cosmic_id, $identity, $uniprot_id, $uniprot_acc, $seq_start, $seq_end\n";
    
    # Map domain boundaries to cosmic_sequence (if this uniprot entry is also in COSMIC)
    
    # In a departure from the method previously used we start at the boundaries of the domain and
    # move inwards until there isn't a gap in the sequence alignment, giving up when we get to the
    # opposite boundary.  What could possibly go wrong.
    
    my ($seq_start_try, $seq_end_try);
    $seq_start_try = 0;
    $seq_end_try = 0;
    
    if ($cosmic_id) {
        if ($identity == 100.0) {
            $seq_start_try = $seq_start;
            $seq_end_try = $seq_end;
        } else {
            
            my $map_res_c;
            $dbq_map_exists->execute($cosmic_id);
            ( $map_res_c ) = $dbq_map_exists->fetchrow_array;
            $dbq_map_exists->finish;
            
            if ($map_res_c == 0) {
                #print "+++ No mappings for $cosmic_id/$uniprot_id ($uniprot_acc)\n";
            } else {
                
                #print "--- $uniprot_acc ($cosmic_id, $uniprot_id) : $seq_start to $seq_end, try ";
                
                
                $seq_start_try = $seq_start;
                
                my ($cosmic_pos, $gap);
                do {
                    #print "$seq_start_try ";
                    
                    $dbq_map->execute($cosmic_id, $uniprot_id, $seq_start_try);
                    ( $cosmic_pos, $gap ) = $dbq_map->fetchrow_array;
                    $dbq_map->finish;
                    
                    if (!$gap) {    # We're looking after the end of both sequences
                        $gap = "both";
                        print "--- $uniprot_acc ($cosmic_id, $uniprot_id) : $seq_start to $seq_end\n";
                    }
                    
                    if ($gap ne "none") {
                        $seq_start_try++;
                    }
                } until ($seq_start_try > $seq_end || $gap eq "none");
                if ($gap eq "none") {
                    #print "...using $seq_start_try\n";
                } else {
                    #print "...given up at $seq_start_try\n";
                    $seq_start_try = 0;
                }
                
                #print "    and try ";
                $seq_end_try = $seq_end;
                
                do {
                    #print "$seq_end_try ";
                    
                    $dbq_map->execute($cosmic_id, $uniprot_id, $seq_end_try);
                    ( $cosmic_pos, $gap ) = $dbq_map->fetchrow_array;
                    $dbq_map->finish;
                    
                    if (!$gap) {    # We're looking after the end of both sequences
                        $gap = "both";
                    }
                    
                    if ($gap ne "none") {
                        $seq_end_try--;
                    }
                } until ($seq_end_try < $seq_start || $gap eq "none");
                if ($gap eq "none") {
                    #print "...using $seq_end_try\n";
                } else {
                    #print "...given up at $seq_end_try\n";
                    $seq_end_try = 0;
                }
            }
        }
    }
    
    my $uniprot_domain_id;
    $dbq_insert_domain->execute($dom_type,  $dom_acc, $dom_id, $uniprot_id, $seq_start, $seq_end, $seq_start_try, $seq_end_try, 0);
    $dbq_domain_id->execute($dom_type,  $dom_acc, $dom_id, $uniprot_id, $seq_start, $seq_end);
    ( $uniprot_domain_id ) = $dbq_domain_id->fetchrow_array;
    $dbq_domain_id->finish;
    
    # Find mutant aggregates within this domain
    
    my ($aggregate_id, $agg_c);
    $agg_c = 0;
    $dbq_aggs->execute($cosmic_id, $seq_start, $seq_end, $seq_start, $seq_end);
    while (($aggregate_id) = $dbq_aggs->fetchrow_array) {
        $dbq_insert_agg->execute($aggregate_id, $uniprot_domain_id);
        $agg_c++;
    }
    
    $dbq_update_mutation->execute($agg_c, $uniprot_domain_id);

    $dbh_titanic->disconnect;
};

