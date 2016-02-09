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

print STDERR "::: Aggregate Cosmic Mutations - Create aggregates from COSMIC complete export\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_drop = $dbh_titanic->prepare(q{DROP TABLE IF EXISTS mut_aggregate});
my $dbq_create = $dbh_titanic->prepare(q{CREATE TABLE mut_aggregate (
    aggregate_id INT UNSIGNED NOT NULL AUTO_INCREMENT,
    cosmic_id INT UNSIGNED NOT NULL,
    mutation_aa VARCHAR(256),
    mutation_description ENUM('Complex - compound substitution','Complex - deletion inframe','Complex - frameshift','Complex - insertion inframe','Deletion - Frameshift','Deletion - In frame','Insertion - Frameshift','Insertion - In frame','No detectable mRNA/protein','Nonstop extension','Substitution - Missense','Substitution - Nonsense','Substitution - coding silent','Unknown','Whole gene deletion'),
    aa_start INT UNSIGNED,
    aa_stop INT UNSIGNED,
    map_aa_start INT UNSIGNED,
    map_aa_stop INT UNSIGNED,
    wt_aa VARCHAR(32),
    mut_aa TEXT,
    PRIMARY KEY(aggregate_id),
    KEY(cosmic_id))});
my $dbq_aggs_by_aa = $dbh_titanic->prepare(q{SELECT distinct gene_name, mutation_aa, mutation_description FROM cosmic_complete_export WHERE mutation_aa != ''});
my $dbq_cosmic_id = $dbh_titanic->prepare(q{SELECT cosmic_id FROM cosmic_hgnc WHERE cosmic_gene_name = ?});
my $dbq_uniprot_id = $dbh_titanic->prepare(q{SELECT uniprot_id, identity FROM cosmic_vs_uniprot WHERE cosmic_id = ?});
my $dbq_mapped_res = $dbh_titanic->prepare(q{SELECT hit_pos, query_res, hit_res FROM alignment_mapping WHERE query_db = 'cosmic' AND hit_db = 'uniprot' AND query_id = ? AND hit_id = ? AND query_pos = ? AND gap = 'none'});
my $dbq_insert = $dbh_titanic->prepare(q{INSERT INTO mut_aggregate VALUES (NULL, ?,?,?,?,?,?,?,?,?)});

# Go

# Drop/create table

$dbq_drop->execute;
$dbq_create->execute;

# Fetch distinct mutations (distinctiveness based on gene name, 'mutation aa' and mutation description)

my ($agg_c, $parsable_c, $identical_c, $mappable_c, $mapped_c, $cosmic_fail_c, $uniprot_fail_c);

$agg_c = 0;
$parsable_c = 0;
$identical_c = 0;
$mappable_c = 0;
$mapped_c = 0;
$cosmic_fail_c = 0;
$uniprot_fail_c = 0;

my ($gene_name, $mutation_aa, $mutation_description);
$dbq_aggs_by_aa->execute;
while (( $gene_name, $mutation_aa, $mutation_description ) = $dbq_aggs_by_aa->fetchrow_array) {
    
    $agg_c++;
    
    # Fetch cosmic_id (should work)
    
    my $cosmic_id;
    $dbq_cosmic_id->execute($gene_name);
    ( $cosmic_id ) = $dbq_cosmic_id->fetchrow_array;
    $dbq_cosmic_id->finish;
    if (!$cosmic_id) {
        die "\n+++ Missing ID for COSMIC gene $gene_name\n";
    }
    
    # Fetch mapped uniprot_id (may not work, mapping may not exist if data is shitty)
    
    my ($uniprot_id, $percentage_identity);
    $dbq_uniprot_id->execute($cosmic_id);
    ( $uniprot_id, $percentage_identity ) = $dbq_uniprot_id->fetchrow_array;
    $dbq_uniprot_id->finish;
    
    # Initialise mappings as zero/empty in case mapping fails at any stage
    
    my ($wt_aa, $mut_pos, $mut_aa);
    $wt_aa = 'NONE';
    $mut_pos = 0;
    $mut_aa = '';
    my ($aa_start, $aa_stop, $map_aa_start, $map_aa_stop);
    $aa_start = 0;
    $aa_stop = 0;
    $map_aa_start = 0;
    $map_aa_stop = 0;
    
    
    if ($mutation_aa =~ /p.([\*A-Z])(\d+)([\*A-Z])/) { # p.A204E p.N49* p.R36R
        
        # Simple missense, silent or nonsense mutation
        
        $wt_aa = $1;
        $mut_pos = $2;
        $mut_aa = $3;
        
        $aa_start = $mut_pos;
        $aa_stop = $mut_pos;
        
        $parsable_c++;
    } elsif ($mutation_aa =~ /p.[\*A-Z](\d+)fs\*>*(\d+)/) { # p.E82fs*12 p.D204fs*>24 p.*710fs*>2
        
        # Frame shift: count as going from start position to new stop codon.
        
        $aa_start = $1;
        $aa_stop = $aa_start + $2 - 1;
        
        $parsable_c++;
    } elsif ($mutation_aa =~ /p.[A-Z](\d+)del[A-Z]*/) { # p.W557delW p.H580del
        $aa_start = $1;
        $aa_stop = $aa_start;
        
        $parsable_c++;
    } elsif ($mutation_aa =~ /p.(\d+)_(\d+)del\d+/) { # p.982_1028del47
        $aa_start = $1;
        $aa_stop = $2;
        
        $parsable_c++;
    } elsif ($mutation_aa =~ /p.[A-Z](\d+)_[\*A-Z](\d+)del[\*A-Z]*/) { # p.P585_R587delPNR p.L262_*730del p.S663_*665delSY*
        $aa_start = $1;
        $aa_stop = $2;
        
        $parsable_c++;
    } elsif ($mutation_aa =~ /p.[A-Z]+(\d+)_(\d+)del/) { # p.TLV1573_1575del
        $aa_start = $1;
        $aa_stop = $2;
        
        $parsable_c++;
    } elsif ($mutation_aa =~ /p[A-Z](\d+)_[A-Z](\d+)del/) { # pQ48_D54del
        $aa_start = $1;
        $aa_stop = $2;
        
        $parsable_c++;
    } elsif ($mutation_aa =~ /p.[\*A-Z]*(\d+)_[\*A-Z]*(\d+)>[\*A-Z]+/) { # p.E746_S752>T p.W236_E237>* p.747_752>S p.T473_*477>R
        $aa_start = $1;
        $aa_stop = $2;
        
        $parsable_c++;
    } elsif ($mutation_aa =~ /p.[\*A-Z](\d+)>[\*A-Z]+/) { # p.A2141>HN p.Y62>*
        $aa_start = $1;
        $aa_stop = $1+1;
        
        $parsable_c++;
    } elsif ($mutation_aa =~ /p.[\*A-Z]*(\d+)_[\*A-Z]*(\d+)ins[A-Z0-9]*/) { # p.L242_L243insRL p.E596_Y597ins12  p.794_795insD p.N1237_*1238insN
        $aa_start = $1;
        $aa_stop = $2;
        
        $parsable_c++;
    } elsif ($mutation_aa eq "p.?"
    || $mutation_aa eq "p.?fs"
    || $mutation_aa eq "p.0?"
    || $mutation_aa eq "p.(=)"
    || $mutation_aa eq "p.?_?ins?"
    || $mutation_aa eq "p.?fs*?"
    || $mutation_aa eq "p.?del"
    || $mutation_aa eq "p.0"
    || $mutation_aa eq "p.>"
    || $mutation_aa eq "p.fs"
    || $mutation_aa eq "p.unknown"
    || $mutation_aa eq "p.fs*?"
    || $mutation_aa eq "?"
    || $mutation_aa eq "p."
    || $mutation_aa eq "p.?*"
    || $mutation_aa eq "p.??>*"
    || $mutation_aa eq "p.?_?del"
    || $mutation_aa eq "p.?ins?"
    || $mutation_aa =~ /p.[A-Z]\d+>\?/          # p.T210>?
    || $mutation_aa =~ /p.[\*A-Z]\d+\?/         # p.Q61? p.*553?
    || $mutation_aa =~ /p.[\*A-Z]\d+fs/         # p.M802fs p.*704fs?
    || $mutation_aa =~ /p.[A-Z]\d+/             # p.R2468
    || $mutation_aa =~ /p.\d+_\d+>\d+/          # p.612_613>17
    || $mutation_aa =~ /p.\([A-Z]\d+\)fs\**\?*/ # p.(N216)fs p.(Y346)fs*?
    || $mutation_aa =~ /p.\(\d+_\d+\)?/         # p.(552_596)?
    || $mutation_aa =~ /p.\(\d+_\d+\)fs\*\?/    # p.(574_1542)fs*?
    || $mutation_aa =~ /p.\d+fs\**\?*\d*/       # p.419fs*? p.553fs*7 p.215fs
    || $mutation_aa =~ /p.\(\d+\)fs\**\?*/      # p.(2287)fs p.(386)fs*?
    || $mutation_aa =~ /p.\([A-Z]\d+\)ins?/     # p.(V769)ins?
    || $mutation_aa =~ /p.\(\d+\)ins\d+/        # p.(1409)ins6
    || $mutation_aa =~ /p. [A-Z]\d+fs\*\d+/     # p. F33fs*107
    || $mutation_aa =~ /p.[\*A-Z]\d+del\*/      # p.*707del*
    || $mutation_aa =~ /p.\d+_\d+del/           # p.1670_1673del
    || $mutation_aa =~ /p.[A-Z]+\d+\?/          # p.GIHS34?
    || $mutation_aa =~ /p.\d+_\d+>/             # p.5_142>
    || $mutation_aa =~ /p.\d+>/                 # p.30>
    || $mutation_aa =~ /p.\?_\?ins\d+/          # p.?_?ins6
    || $mutation_aa =~ /p.\?_\?ins[A-Z]+/       # p.?_?insXXXX
    || $mutation_aa =~ /p.\d+_\d+dup/           # p.556_871dup
    || $mutation_aa =~ /p.\?ins\d+/             # p.?ins22
    || $mutation_aa =~ /p.\?[A-Z]>[A-Z]/        # p.?A>V
    || $mutation_aa =~ /p.\?fs\*\(\d+_\d+\)/    # p.?fs*(46_47)
    || $mutation_aa =~ /p.\([A-Z]\d+_[A-Z]\d+\)ins[A-Z]+/           # p.(I682_E684)insQG
    || $mutation_aa =~ /p.\([A-Z]\d+_[A-Z]\d+\)[A-Z]/               # p.(K1302_K1303)L
    || $mutation_aa =~ /p.\([A-Z]\d+_[A-Z]\d+\)fs\*\(\d+[\-_]\d+\)/ # p.(H506_L508)fs*(143-145)
                                                                    # p.(H667_H668)fs*(13_14)
    || $mutation_aa =~ /p.\([A-Z]\d+_[A-Z]\d+\)fs\**\?*/            # p.(S354_G355)fs*? p.(G67_A68)fs
    || $mutation_aa =~ /p.\([A-Z]\d+_[A-Z]\d+\)\?/                  # p.(S505_W515)?
    ) {
        
        # Discard this: unknown, uninterpretable, or woefully broken mutation string.
        
    } else {
        print "::: $mutation_aa\n";
    }
    
    # If we have a uniprot mapping and a mutation position (which will be both start and stop until
    # the code is extended to map longer mutations) attempt to map the cosmic mutation position onto
    # the uniprot mutation position
    
    if ($uniprot_id && $aa_start) {
        if ($percentage_identity == 100.0) {
            
            $identical_c++;
            
            # If we know the mapping is at 100% identity, we don't need to be clever
            
            $map_aa_start = $aa_start;
            $map_aa_stop = $aa_stop;
            
            # Let's hope the wt_aa is correct, eh?
            
        } else {
            
            $mappable_c++;
            
            # We need to use the stored COSMIC/uniprot alignment to map residue numbers
            
            my ($map_pos, $cosmic_res, $uniprot_res);
            
            $dbq_mapped_res->execute($cosmic_id, $uniprot_id, $aa_start);
            ( $map_pos, $cosmic_res, $uniprot_res ) = $dbq_mapped_res->fetchrow_array;
            $dbq_mapped_res->finish;
            
            if ($map_pos) {
                
                $mapped_c++;
                
                # If we actually have a mapping...
                
                if ($wt_aa ne "NONE") {
                    
                    # Residue sanity check
                    
                    if ($cosmic_res ne $wt_aa) {
                        print STDERR "\n+++ $gene_name $mutation_aa: COSMIC wt actually $cosmic_res\n";
                        $cosmic_fail_c++;
                        $tick_c = 0;
                    } elsif ($uniprot_res ne $wt_aa) {
                        print STDERR "\n+++ $gene_name $mutation_aa: UNIPROT wt actually $uniprot_res\n";
                        $uniprot_fail_c++;
                        $tick_c = 0;
                    }
                }
                
                # Should a failed sanity check stop us making the mapping?
                
                $map_aa_start = $map_pos;
            }
            
            if ($aa_stop) {
                $dbq_mapped_res->execute($cosmic_id, $uniprot_id, $aa_stop);
                ( $map_pos, $cosmic_res, $uniprot_res ) = $dbq_mapped_res->fetchrow_array;
                $dbq_mapped_res->finish;
                
                if ($map_pos) {
                    $map_aa_stop = $map_pos;
                }
            }
        }
    }
    
    $dbq_insert->execute($cosmic_id, $mutation_aa, $mutation_description, $aa_start, $aa_stop, $map_aa_start, $map_aa_stop, $wt_aa, $mut_aa);
    
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
print STDERR "--- $agg_c aggregates\n    $parsable_c parsed\n    $identical_c don't need mapping, $mappable_c do\n    $mapped_c were mapped\n    $cosmic_fail_c COSMIC residue mismatches, $uniprot_fail_c UNIPROT residue mismatches\n";


$dbh_titanic->disconnect;
