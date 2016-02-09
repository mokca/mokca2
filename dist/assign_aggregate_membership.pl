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

print STDERR "::: Assign Aggregate Membership - Assign individual mutations to established aggregates\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_aggregates = $dbh_titanic->prepare(q{SELECT aggregate_id, cosmic_gene_name, mutation_aa, mutation_description FROM mut_aggregate AS ma, cosmic_hgnc AS ch WHERE ch.cosmic_id = ma.cosmic_id});
my $dbq_member_c = $dbh_titanic->prepare(q{SELECT COUNT(*) FROM aggregate_membership WHERE aggregate_id = ?});
my $dbq_members = $dbh_titanic->prepare(q{SELECT complete_id FROM cosmic_complete_export WHERE gene_name = ? AND mutation_aa = ? AND mutation_description = ?});
my $dbq_insert = $dbh_titanic->prepare(q{INSERT INTO aggregate_membership VALUES(?,?)});

# Go

my ($aggregate_id, $gene_name, $mutation_aa, $mutation_description);
$dbq_aggregates->execute;
while (($aggregate_id, $gene_name, $mutation_aa, $mutation_description) = $dbq_aggregates->fetchrow_array) {
    
    my $member_c;
    $dbq_member_c->execute($aggregate_id);
    ( $member_c ) = $dbq_member_c->fetchrow_array;
    $dbq_member_c->finish;
    
    if ($member_c == 0) {
        my $complete_id;
        $dbq_members->execute($gene_name, $mutation_aa, $mutation_description);
        while (($complete_id) = $dbq_members->fetchrow_array) {
            $dbq_insert->execute($complete_id, $aggregate_id);
        }
        $dbq_members->finish;
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
