#!/usr/bin/perl -w

use strict;

use Getopt::Long;
use FileHandle;
use DBI();
use Parallel::ForkManager;

# Autoflush STDERR

STDERR->autoflush(1);

# Ticker

my $tick_c = 0;
my $col_wrap = 78;
my $tick_step = 10;
my $tick_char = '.';

# Explain purpose

print STDERR "::: Generate Aggregate Sites - tally aggregate tissue occurences, generate web buttons\n";

# Setup universal paths, database details

use path_setup;
use database_setup;

# Connect to databases

my $dbh_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});

# Prepare queries

my $dbq_agg_idents = $dbh_titanic->prepare(q{SELECT ma.aggregate_id, ch.cosmic_gene_name, ma.mutation_aa, ma.mutation_description FROM mut_aggregate AS ma, cosmic_hgnc AS ch WHERE ma.cosmic_id = ch.cosmic_id});
my $dbq_sites = $dbh_titanic->prepare(q{SELECT primary_site, COUNT(*) FROM cosmic_complete_export WHERE gene_name = ? AND mutation_aa = ? AND mutation_description = ? GROUP BY primary_site});
my $dbq_insert = $dbh_titanic->prepare(q{INSERT INTO agg_sites VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)});

# Go

# All mosaic commands will be written to a text file for reasons lost in the mists of time

my $fm = new Parallel::ForkManager(10); # Parallelise all the things

my $button_fn = $mokcatanic_data_root . "agg_buts.txt";
open(BUTTS, ">$button_fn") || die "+++ Cannot open $button_fn to write data: $?\n";
my ($aggregate_id, $gene_name, $mutation_aa, $mutation_description);
$dbq_agg_idents->execute;
while (( $aggregate_id, $gene_name, $mutation_aa, $mutation_description ) = $dbq_agg_idents->fetchrow_array) {
    
    $fm->start and next;

    agg_sites($aggregate_id, $gene_name, $mutation_aa, $mutation_description);
    
    # Pointless ticker
    
    if ($tick_c % $tick_step == 0) {
        print STDERR $tick_char;
        if ($tick_c % ($tick_step * $col_wrap) == 0 && $tick_c) {
            print STDERR "\n";
        }
    }
    $tick_c++;
    
    $fm->finish;
}
close(BUTTS) || die "+++ Cannot close $button_fn: $?\n";

# Tick finish

print STDERR "\n";

$dbh_titanic->disconnect;

sub agg_sites{
    
    my $aggregate_id = $_[0];
    my $gene_name = $_[1];
    my $mutation_aa = $_[2];
    my $mutation_description = $_[3];
    
    my $png_path = $mokcatanic_web_root . "/tissbuts/agg/" . ($aggregate_id % 100) . "/";
    my $png_fn = $png_path . "agg_$aggregate_id.png";
    
    if (-f $png_fn) {
        return;
    }
    
    my $dbh_sub_titanic = DBI->connect("DBI:mysql:database=$mokca_db;host=$mokca_host",$mokca_user,$mokca_pw, {'RaiseError' => 1});
    my $dbq_insert = $dbh_sub_titanic->prepare(q{INSERT INTO agg_sites VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)});
    my $dbq_sites = $dbh_sub_titanic->prepare(q{SELECT primary_site, COUNT(*) FROM cosmic_complete_export WHERE gene_name = ? AND mutation_aa = ? AND mutation_description = ? GROUP BY primary_site});
    
    my %site_counts;
    
    $site_counts{'autonomic_ganglia'} = 0;
    $site_counts{'central_nervous_system'} = 0;
    $site_counts{'eye'} = 0;
    $site_counts{'meninges'} = 0;
    $site_counts{'pituitary'} = 0;
    $site_counts{'adrenal_gland'} = 0;
    $site_counts{'parathyroid'} = 0;
    $site_counts{'salivary_gland'} = 0;
    $site_counts{'thyroid'} = 0;
    $site_counts{'biliary_tract'} = 0;
    $site_counts{'gastrointestinal_tract_(site_indeterminate)'} = 0;
    $site_counts{'large_intestine'} = 0;
    $site_counts{'oesophagus'} = 0;
    $site_counts{'pancreas'} = 0;
    $site_counts{'small_intestine'} = 0;
    $site_counts{'stomach'} = 0;
    $site_counts{'upper_aerodigestive_tract'} = 0;
    $site_counts{'bone'} = 0;
    $site_counts{'skin'} = 0;
    $site_counts{'soft_tissue'} = 0;
    $site_counts{'breast'} = 0;
    $site_counts{'cervix'} = 0;
    $site_counts{'endometrium'} = 0;
    $site_counts{'fallopian_tube'} = 0;
    $site_counts{'genital_tract'} = 0;
    $site_counts{'female_genital_tract_(site_indeterminate)'} = 0;
    $site_counts{'ovary'} = 0;
    $site_counts{'penis'} = 0;
    $site_counts{'placenta'} = 0;
    $site_counts{'prostate'} = 0;
    $site_counts{'testis'} = 0;
    $site_counts{'paratesticular_tissues'} = 0;
    $site_counts{'vagina'} = 0;
    $site_counts{'vulva'} = 0;
    $site_counts{'haematopoietic_and_lymphoid_tissue'} = 0;
    $site_counts{'thymus'} = 0;
    $site_counts{'kidney'} = 0;
    $site_counts{'liver'} = 0;
    $site_counts{'lung'} = 0;
    $site_counts{'pleura'} = 0;
    $site_counts{'urinary_tract'} = 0;
    $site_counts{'mediastinum'} = 0;
    $site_counts{'midline_organs'} = 0;
    $site_counts{'pericardium'} = 0;
    $site_counts{'peritoneum'} = 0;
    $site_counts{'retroperitoneum'} = 0;
    $site_counts{'NS'} = 0;
    
    $dbq_sites->execute($gene_name, $mutation_aa, $mutation_description);
    while (( my $primary_site, my $agg_c ) = $dbq_sites->fetchrow_array) {
        # Tally all tissue counts for this mutation
        $site_counts{$primary_site} = $agg_c;
    }
    
    # Insert total counts into database
    
    $dbq_insert->execute($aggregate_id, $site_counts{'autonomic_ganglia'} || 0,
    $site_counts{'central_nervous_system'} || 0,
    $site_counts{'eye'} || 0,
    $site_counts{'meninges'} || 0,
    $site_counts{'pituitary'} || 0,
    $site_counts{'adrenal_gland'} || 0,
    $site_counts{'parathyroid'} || 0,
    $site_counts{'salivary_gland'} || 0,
    $site_counts{'thyroid'} || 0,
    $site_counts{'biliary_tract'} || 0,
    $site_counts{'gastrointestinal_tract_(site_indeterminate)'} || 0,
    $site_counts{'large_intestine'} || 0,
    $site_counts{'oesophagus'} || 0,
    $site_counts{'pancreas'} || 0,
    $site_counts{'small_intestine'} || 0,
    $site_counts{'stomach'} || 0,
    $site_counts{'upper_aerodigestive_tract'} || 0,
    $site_counts{'bone'} || 0,
    $site_counts{'skin'} || 0,
    $site_counts{'soft_tissue'} || 0,
    $site_counts{'breast'} || 0,
    $site_counts{'cervix'} || 0,
    $site_counts{'endometrium'} || 0,
    $site_counts{'fallopian_tube'} || 0,
    $site_counts{'genital_tract'} + $site_counts{'female_genital_tract_(site_indeterminate)'},
    $site_counts{'ovary'} || 0,
    $site_counts{'penis'} || 0,
    $site_counts{'placenta'} || 0,
    $site_counts{'prostate'} || 0,
    $site_counts{'testis'} + $site_counts{'paratesticular_tissues'} || 0,
    $site_counts{'vagina'} || 0,
    $site_counts{'vulva'} || 0,
    $site_counts{'haematopoietic_and_lymphoid_tissue'} || 0,
    $site_counts{'thymus'} || 0,
    $site_counts{'kidney'} || 0,
    $site_counts{'liver'} || 0,
    $site_counts{'lung'} || 0,
    $site_counts{'pleura'} || 0,
    $site_counts{'urinary_tract'} || 0,
    $site_counts{'mediastinum'} + $site_counts{'midline_organs'} + $site_counts{'pericardium'} + $site_counts{'peritoneum'} + $site_counts{'retroperitoneum'},
    $site_counts{'NS'} || 0);
    
    # Generate mosaic string
    
    my $mosaic = "montage ";
    if ($site_counts{'autonomic_ganglia'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Agl.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'central_nervous_system'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/CNS.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'eye'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Eye.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'meninges'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Mng.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'pituitary'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Pty.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'adrenal_gland'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Adr.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'parathyroid'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Pth.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'salivary_gland'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Sal.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'thyroid'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Thy.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'biliary_tract'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Bil.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'gastrointestinal_tract_(site_indeterminate)'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/GIT.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'large_intestine'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/LIn.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'oesophagus'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Oes.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'pancreas'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Pnc.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'small_intestine'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/SIn.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'stomach'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Sto.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'upper_aerodigestive_tract'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/UAT.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'bone'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Bon.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'skin'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Ski.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'soft_tissue'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/SoT.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'breast'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Bre.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'cervix'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Cvx.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'endometrium'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/End.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'fallopian_tube'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/FTu.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'genital_tract'} + $site_counts{'female_genital_tract_(site_indeterminate)'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/GTr.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'ovary'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Ova.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'penis'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Pen.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'placenta'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Pla.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'prostate'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Pro.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'testis'} + $site_counts{'paratesticular_tissues'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Tes.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'vagina'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Vag.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'vulva'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Vul.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'haematopoietic_and_lymphoid_tissue'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/HLT.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'thymus'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Thm.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'kidney'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Kid.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'liver'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Lvr.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'lung'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Lng.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'pleura'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/Plr.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'urinary_tract'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/UTr.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'mediastinum'} + $site_counts{'midline_organs'} + $site_counts{'pericardium'} + $site_counts{'peritoneum'} + $site_counts{'retroperitoneum'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/oth.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    if ($site_counts{'NS'} > 0 ) { $mosaic .= " $mokcatanic_web_root/tissbuts/NS.png" } else { $mosaic .= " $mokcatanic_web_root/tissbuts/blank.png" };
    
    $mosaic .= " -tile 14x -geometry +2+2 -background none " . $png_fn;
    print BUTTS "   $mosaic\n";
    
    # Run mosaic string to generate mosaic
    
    system($mosaic);
    
    $dbh_sub_titanic->disconnect;
}
