package path_setup;
use warnings;
use parent 'Exporter';

use Getopt::Long;

# Package to set up universal paths

our $cosmic_release;
our $xml_chunk = "unset";
our $default_cosmic_release = "v71";
our $mokcatanic_root;
our $username = getpwuid( $< );
if ( $username eq "foop") {
	$mokcatanic_root = "/home/foop/bioinformatics/mokca2/";
	$mokcatanic_local_root = "/local/foop/mokca2/";
} else {
	$mokcatanic_root = "/home/mokca/bioinformatics/m_rebuild/";
	$mokcatanic_local_root = "/local/mokca/m_rebuild/";
}
our $mokcatanic_data_root = $mokcatanic_root . "data/";
our $mokcatanic_web_root = $mokcatanic_root . "web/";
our $tmp_root = "/tmp/";

GetOptions( "release=s" => \$cosmic_release,
            "chunk=s" => \$xml_chunk );

if ( !$cosmic_release ) {
	print STDERR "--- No release specified; using default: $default_cosmic_release\n";
	$cosmic_release = $default_cosmic_release;
}

our @EXPORT = qw($cosmic_release $xml_chunk $default_cosmic_release $mokcatanic_root $username $mokcatanic_data_root $mokcatanic_local_root $mokcatanic_web_root $tmp_root);

1;
