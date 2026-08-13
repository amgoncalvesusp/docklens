#!/usr/bin/perl -w
use strict;
use warnings;
use MdmDiscoveryScript;

my $document = Mdm::Document::Create({ NoDisplay => True })
    or die "Could not create an in-memory Discovery Studio document\n";
my $molecule = $document->CreateMolecule()
    or die "Could not create an in-memory molecule\n";
$molecule->CreateAtom("O")
    or die "Could not create an in-memory atom\n";
my $monitor = $document->CreateNonbondMonitor($document->Atoms);

print "# Discovery Studio non-bond defaults\n";
print "# criterion\tcurrent\tdefault\n";
my $criteria = Mdm::NonbondMonitor::AllNonbondCriterionTypes();
for my $criterion (@$criteria) {
    printf "%s\t%.6f\t%.6f\n",
        $criterion,
        $monitor->NonbondCriterion($criterion),
        $monitor->DefaultNonbondCriterion($criterion);
}

print "# nonbond-types\n";
my $types = Mdm::NonbondMonitor::AllNonbondTypes();
print join("\n", @$types), "\n";
