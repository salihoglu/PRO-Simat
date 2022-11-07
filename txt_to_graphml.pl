######################################################################
# Synopsis:
#	This script performs the extraction of interaction information out of a .txt-file 
#	Text file should contain the information in the layout: 
#	'input_node	interaction	output_node'
#	whereas interactions are defined eather as 'activation' or 'inhibition' (make sure to use lower case!)
#	columns have to be separated by '\t'
#
#	The output file will be saved in the active directory with the same filename as the input file just changing the file type:
#	'exampleFile.txt'	-->	'exampleFile.graphml'
#
#	DEPENDENCIES
#	
#	NONE
#
#	COPYRIGHT AND LICENCE
#
#	Copyright (C) 2016 Martin Kaltdorf <martin dot kaltdorf at uni-wuerzburg dot de>
#
#	All Rights Reserved. This module is free software. It may be used,
#	redistributed and/or modified under the same terms as Perl itself.
######################################################################
use strict;
use warnings;
use Data::Dumper;

my %nodes_list_h;
my %edge_list_h;
my @reactions_sort;
my @protlist_sort;

#	definition of input .txt- file. Please change the name according to your input file!
my $file="input.txt";

#	access of input file
open (TXT,$file) or die $!;

#	defining ids for nodes and reactions
my $node_id=0;
my $re_id=0;

#	extract all nodes contained in input file and storage of data in hash '%node_list_h'
while (my $txt=<TXT>)
  {
  chomp $txt;
    #	detection of nodes already contained by '%nodes_list_h' and selection of fitting procedure
    if ($txt=~/^(.+)?\sactivation\s(.+)$/ || $txt=~/^(.+)?\sinhibition\s(.+)$/){
      if ($nodes_list_h{$1} && $nodes_list_h{$2}){
      } 
      elsif($nodes_list_h{$1}){
	$node_id++;
	$nodes_list_h{$2}=$node_id;
      }
      elsif($nodes_list_h{$2}){
	$node_id++;
	$nodes_list_h{$1}=$node_id;
      }
      else{
	$node_id++;
	$nodes_list_h{$2}=$node_id;
	$node_id++;
	$nodes_list_h{$1}=$node_id;
      }
    }
  }
close TXT;

#	creation of list of all reactions/interactions extracted from input file and connected with nodes contained in '%node_list_h'
open (TXT,$file) or die $!."\n File could not be found!\n";
while (my $txt2=<TXT>){
  my $edge;
  my @reac_type_a;
  my $reactand;
  my $product;
  if ($txt2=~/^(.+)?\sactivation\s(.+)/){
    $re_id++;
    $edge=$re_id;
    #	interaction type 'standard' is used for activating edges
    @reac_type_a=("standard","#008000");
    $reactand=$nodes_list_h{$1};
    $product=$nodes_list_h{$2};
  }
  elsif ($txt2=~/^(.+)?\sinhibition\s(.+)/){
    $re_id++;
    $edge=$re_id;
    #	interaction type 't_shape' is used for activating edges
    @reac_type_a=("t_shape","#FF0000");
    $reactand=$nodes_list_h{$1};
    $product=$nodes_list_h{$2};
  }
  #	control of existence of all required information befor saving into '%edge_list_h'
  if (defined $edge && defined $reactand && defined $product){
    my @reaction_complete=($reactand,$product,@reac_type_a);
    @{$edge_list_h{$edge}}=@reaction_complete;
  }
}
close TXT or die $!;

#	alpha-numerical sorting of '%edge_list_h' according to '$edge_id'
@reactions_sort = sort { $a <=> $b } (keys %edge_list_h);
#	alpha-numerical sorting of '%edge_list_h' according to '$node_id'
@protlist_sort = sort { $nodes_list_h{$a} <=> $nodes_list_h{$b} } (keys %nodes_list_h);

#	extraction of filename without file type specifiing ending ('.txt') in order to create the output file accordingly
my $f1;
if ($file =~ /(\w.+)?(\.txt)$/){
  $f1 = $1;
}

#	initialise output file
open (OUT, "> $f1.graphml") or die $!;

#	create output file 
header();
Prot_sub();
Reac_sub();
close_out();
close OUT or die $!;

#	terminal notification
print "network nodes:\t\t".@protlist_sort."\nnetwork edges:\t\t".@reactions_sort."\nnetwork density:\t".(@reactions_sort/@protlist_sort)."\n";
print "\nOutput was saved to \"> $f1.graphml\".\n\n";


###############################
##	subroutines for txt-to-graphml-converter
###############################

#	creation of yEd conform file header
sub header {
  print OUT '<?xml version="1.0" encoding="UTF-8" standalone="no"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:y="http://www.yworks.com/xml/graphml" xmlns:yed="http://www.yworks.com/xml/yed/3" xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns http://www.yworks.com/xml/schema/graphml/1.1/ygraphml.xsd">
  <!--Created by yFiles for Java 2.11-->
  <key for="graphml" id="d0" yfiles.type="resources"/>
  <key for="port" id="d1" yfiles.type="portgraphics"/>
  <key for="port" id="d2" yfiles.type="portgeometry"/>
  <key for="port" id="d3" yfiles.type="portuserdata"/>
  <key attr.name="url" attr.type="string" for="node" id="d4"/>
  <key attr.name="description" attr.type="string" for="node" id="d5"/>
  <key for="node" id="d6" yfiles.type="nodegraphics"/>
  <key attr.name="Beschreibung" attr.type="string" for="graph" id="d7"/>
  <key attr.name="url" attr.type="string" for="edge" id="d8"/>
  <key attr.name="description" attr.type="string" for="edge" id="d9"/>
  <key for="edge" id="d10" yfiles.type="edgegraphics"/>
  <graph edgedefault="directed" id="G">
    <data key="d7"/>';
}

#	creation of yEd conform node list
sub Prot_sub {  
foreach my $prot_key(@protlist_sort){
  print OUT '<node id="';
  print OUT 'n'.$nodes_list_h{$prot_key}.'">';
  print OUT '<data key="d6">
        <y:ShapeNode>
          <y:Geometry height="30.0" width="30.0" x="275.0" y="271.0"/>
          <y:Fill hasColor="false" transparent="false"/>
          <y:BorderStyle color="#000000" type="line" width="1.0"/>
          <y:NodeLabel alignment="center" autoSizePolicy="content" fontFamily="Dialog" fontSize="12" fontStyle="plain" hasBackgroundColor="false" hasLineColor="false" height="17.96875" modelName="custom" textColor="#000000" visible="true" width="26.88671875" x="1.556640625" y="6.015625">'.$prot_key.'<y:LabelModel>
              <y:SmartNodeLabelModel distance="4.0"/>
            </y:LabelModel>
            <y:ModelParameter>
              <y:SmartNodeLabelModelParameter labelRatioX="0.0" labelRatioY="0.0" nodeRatioX="0.0" nodeRatioY="0.0" offsetX="0.0" offsetY="0.0" upX="0.0" upY="-1.0"/>
            </y:ModelParameter>
          </y:NodeLabel>
          <y:Shape type="rectangle"/>
        </y:ShapeNode>
      </data>
    </node>';
  }
}

#	creation of yEd conform reaction list connected to according node-ids
sub Reac_sub {
foreach my $key_r(@reactions_sort)
  {
  my @re_data=@{$edge_list_h{$key_r}};
  if (defined $re_data[2])
    {
#     print "$key_r\n";
    print OUT '<edge id="e'.$key_r.'" source="n'.$re_data[0].'" target="n'.$re_data[1].'">
      <data key="d9"/>
      <data key="d10">
        <y:PolyLineEdge>
          <y:Path sx="0.0" sy="0.0" tx="0.0" ty="0.0">
            <y:Point x="552.0" y="365.0"/>
          </y:Path>
          <y:LineStyle color="'.$re_data[3].'" type="line" width="1.0"/>
          <y:Arrows source="none" target="'.$re_data[2].'"/>
          <y:BendStyle smoothed="false"/>
        </y:PolyLineEdge>
      </data>
    </edge>';
      }
    }
}

#	creation of yEd conform file ending
sub close_out {
print OUT '  </graph>
  <data key="d0">
    <y:Resources/>
  </data>
</graphml>';

}
