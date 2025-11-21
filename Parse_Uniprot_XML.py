from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import xml.etree.ElementTree as ET
import pandas as pd
import argparse
import time

#Xample xml entry from uniprot:
# <?xml version="1.0" encoding="UTF-8"?>
# <uniprot xmlns="https://uniprot.org/uniprot"
#  xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
#  xsi:schemaLocation="http://uniprot.org/uniprot http://www.uniprot.org/docs/uniprot.xsd">
# <entry dataset="Swiss-Prot" created="2011-06-28" modified="2025-04-09" version="36" xmlns="https://uniprot.org/uniprot">
#     <accession>Q6GZQ5</accession>
#     <name>070R_FRG3G</name>
#     <protein>
#         <recommendedName>
#             <fullName>Uncharacterized protein 070R</fullName>
#         </recommendedName>
#     </protein>
#     <gene>
#         <name type="ORF">FV3-070R</name>
#     </gene>
#     <organism>
#         <name type="scientific">Frog virus 3 (isolate Goorha)</name>
#         <name type="common">FV-3</name>
#         <dbReference type="NCBI Taxonomy" id="654924"/>
#         <lineage>
#             <taxon>Viruses</taxon>
#             <taxon>Varidnaviria</taxon>
#             <taxon>Bamfordvirae</taxon>
#             <taxon>Nucleocytoviricota</taxon>
#             <taxon>Megaviricetes</taxon>
#             <taxon>Pimascovirales</taxon>
#             <taxon>Iridoviridae</taxon>
#             <taxon>Alphairidovirinae</taxon>
#             <taxon>Ranavirus</taxon>
#             <taxon>Frog virus 3</taxon>
#         </lineage>
#     </organism>
#     <organismHost>
#         <name type="scientific">Dryophytes versicolor</name>
#         <name type="common">chameleon treefrog</name>
#         <dbReference type="NCBI Taxonomy" id="30343"/>
#     </organismHost>
#     <organismHost>
#         <name type="scientific">Lithobates pipiens</name>
#         <name type="common">Northern leopard frog</name>
#         <name type="synonym">Rana pipiens</name>
#         <dbReference type="NCBI Taxonomy" id="8404"/>
#     </organismHost>
#     <organismHost>
#         <name type="scientific">Lithobates sylvaticus</name>
#         <name type="common">Wood frog</name>
#         <name type="synonym">Rana sylvatica</name>
#         <dbReference type="NCBI Taxonomy" id="45438"/>
#     </organismHost>
#     <organismHost>
#         <name type="scientific">Notophthalmus viridescens</name>
#         <name type="common">Eastern newt</name>
#         <name type="synonym">Triturus viridescens</name>
#         <dbReference type="NCBI Taxonomy" id="8316"/>
#     </organismHost>
#     <reference key="1">
#         <citation type="journal article" date="2004" name="Virology" volume="323" first="70" last="84">
#             <title>Comparative genomic analyses of frog virus 3, type species of the genus Ranavirus (family Iridoviridae).</title>
#             <authorList>
#                 <person name="Tan W.G."/>
#                 <person name="Barkman T.J."/>
#                 <person name="Gregory Chinchar V."/>
#                 <person name="Essani K."/>
#             </authorList>
#             <dbReference type="PubMed" id="15165820"/>
#             <dbReference type="DOI" id="10.1016/j.virol.2004.02.019"/>
#         </citation>
#         <scope>NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA]</scope>
#     </reference>
#     <dbReference type="EMBL" id="AY548484">
#         <property type="protein sequence ID" value="AAT09730.1"/>
#         <property type="molecule type" value="Genomic_DNA"/>
#     </dbReference>
#     <dbReference type="RefSeq" id="YP_031649.1">
#         <property type="nucleotide sequence ID" value="NC_005946.1"/>
#     </dbReference>
#     <dbReference type="KEGG" id="vg:2947770"/>
#     <dbReference type="Proteomes" id="UP000008770">
#         <property type="component" value="Segment"/>
#     </dbReference>
#     <proteinExistence type="predicted"/>
#     <keyword id="KW-1185">Reference proteome</keyword>
#     <feature type="chain" id="PRO_0000410523" description="Uncharacterized protein 070R">
#         <location>
#             <begin position="1"/>
#             <end position="124"/>
#         </location>
#     </feature>
#     <feature type="region of interest" description="Disordered" evidence="1">
#         <location>
#             <begin position="93"/>
#             <end position="124"/>
#         </location>
#     </feature>
#     <feature type="compositionally biased region" description="Polar residues" evidence="1">
#         <location>
#             <begin position="102"/>
#             <end position="111"/>
#         </location>
#     </feature>
#     <evidence type="ECO:0000256" key="1">
#         <source>
#             <dbReference type="SAM" id="MobiDB-lite"/>
#         </source>
#     </evidence>
#     <sequence length="124" mass="13371" checksum="5DD06098A99407B7" modified="2004-07-19" version="1">MASHYYSKRPERPSDGELASIVAEAAARVLSKYGLKVRDPPAFSAAASASLSRADSDPSTIPMGMNRRQTAVYFTMKGMLADASARAVVVQPRSVHPAHPSTHFNGTSSAVRPSRHYNAPGRFR</sequence>
# </entry>

parser = argparse.ArgumentParser()
parser.add_argument("--xml_files", help="Xml files to parse (accepts multiple)", type =str, nargs="+")
parser.add_argument("--info_output", help="The output table file with the generated data for the sequences", default="Seq_Info.tsv", type =str)
parser.add_argument("--fasta_output", help="The output fasta file", default="Sequences.fasta", type =str)
args = parser.parse_args()

    
def parseXML(xml_files):
     # create empty list for entries
    seq_info = defaultdict(dict)
    OUTSEQ = open(args.fasta_output,'w', newline='')
    OUTINFO = open(args.info_output,'w', newline='')
    #print the header to outinfo
    OUTINFO.write("Accession\tName\tDescription\tLineage\tHost\tPfam_Domains\n")
    for xmlfile in xml_files:
        print(f"Parsing {xmlfile}")
        root = None
        for event, elem in ET.iterparse(xmlfile, events=("start", "end")):
            if event =='start':
                if root is None: root = elem
            elif event == 'end':
                if len(root) > 0 and elem == root[0]:
                    del root[0]
                # process the tag
                # if elem.tag == '{https://uniprot.org/uniprot}entry' and len(path) == 2:
                if elem.tag == '{https://uniprot.org/uniprot}entry':
                    # print(f"Found entry: {elem.tag} with attributes {elem.attrib}")
                    accession = elem.find('./uniprot:accession', namespaces={'uniprot': 'https://uniprot.org/uniprot'})
                    name = elem.find('./uniprot:name', namespaces={'uniprot': 'https://uniprot.org/uniprot'})
                    #Get the amino acid sequence from the item
                    sequence = elem.find('./uniprot:sequence', namespaces={'uniprot': 'https://uniprot.org/uniprot'})
                    #get protein, <recommendedName>. fullName
                    description = elem.find('./uniprot:protein/uniprot:recommendedName/uniprot:fullName', namespaces={'uniprot': 'https://uniprot.org/uniprot'})
                    #get all the taxa in lineage field
                    taxa_list = []
                    lineage = elem.findall('./uniprot:organism/uniprot:lineage/uniprot:taxon', namespaces={'uniprot': 'https://uniprot.org/uniprot'})
                    for element in lineage:
                        taxon = element.text
                        taxa_list.append(taxon)
                    #get the scientific name of the organism (last level of the taxonomy)
                    tax_name = elem.find('./uniprot:organism/uniprot:name[@type="scientific"]', namespaces={'uniprot': 'https://uniprot.org/uniprot'})
                    if tax_name is not None:
                        taxa_list.append(tax_name.text)
                    #get the tax scientific names and tax ids of the hosts
                    host_names = []
                    host_tids =[]
                    host_elems = elem.findall('./uniprot:organismHost', namespaces={'uniprot': 'https://uniprot.org/uniprot'})
                    if host_elems:
                        for helement in host_elems:
                            host_name = helement.find('./uniprot:name[@type="scientific"]', namespaces={'uniprot': 'https://uniprot.org/uniprot'})
                            #host_name = helement
                            if host_name is not None:
                                host_names.append(host_name.text)
                    #get Pfam info
                    pfam_doms = set()
                    pfam_elems = elem.findall('./uniprot:dbReference[@type="Pfam"]', namespaces={'uniprot': 'https://uniprot.org/uniprot'})
                    if pfam_elems:
                        for pelement in pfam_elems:
                            pfam_id = pelement.attrib.get('id')
                            if pfam_id is not None:
                                pfam_doms.add(pfam_id)
                    #paste the elements of taxa_list into a tring separated by ;
                    if (sequence.text is not None and sequence.text.strip() != "") and (accession is not None and accession.text is not None):
                        # print(f"adding info for {accession.text}")
                        final_description = description.text if description is not None else "Uncharacterized protein"
                        # seq_info["Description"][accession.text] = final_description
                        # seq_info["Name"][accession.text] = name.text
                        # if name is not None:
                        #     seq_info["Name"][accession.text] = name.text
                        # # if description is not None:
                        # # seq_info["Description"][accession.text] = description.text
                        # if taxa_list:
                        #     seq_info["Lineage"][accession.text] = ';'.join(taxa_list)
                        # if host_names:
                        #     seq_info["Host"][accession.text] = ';'.join(host_names)
                        #seq_info["Host_Tax_IDs"][accession.text] = ';'.join(host_tids)
                        tax_string = ';'.join(taxa_list)
                        host_string = ';'.join(host_names) 
                        full_desc_string = "|".join([final_description, tax_string]) #, host_string
                        pfam_string = ';'.join(pfam_doms) if pfam_doms else 'NA'
                        # add pfam info to description
                        seqobj = SeqRecord(id= accession.text,seq=Seq(sequence.text), description=full_desc_string) #description=final_description
                        SeqIO.write(seqobj, OUTSEQ, "fasta")
                        # write the info to the tsv file
                        OUTINFO.write(f"{accession.text}\t{name.text if name is not None else 'NA'}\t{final_description}\t{tax_string}\t{host_string}\t{pfam_string}\n")

    return seq_info



def savetoTSV(seq_info=None, filename=args.info_output):
    info_df = pd.DataFrame.from_dict(seq_info)
    info_df.index.name = 'Sequence'
    info_df.to_csv(filename,sep="\t",na_rep='NA')

    
def main():
    # parse xml file
    seq_info = parseXML(args.xml_files)
    # store item_dict items in a tsv file
    #savetoTSV(seq_info=seq_info)
    
    
if __name__ == "__main__":
    # calling main function
    main()

