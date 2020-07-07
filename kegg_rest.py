#!/usr/bin/env python

import io
import os
import sys, argparse
import pandas as pd

from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
#from Bio.Graphics.KGML_vis import KGMLCanvas

##### KEGG REST FUNCTIONS
### kegg_conv(target_db, source_db, option=None)
# convert KEGG identifiers to/from outside identifiers.
# Option - Can be "turtle" or "n-triple" (string).
### kegg_find(database, query, option=None)
# Finds entries with matching query keywords or other query data in a given database.
# db - database or organism (string); query - search terms (string); option - search option (string)
### kegg_get(dbentries, option=None)
# Data retrieval
# Identifiers (single string, or list of up to 10 strings)
# Option - One of "aaseq", "ntseq", "mol", "kcf", "image", "kgml" (string)
### kegg_info(database)
# Displays the current statistics of a given database.
### kegg_link(target_db, source_db, option=None)
# find related entries by using database cross-references
# Option - Can be "turtle" or "n-triple" (string).
### kegg_list(database, org=None)
# Entry list for database, or specified database entries.
# db - database or organism (string); org - optional organism (string)

def PDF(filename):
    return HTML('<iframe src=%s width=700 height=350></iframe>' % filename)

def to_df(result):
    return pd.read_table(io.StringIO(result), header=None)

def script_usage():
    print('kegg_rest.py arguments:')
    print('\n Version:')
    print('-v | --version returns statistics for KEGG database')
    print('\n Database Info:')
    print('-i | --info returns info for a given Kegg <database/dbentries>')
    print('usage example: kegg_rest.py --info <pathwayID | speciesID | etc >')
    print('\n Retrieve List:')
    print('-l | --list returns a <database> list for a given -s | --species <speciesID> ')
    print('usage example: kegg_rest.py --list <database; eg. "pathwayID"> --species <speciesID; optional>')
    print('\n Search database:')
    print('At -d | --db <database> search -q | --query <query> with -o | --option <search parameters>')
    print('usage example: kegg_rest.py --db <database; eg. pathway> --query <query; eg. "Repair"> --option <parameters; optional; eg. mol_weight <= n>')
    print('\n Retrieve data:')
    print('-r | --retrieve data for a given keggID <dbentries> and -f | --format <optional; eg. "aaseq", "ntseq", "mol", "kcf", "image", "kgml"> ')
    print('usage example: kegg_rest.py --retrieve <pathwayID>')
    print('usage example: kegg_rest.py --retrieve <pathwayID> --format <"image">')
    print('-g | --genes returns parsed list of genes for a Kegg pathway')
    print('usage example: kegg_rest.py --retrieve <pathwayID> --genes')
    print('\n Convert Identifiers:')
    print('-c | --convert <IDs; database> from/to -t | --targetdb <database>')
    print('For genes <IDs = organism or dbentries for KEGG> AND <ncbi-geneid | ncbi-proteinid | uniprot for outside dbs>')
    print('usage example: kegg_rest.py --convert "crg" --targetdb "uniprot"')
    print('usage example: kegg_rest.py --convert "K1PNA0" --targetdb "uniprot"')
    print('\n Find related entries by dbs cross-references:')
    print('-t | --targetdb <database> -x | --xreference <database | dbentries>')
    print('usage example: kegg_rest.py --targetdb <database; eg. pathway> --xreference <database; eg. enzyme>')
    print('usage example: kegg_rest.py --targetdb <database; eg. pathway> --xreference <dbentries; eg. ec:1.1.1.1>')

def keggInfo(database):
    kegg_info = REST.kegg_info(database).read()
    print(kegg_info)

def keggList(database,species):
    kegg_list = REST.kegg_list(database, species).read()
    print(kegg_list)

def keggSearch(database,query,option):
    kegg_search = REST.kegg_find(database, query, options).read()
    print(kegg_search)

def keggGetData(dbentries,option):
    kegg_getdata = REST.kegg_get(dbentries,option).read()
    print(kegg_getdata)
    return(kegg_getdata)

def keggConvert(targetdb,sourcedb):
    kegg_convert = REST.kegg_conv(targetdb,sourcedb).read()
    print(kegg_convert)

def keggLink(targetdb,sourcedb):
    kegg_link = REST.kegg_link(targetdb,sourcedb).read()
    print(kegg_link)

def parsePathway(pathway):
    print("Genes for -> {0:s}, {1:s}".format(pathway,pathwayName))

def keggPathwayGenes(dbentries):
    kegg_pathway = keggGetData(dbentries)
    gene_list = parsePathway(kegg_pathway)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='store_true')
    parser.add_argument('-i', '--info', help='returns info for a given <database/dbentries>')
    parser.add_argument('-l', '--list', help='obtain a list of entry identifiers and associated definition <database/dbentries>; --species is optional')
    parser.add_argument('-s', '--species', help='narrows a list by organism/species <kegg speciesID>; optional')
    parser.add_argument('-d', '--db', help='find entries in database <db>', choices=['pathway', 'brite', 'module', 'ko', 'genome', 'genes', '<org>', 'vg', 'ag','ligand', 'compound', 'glycan', 'reaction', 'rclass', 'enzyme', 'network', 'variant','disease','drug','dgroup', 'environ','<medicus>'])
    parser.add_argument('-q', '--query', help='find entries with matching query <keyword> or other <query data>')
    parser.add_argument('-o', '--option', help='for compound and drug db, option narrows search by <formula | exact_mass | mol_weight | nop>')
    parser.add_argument('-r', '--retrieve', help='retrieve given database entries <database/dbentries>; --genes and --format are optional')
    parser.add_argument('-g', '--genes', action='store_true', help='parse genes from a given pathway; requires --retrieve')
    parser.add_argument('-f', '--format', choices=["aaseq", "ntseq", "mol", "kcf", "image", "kgml"])
    parser.add_argument('-c', '--convert', help='convert KEGG identifiers to/from outside identifiers <db/dbentries>; kegg db = <org> OR <keggID>; outside db = ncbi-geneid | ncbi-proteinid | uniprot')
    parser.add_argument('-t', '--targetdb', help='target database <database>; use with --convert OR --xreference')
    parser.add_argument('-x', '--xreference', help='find related entries by using database cross-references <sourcedb/dbentries>')
    if len(sys.argv)==1:
        parser.script_usage(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
