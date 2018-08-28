#!/usr/bin/python3
# -*- coding: UTF-8 -*-
#
# INTRODUCTION
# This script is designed to run CobraPY analysis automatically.
#
# DEPENDENCIES
# =============================================================================
# Python (>=2.7, 3.x is recommended for its speed.)
# o CobraPY (>=0.13)
# =============================================================================
#
# USAGE
# =============================================================================
# calculate_ani.py [options]
#
# Options:
#   -h, --help            show this help message and exit
#   -m, --model           input model
#   -r, --reactions       target reaction(s)
#   -v, --version         program version
#   -p, --ppi             species protein-protein interaction (PPI) file
#   -a, --annotation      gene annotation file (from BioMart)
#   -t, --threads         How many threads will be used? [default all]
# ==============================================================================
#
# The MIT License
#
# Copyright (c) 2018 Northwest A&F University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ===============================================================================

import os
import re
import shutil
import cobra
import requests
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib as mpl
import warnings
import svgutils.compose as sc
from bioservices.kegg import KEGG
from urllib.request import urlretrieve
from pybiomart import Dataset, Server
from argparse import ArgumentParser
from datetime import datetime
from collections import defaultdict, OrderedDict
from multiprocessing import Pool, cpu_count


def parse_cmdline():
    """
    Parse command-line arguments for script.
    :return: Input command-line arguments
    """
    parser = ArgumentParser(prog="cobra_analysis.py")
    parser.add_argument("-m", "--model", dest="model", action="store", default=None,
                        required=True, help="Input model file")
    parser.add_argument("-r", "--reactions", dest="reactions", action="append",
                        default=None, required=True, help="Input reaction")
    parser.add_argument("-d", "--dataset", dest="dataset", action="store",
                        default=None, required=True, help="BioMart organism dataset")
    parser.add_argument("-p", "--ppi", dest="ppi", action="store", default=None,
                        required=True, help="PPI file")
    parser.add_argument("-c", "--metacyc", dest="metacyc", action="store", default=None,
                        required=False, help="Metacyc pathway file")
    parser.add_argument("-s", "--species", dest="species", action="store", default=None,
                        required=True, help="KEGG species abbreviation")
    parser.add_argument("-v", "--version", action="version", version="%(prog)s v1.0",
                        default=None, help="Show program’s version number and exit")
    parser.add_argument("-n", "--cpu", dest="cpu_num", action="store", required=False,
                        default=cpu_count(), help="How many processes will be used")
    parser.add_argument("-t", "--taxon", dest="taxon", action="store", required=True,
                        default=None, help="Input STRING organisms taxon ID")
    parser.add_argument("-o", "--output", dest="output", action="store", required=False,
                        default=None, help="Output dir name")
    parser.add_argument("-e", "--essential", dest="essential", action="store", required=False,
                        default=None, help="DEG essential gene file")
    parsers = parser.parse_args()
    return parsers


def get_biomart(species, meta):
    tmp_host = 'http://asia.ensembl.org'
    server = Server(host=tmp_host)
    query_set = None
    try:
        dataset = Dataset(name=species, host=tmp_host)
        if meta:
            query_set = dataset.query(attributes=['ensembl_gene_id',
                                                  'external_gene_name',
                                                  'description',
                                                  'uniprotswissprot',
                                                  'kegg_enzyme',
                                                  'metacyc'])
        else:
            query_set = dataset.query(attributes=['ensembl_gene_id',
                                                  'external_gene_name',
                                                  'description',
                                                  'uniprotswissprot',
                                                  'kegg_enzyme'])
    except IndexError:
        mart = server['ENSEMBL_MART_ENSEMBL']
        print('Invalid dataset in BioMart')
        print(mart.list_datasets())
    return query_set


def get_flux(index, model):
    try:
        flux_value = float(model.reactions[index].x)
    except IndexError:
        flux_value = None
    return round(flux_value, 4)


def essential_gene_knockout(model_file):
    re_pattern = re.compile(r"frozenset\({'(.*)'\}\)")
    model = cobra.io.read_sbml_model(model_file)
    model.optimize()
    deletion_results = cobra.flux_analysis.deletion.single_gene_deletion(model)
    tmp_result_file = os.path.join(my_path, 'tmp.csv')
    deletion_results.to_csv(tmp_result_file)
    knockout_dict = OrderedDict()
    with open(tmp_result_file, 'r') as f:
        for each_line in f.readlines()[1:]:
            a_list = each_line.strip().split(',')
            m = re.search(pattern=re_pattern, string=a_list[0])
            gene_name = m.groups()[0]
            flux_value = round(float(a_list[1]), 4)
            knockout_dict[gene_name] = str(flux_value)
    os.remove(tmp_result_file)
    return knockout_dict


def load_gene_annotation(query_set=pd.DataFrame(), meta=None):
    gene_dict = defaultdict()
    uniprot_dict = defaultdict()
    meta_dict = defaultdict()
    for index, row in query_set.iterrows():
        table = row['Gene stable ID']
        name = row['Gene name']
        description = row['Gene description']
        uniprot_id = row['UniProtKB/Swiss-Prot ID']
        # kegg_list = row['KEGG Pathway and Enzyme ID'].split('+')
        # kegg_enzymes = uniprot_id
        uniprot_dict[table] = uniprot_id
        if meta:
            metacyc_id = row['MetaCyc ID']
            meta_dict.setdefault(table, []).append(metacyc_id)
        gene_dict[table] = [name, description]
    return gene_dict, uniprot_dict, meta_dict


def map_values_to_colors(values, non_negative_colormap='inferno',
                         divergent_color_map='RdBu_r'):
    """ if values are non-negative, a sequencial colormap is used,
    if the values are negative and positive, a symetrical Red to blue colormap is used"""
    if all(np.array(values) >= 0):
        norm = None
        cmap = non_negative_colormap
    else:
        extreme = np.abs(values).max()
        norm = mpl.colors.Normalize(vmin=-extreme, vmax=extreme)
        cmap = divergent_color_map
    color_scaler = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    rgb = color_scaler.to_rgba(values)[:, :3]
    color_codes = ['RGB({:.0f},{:.0f},{:.0f})'.format(*tuple(row)) for row in rgb * 255]
    return color_codes


def map_values_to_range(values, output_range, vmin, vmax):
    if output_range is None:
        output_range = [0, 1]
    np_values = np.array(values)
    if vmin is None:
        vmin = np_values.min()
    if vmax is None:
        vmax = np_values.max()
    if vmin < 0:
        warnings.warn('minimum of values is negative, '
                      'data will be mapped from [{},{}] -> [{},{}]'.format(vmin, vmax,
                                                                           *tuple(output_range)))
    diff = vmax - vmin
    # clip values
    values[np.where(values > vmax)[0]] = vmax
    values[np.where(values < vmin)[0]] = vmin
    normed_values = (values - vmin) / diff
    assert len(output_range) == 2
    out_diff = output_range[1] - output_range[0]
    assert out_diff > 0
    return normed_values * out_diff + output_range[0]


def create_selection(data, color_column=None, width_column=None, opacity_column=None,
                     color_kws=None, width_kws=None, opacity_kws=None):
    """transforms a pandas dataframe to selection which can be used in ipath2:
    # Supported data types for indexes
    The following data types are supported by iPath, and can be used to customize the maps.
    Make sure you use the required prefix for each data type, as shown in the examples below:
    Data type	Prefix	Example ID
    KEGG Pathways	-	00650
    KEGG Compounds	-	C00003
    KEGG KOs	-	K01000
    STRING proteins	-	224324.AQ_626
    KEGG proteins	-	aae:aq_626
    COGs/eggNOG OGs	-	COG0007
    Enzyme EC numbers	E or EC	E2.4.1.82
    Uniprot IDs/ACCs	UNIPROT:	UNIPROT:Q93015
    IPI IDs	-	IPI00745889
    NCBI GI IDs	ncbi-gi:	ncbi-gi:326314893
    """
    assert type(data) == pd.DataFrame
    output_data = pd.DataFrame(index=data.index)
    if color_column:
        if color_kws:
            color_kws = defaultdict()
        output_data['color'] = map_values_to_colors(data[color_column], **color_kws)
    if width_column:
        width_default_param = dict(output_range=[0, 50])
        if width_kws:
            width_default_param.update(width_kws)
        output_data['Width'] = map_values_to_range(data[width_column],
                                                   **width_default_param).astype(str)
        output_data['Width'] = "W" + output_data['Width']
    if opacity_column:
        if opacity_kws is None:
            opacity_kws = defaultdict()
        output_data['Opacity'] = map_values_to_range(data[opacity_column],
                                                     **opacity_kws).astype(str)
    output_str = ''
    for i, row in output_data.iterrows():
        output_str += ' '.join([str(i)] + list(row)) + '\n'
    return output_str


def get_ppi_image(taxon_id, my_gene, reaction_dir):
    string_api_url = 'https://string-db.org/api'
    output_format = 'image'
    method = 'network'
    my_app = 'www.awesome_app.org'
    # Construct the request
    request_url = string_api_url + '/' + output_format + '/' + method + '?'
    request_url += 'identifiers=%s'
    request_url += '&' + 'species=' + taxon_id
    request_url += '&' + 'add_white_nodes=15'
    request_url += '&' + 'caller_identity=' + my_app
    result_file = os.path.join(reaction_dir, '{0}_ppi.png'.format(my_gene))
    # For each gene call STRING
    urlretrieve(request_url % "%0d".join(list(my_gene)), result_file)


def to_parameters(selection, export_type='svg', include_metabolic=True, include_secondary=False,
                  include_antibiotic=False, include_microbial=False, whole_modules=False,
                  whole_pathways=False, keep_colors=False, default_opacity=1, default_width=3,
                  default_radius=7, default_color='#666666', query_reactions=False, tax_filter='',
                  export_dpi=1200):
    allowed_export_types = ['svg', 'png', 'pdf', 'eps']
    assert export_type in allowed_export_types, \
        "export_type {0} needs to be one of {1}".format(export_type,
                                                        allowed_export_types)
    result_dict = dict(selection=selection,
                       export_type=export_type,
                       keep_colors=int(keep_colors),
                       include_metabolic=int(include_metabolic),
                       include_secondary=int(include_secondary),
                       include_antibiotic=int(include_antibiotic),
                       include_microbial=int(include_microbial),
                       whole_modules=int(whole_modules),
                       default_opacity=default_opacity,
                       whole_pathways=int(whole_pathways),
                       default_width=default_width,
                       default_color=default_color,
                       default_radius=default_radius,
                       query_reactions=int(query_reactions),
                       tax_filter=tax_filter,
                       export_dpi=export_dpi)
    return result_dict


def get_ipath_map(selection, reaction_dir, map_name='map', **param):
    url = 'https://pathways.embl.de/mapping.cgi'
    parameters = to_parameters(selection, **param)
    r = requests.post(url, data=parameters)
    assert r.ok, r.text
    result_file = os.path.join(reaction_dir, '{0}.svg'.format(map_name))
    with open(result_file, 'w') as file:
        file.write(r.text)


def inspect_online(selection, **param):
    url = 'https://pathways.embl.de/ipath.cgi'
    parameters = to_parameters(selection, **param)
    r = requests.post(url, data=parameters)
    return r


def load_kegg(gene, organism):
    k = KEGG()
    result_line = ''
    try:
        a = k.get_pathway_by_gene(gene, organism)
        if a:
            k_list = list(a.values())
            result_line = ', '.join(k_list)
    except:
        pass
    return result_line


def scale_map(map_name):
    sc.Figure("26cm", "16cm",
              sc.Panel(sc.SVG(map_name + ".svg").scale(0.265)),
              ).save(map_name + '_scaled.svg')


def run_cobra(reaction_index, model_file, index):
    # normal
    model = cobra.io.read_sbml_model(model_file)
    model.optimize()
    normal_flux = get_flux(reaction_index, model)
    # knockout
    model = cobra.io.read_sbml_model(model_file)
    cobra.manipulation.delete.remove_genes(model, [model.genes[index]])
    model.optimize()
    knockout_flux = get_flux(reaction_index, model)
    # stat
    result = []
    if normal_flux:
        if normal_flux != knockout_flux:
            diff = round(normal_flux - knockout_flux, 4)
            result = [str(model.genes[index].id), str(normal_flux),
                      str(knockout_flux), str(diff)]
    return result


def load_ppi(ppi_file):
    ppi_dict = defaultdict()
    edge_list = []
    ppi_graph = nx.Graph()
    with open(ppi_file, 'r') as f:
        for each_line in f.readlines():
            p_list = each_line.strip('\r|\n').split(' ')
            gene_1 = p_list[1]
            gene_2 = p_list[3]
            ppi_value = int(p_list[4])
            ppi_dict.setdefault(gene_1, []).append(gene_2)
            edge_list.append((gene_1, gene_2, ppi_value))
    ppi_graph.add_weighted_edges_from(edge_list)
    return ppi_dict, ppi_graph


def load_essential(essential_file):
    e_dict = defaultdict()
    with open(essential_file, 'r') as f:
        for each_line in f.readlines()[1:]:
            e_list = each_line.strip().split('\t')
            e_dict[e_list[2]] = e_list[3]
    return e_dict


def flux_stat(model_file, reaction_name, knockout_dict, organism, gene_dict,
              ppi_dict, taxon_id, cpu, uniprot_dict, meta_dict, essential_dict):
    result_list = []
    model = cobra.io.read_sbml_model(model_file)
    reaction_index = model.reactions.index(reaction_name)
    reaction_dir = os.path.join(output_dir, reaction_name)
    if os.path.exists(reaction_dir):
        shutil.rmtree(reaction_dir)
    os.makedirs(reaction_dir)
    p = Pool(cpu)
    for gene_index in range(len(model.genes)):
        line = p.apply_async(run_cobra, args=(reaction_index, model_file, gene_index))
        result_list.append(line)
    p.close()
    p.join()
    result_file = os.path.join(reaction_dir, '{0}_results.txt'.format(reaction_name))
    valid_gene_list = []
    uniprot_list = []
    print('>>> Loading KEGG pathway')
    with open(result_file, 'w') as f:
        f.write('gene_table\tgene_name\tnormal_knockout_flux\tnormal_reaction_flux\t'
                'reaction_knockout_flux\tdifference\tessentiality\tgene_annotation\t'
                'KEGG\tMetaCyc\tPPI\n')
        for m in result_list:
            m_list = m.get()
            if m_list:
                table = m_list[0]
                k_flux = knockout_dict[table]
                annotation = ''
                name = ''
                meta_line = ''
                if table in gene_dict:
                    annotation = gene_dict[table][1]
                    name = gene_dict[table][0]
                    if meta_dict:
                        if table in meta_dict:
                            meta_line = ', '.join(meta_dict[table])
                n_flux = m_list[1]
                r_flux = m_list[2]
                diff = m_list[3]
                ppi_line = ''
                if table in ppi_dict:
                    ppi_line = ', '.join(ppi_dict[table])
                uniprot_id = ''
                if table in uniprot_dict:
                    uniprot_id = 'UNIPROT:{0}'.format(str(uniprot_dict[table]))
                if float(k_flux) != 0 and float(diff) >= 0.1 * float(n_flux):
                    essentiality = ''
                    if essential_dict:
                        if table in essential_dict:
                            essentiality = essential_dict[table]
                            if essentiality == 'NE':
                                valid_gene_list.append(table)
                                if table in uniprot_dict:
                                    if uniprot_id:
                                        uniprot_list.append(uniprot_id)
                        valid_gene_list.append(table)
                        if table in uniprot_dict:
                            uniprot_list.append(uniprot_id)
                    result_line = ''
                    kegg_pathway_line = load_kegg(table, organism)
                    for i in [table, name, k_flux, n_flux, r_flux, diff, essentiality,
                              annotation, kegg_pathway_line, meta_line]:
                        result_line += '{0}\t'.format(str(i))
                    f.write(result_line.strip('\r|\n') + '\t' + ppi_line + '\n')
    ppi_out_dir = os.path.join(reaction_dir, 'PPI')
    if os.path.exists(ppi_out_dir):
        shutil.rmtree(ppi_out_dir)
    os.makedirs(ppi_out_dir)
    print('>>> Building PPI graph')
    for gene in valid_gene_list:
        get_ppi_image(taxon_id, gene, ppi_out_dir)
    print('>>> Building iPath pathway graph')
    get_ipath_map('\n'.join(uniprot_list), keep_colors=True,
                  reaction_dir=reaction_dir, map_name=reaction_name)


# Run as script
if __name__ == '__main__':
    # Parse command-line
    args = parse_cmdline()
    start_time = datetime.now()
    reactions = args.reactions
    model_file_name = args.model
    thread_num = args.cpu_num
    ppi_file_name = args.ppi
    species_dataset = args.dataset
    metacyc = args.metacyc
    kegg_organism = args.species
    my_taxon = args.taxon
    output = args.output
    essential = args.essential
    # Begin
    my_path = os.getcwd()
    print('> Begin')
    # Print input model file
    print('>> Input model file: {0}'.format(model_file_name))
    my_model_file = os.path.join(my_path, model_file_name)
    print('>> Input reaction(s): {0}'.format(', '.join(reactions)))
    if output is None:
        print('>> Output dir name: all_results')
        output_dir = os.path.join(my_path, 'all_results')
    else:
        print('>> Output dir name: {0}'.format(output))
        output_dir = os.path.join(my_path, output)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)
    # Load all genes annotation
    print('>> Loading BioMart dataset: {0}'.format(species_dataset))
    my_biomart_set = get_biomart(species=species_dataset, meta=metacyc)
    gene_d, uniprot_d, metacyc_d = load_gene_annotation(query_set=my_biomart_set,
                                                        meta=metacyc)
    # Load all PPI data
    print('>> Loading PPI network: {0}'.format(ppi_file_name))
    input_ppi_file = os.path.join(my_path, ppi_file_name)
    input_ppi_dict, input_ppi_graph = load_ppi(input_ppi_file)
    # Load essential genes in DEG database
    essential_gene_dict = None
    if essential:
        print('>> Loading DEG essential genes: {0}'.format(essential))
        essential_gene_file = os.path.join(my_path, essential)
        essential_gene_dict = load_essential(essential_gene_file)
    print('>> Step 1: All genes knockout')
    gene_knockout_dict = essential_gene_knockout(my_model_file)
    print('>> Step 2: gene knockout in inputted reaction(s)')
    for each_reaction in reactions:
        print('>>> Performing reaction {0}'.format(each_reaction))
        flux_stat(model_file=my_model_file,
                  reaction_name=each_reaction,
                  knockout_dict=gene_knockout_dict,
                  gene_dict=gene_d,
                  ppi_dict=input_ppi_dict,
                  taxon_id=my_taxon,
                  cpu=thread_num,
                  organism=kegg_organism,
                  meta_dict=metacyc_d,
                  uniprot_dict=uniprot_d,
                  essential_dict=essential_gene_dict)
        print('>>> done')
    end_time = datetime.now()
    run_time = (end_time - start_time).seconds
    print('> End: {0}s'.format(str(run_time)))
