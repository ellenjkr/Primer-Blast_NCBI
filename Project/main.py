import requests
import time
import re
import pandas as pd
import os
import pymysql.cursors

from bs4 import BeautifulSoup
from itertools import product

from data_retriever import DataRetriever


def extend_ambiguous_dna(seq):
	degenerated_table = {
		'A': ['A'],
		'C': ['C'],
		'G': ['G'],
		'T': ['T'],
		'W': ['A', 'T'],
		'S': ['C', 'G'],
		'M': ['A', 'C'],
		'K': ['G', 'T'],
		'R': ['A', 'G'],
		'Y': ['C', 'T'],
		'B': ['C', 'G', 'T'],
		'D': ['A', 'G', 'T'],
		'H': ['A', 'C', 'T'],
		'V': ['A', 'C', 'G'],
		'N': ['A', 'C', 'G', 'T'],
	}

	r = []

	# Aplica produto cartesiano nos conjuntos (conjunto = possíveis bases para cada letra)
	for i in product(*[degenerated_table[j] for j in seq]):
		r.append("".join(i))

	return r


def sep_degenerated(path):
	for input_file in os.listdir(path):
		if '.csv' in input_file:
			primer_set = input_file.replace('.csv', '')
			primers_df = pd.read_csv(f"{path}/{input_file}", sep=';')
			with open(f"{path}/{input_file}", "w") as non_degenerated:
				non_degenerated.write("NAME;LEFT_PRIMER;RIGHT_PRIMER\n")
			
			for index, row in primers_df.iterrows():
				left_results = extend_ambiguous_dna(row['LEFT_PRIMER'])
				right_results = extend_ambiguous_dna(row['RIGHT_PRIMER'])
				all_possible_pairs = product(*[left_results, right_results])

				all_possible_pairs = [i for i in all_possible_pairs]

				if len(all_possible_pairs) == 1:
					with open(f"{path}/{input_file}", "a") as non_degenerated:
						non_degenerated.write(f"{row['NAME']};{row['LEFT_PRIMER']};{row['RIGHT_PRIMER']}\n")
				else:
					with open(f"{path}/degenerated_{row['NAME']}.csv", "w") as f:
						f.write("NAME;LEFT_PRIMER;RIGHT_PRIMER\n")
						for pos, item in enumerate(all_possible_pairs):
							f.write(f"{pos + 1};{item[0]};{item[1]}\n")



def run_from_input(primers_dir):

	sep_degenerated(primers_dir)

	for input_file in os.listdir(primers_dir):
		if '.csv' in input_file:
			primer_set = input_file.replace('.csv', '')
			primers_df = pd.read_csv(f"{primers_dir}/{input_file}", sep=';')

			if os.path.exists(f'{primers_dir}/{primer_set}') is False:
				os.mkdir(f'{primers_dir}/{primer_set}')

				if 'degenerated' in input_file and primers_df.shape[0] > 5:
					primers_df = primers_df.sample(n=10).reset_index()

				for index, row in primers_df.iterrows():

					left_primer = str(row.LEFT_PRIMER).strip()
					right_primer = str(row.RIGHT_PRIMER).strip()

					if len(left_primer) < 36 and len(right_primer) < 36:
						
						print('Par de primers:', input_file, row['NAME'])
						# Realiza a pesquisa e obtem o job key para redirecionar para a página de resultados
						parameters = {
							'PRIMER_LEFT_INPUT': left_primer,
							'PRIMER_RIGHT_INPUT': right_primer,
							'PRIMER_SPECIFICITY_DATABASE': 'nt', # taxids_sg
							'SEARCH_SPECIFIC_PRIMER': True,
							'PRIMER_PRODUCT_MIN': 70,
							'PRIMER_PRODUCT_MAX': 1000,
							'PRIMER_NUM_RETURN': 10,
							'PRIMER_MIN_TM': 57.0,
							'PRIMER_OPT_TM': 60.0,
							'PRIMER_MAX_TM': 63.0,
							'PRIMER_MAX_DIFF_TM': 3,
							'PRIMER_ON_SPLICE_SITE': 0,
							'SPLICE_SITE_OVERLAP_5END': 7,
							'SPLICE_SITE_OVERLAP_3END': 4,
							'SPLICE_SITE_OVERLAP_3END_MAX': 8,
							'MIN_INTRON_SIZE': 1000,
							'MAX_INTRON_SIZE': 1000000,
							'SEARCH_SPECIFIC_PRIMER': 'checked',
							'SEARCHMODE': 0,
							'PRIMER_SPECIFICITY_MISMATCH': 1,
							'PRIMER_3END_SPECIFICITY_MISMATCH': 1,
							'MISMATCH_REGION_LENGTH': 5,
							'TOTAL_MISMATCH_IGNORE': 6,
							'MAX_TARGET_SIZE': 4000,
							'HITSIZE': 100000,
							'EVALUE': 30000,
							'WORD_SIZE': 7,
							'MAX_CANDIDATE_PRIMER': 500,
							'NUM_TARGETS': 20,
							'NUM_TARGETS_WITH_PRIMERS': 1000,
							'MAX_TARGET_PER_TEMPLATE': 100000,
							'PRIMER_MIN_SIZE': 15,
							'PRIMER_OPT_SIZE': 20,
							'PRIMER_MAX_SIZE': 25,
							'PRIMER_MIN_GC': 20.0,
							'PRIMER_MAX_GC': 80.0,
							'GC_CLAMP': 0,
							'POLYX': 5,
							'PRIMER_MAX_END_STABILITY': 9,
							'PRIMER_MAX_END_GC': 5,
							'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': 40.00,
							'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH': 70.00,
							'PRIMER_MAX_SELF_ANY_TH': 45.0,
							'PRIMER_MAX_SELF_END_TH': 35.0,
							'PRIMER_PAIR_MAX_COMPL_ANY_TH': 45.0,
							'PRIMER_PAIR_MAX_COMPL_END_TH': 35.0,
							'PRIMER_MAX_HAIRPIN_TH': 24.0,
							'PRIMER_MAX_TEMPLATE_MISPRIMING': 12.00,
							'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING': 24.00,
							'SELF_ANY': 8.00,
							'SELF_END': 3.00,
							'PRIMER_PAIR_MAX_COMPL_ANY': 8.00,
							'PRIMER_PAIR_MAX_COMPL_END': 3.00,
							'OVERLAP_5END': 7,
							'OVERLAP_3END': 4,
							'MONO_CATIONS': 50.0,
							'DIVA_CATIONS': 1.5,
							'CON_DNTPS': 0.6,
							'SALT_FORMULAR': 1,
							'TM_METHOD': 1,
							'CON_ANEAL_OLIGO': 50.0,
							'PRIMER_MISPRIMING_LIBRARY': 'AUTO',
							'LOW_COMPLEXITY_FILTER': 'checked',
							'PRIMER_INTERNAL_OLIGO_MIN_SIZE': 18,
							'PRIMER_INTERNAL_OLIGO_OPT_SIZE': 20,
							'PRIMER_INTERNAL_OLIGO_MAX_SIZE': 27,
							'PRIMER_INTERNAL_OLIGO_MIN_TM': 57.0,
							'PRIMER_INTERNAL_OLIGO_OPT_TM': 60.0,
							'PRIMER_INTERNAL_OLIGO_MAX_TM': 63.0,
							'PRIMER_INTERNAL_OLIGO_MIN_GC': 20.0,
							'PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT': 50,
							'PRIMER_INTERNAL_OLIGO_MAX_GC': 80.0,
						}
						#forward:
							#length = 20
							#Tm = 56.92
							#GC 50.00
							#selft complementarity 7
							#Self 3' complementarity 5

						print('Forward primer:', left_primer)
						print('Reverse primer:', right_primer)
						data_retriever = DataRetriever(parameters, left_primer, right_primer, f"{primers_dir}/{primer_set}/{row['NAME']}")
						data_retriever.retrieve_data()
						print()


if __name__ == '__main__':
	primers_dir = 'Primers_Plantas'
	run_from_input(primers_dir)

# tem uns primers com um traço na frente, o que é esse traço?
# isso aqui também acontece:
# Reverse primer  1        CACCTCA-GGGTGTCCGAAGAACCAGAA  27
# Template        1483744  G......G.....G..............  1483771

# Verificar primers que vieram vazios
# Eliminar os que so tem o reverse primer?
# Obter identificação pelo accession id e nao pelo nome
# https://realpython.com/python-concurrency/
# https://www.toptal.com/python/beginners-guide-to-concurrency-and-parallelism-in-python
# https://www.infoworld.com/article/3632284/python-concurrency-and-parallelism-explained.html
# https://docs.python.org/3/library/multiprocessing.html
# https://www.section.io/engineering-education/mysql-query-optimization-using-indexes-with-examples/