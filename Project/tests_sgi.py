import requests
import re
import pandas as pd
from bs4 import BeautifulSoup
import subprocess
import pymysql.cursors
import concurrent.futures


# # set the url to perform the get request 
# URL = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?ctg_time=1697126947&job_key=QEqeLsvdxnXhT1ZKWypyeCExY0oMInhXDQ'
# page = requests.get(URL) 

# # load the page content 
# text = page.content 

# # make a soup object by using beautiful 
# # soup and set the markup as html parser 
# soup = BeautifulSoup(text, "html.parser") 

# # open the file in w mode 
# # set encoding to UTF-8 
# with open("aux3.html", "w", encoding = 'utf-8') as file: 
	
# 	# prettify the soup object and convert it into a string 
# 	file.write(str(soup.prettify()))



with open("aux3.html") as fp:
	parser = BeautifulSoup(fp, 'html.parser')

primerblast_results = parser.find_all("div", class_="prPairDtl") 


# Obtém organismos
organisms = []
for div in primerblast_results:
	if div.text.strip() != '':
		# Cada organismo é separado por um título que começa com ">"
		organisms.extend(div.text.strip().split('>'))

if '' in organisms:
	organisms.remove('') # Remove espaços vazios

# Lista de organismos
unique_organisms = {'ACC': [], 'OcorrênciaAcc': []}

# Lista de resultados para cada organismo
binding_results = {}

mismatched_positions = {'f': [], 'r': []}

def get_mismatch_info(binded_primer, original_primer):
	# print(binded_primer)
	# print(original_primer)
	mismatched_primer = ''
	mismatches = 0
	mismatched_positions = []

	for pos, i in enumerate(binded_primer):
		if i == '.':
			mismatched_primer += original_primer[pos]
		elif i != '-':
			mismatched_primer += i
			mismatches += 1
			mismatched_positions.append(pos)

	return (mismatched_primer, mismatches, mismatched_positions)


def get_gap_info(binded_primer, original_primer):
	# print(binded_primer)
	# print(original_primer)
	gapped_primer = ''
	gaps = 0
	gaps_positions = []

	for pos, i in enumerate(binded_primer):
		if i == '.' or i != '-':
			gapped_primer += original_primer[pos]
		else:
			gapped_primer += i
			gaps += 1
			gaps_positions.append(pos)

	return (gapped_primer, gaps, gaps_positions)

def get_binding_info(organism_binding):	
	valid_products = primers_f = re.findall(r'(product length = )(\d*)\n(Forward primer\s*\d*\s*)([\.A-Z-]*)\s*\d*\n(Template \s*)(\d*\s*)([\.A-Z-]*)(\s*\d*)\n\n(Reverse primer\s*\d*\s*)([\.A-Z-]*)\s*\d*\n(Template \s*)(\d*\s*)([\.A-Z-]*)' , organism_binding)
	if valid_products != []:
		valid_product = valid_products[0]
		length = valid_product[1]

		primers_f = re.findall(r'(Forward.*)\n(Template \s*)(\d*\s*)([\.A-Z-]*)' , organism_binding)

		forward_primer = valid_product[3]
		template_f = valid_product[5].strip()
		primer_f = valid_product[6]
		# Alguns tem traços, precisa de tratamento, por enquanto ignoramos
		if '-' in forward_primer:
			print('"-" in forward primer')


		mismatched_primer_f, mismatches_f, mismatched_positions_f = get_mismatch_info(primer_f, 'CGAATAAATAATATAAGATTTTG')
		gapped_primer_f, gaps_f, gaps_positions_f = get_gap_info(primer_f, 'CGAATAAATAATATAAGATTTTG')
		print(gapped_primer_f, gaps_f, gaps_positions_f)

		reverse_primer = valid_product[9]
		template_r = valid_product[11].strip()
		primer_r = valid_product[12]

		if '-' in reverse_primer:
			print('"-" in reverse primer')
		mismatched_primer_r, mismatches_r, mismatched_positions_r = get_mismatch_info(primer_r, 'CCATCTAAAAATTTTAATTCCAGT')
		# print(primer_r)
		# print(mismatches_r)

		return {
			'Tamanho': [length.strip()],
			'Template F': [template_f],
			'Primer F': [mismatched_primer_f],
			'Mismatches F': [mismatches_f],
			'Mismatches Pos F': [mismatched_positions_f],
			'Template R': [template_r],
			'Primer R': [mismatched_primer_r],
			'Mismatches R': [mismatches_r],
			'Mismatches Pos R': [mismatched_positions_r]
		}
	else:
		return None
	


def get_taxid(arguments):
	connection = pymysql.connect(
		host='localhost',
		user='root',
		password='Amora#1000',
		database='ncbi_data',
		cursorclass=pymysql.cursors.DictCursor
	)


	with connection:
		with connection.cursor() as cursor:
			try:
				acc = arguments[0]
				occurrence = arguments[1]

				cursor.execute(f"SELECT TaxId FROM accession2taxid WHERE Accession = '{acc}' LIMIT 1")
				result = cursor.fetchone()
				if result is None:
					taxid = str(subprocess.run(
						f'esearch -db nucleotide -query "{acc}" | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId',
						shell=True,
						executable='/bin/bash',
						capture_output=True
					).stdout)
					if taxid != "b''":
						taxid = re.findall(r'\d+', taxid)[0]
						result = {'TaxId': taxid}
					else:
						result = {'TaxId': 'unclassified'}

				output_dict = {'ACC': acc, 'OcorrênciaAcc': occurrence}

				output_dict.update(result)

				return output_dict
			except Exception as e:
				print(f"Não foi possível encontrar a taxonomia de '{acc}'. Erro: {e}")


for organism in organisms:
	organism_acc = organism.strip().split('\n')[0].split('.')[0]

	if organism_acc not in unique_organisms['ACC']:
		binding_info = get_binding_info(organism)
		if binding_info is not None:
			unique_organisms['ACC'].append(organism_acc)
			# Contabiliza uma ocorrência
			unique_organisms['OcorrênciaAcc'].append(1)

			binding_results[organism_acc] = binding_info

	# Se for um nome repetido contabiliza mais uma ocorrência
	else:
		binding_info = get_binding_info(organism)
		if binding_info is not None:
			unique_organisms['OcorrênciaAcc'][unique_organisms['ACC'].index(organism_acc)] += 1

			for key in binding_results[organism_acc].keys():
				binding_results[organism_acc][key].extend(binding_info[key])


with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
	results = list(executor.map(get_taxid, zip(unique_organisms['ACC'], unique_organisms['OcorrênciaAcc'])))

output_df = pd.json_normalize(results)
taxid_count = output_df['TaxId'].value_counts().to_frame().reset_index()


output_df = output_df.merge(taxid_count, on='TaxId', how='left')
output_df = output_df.rename(columns={'count': 'Ocorrência'})
output_df['Ocorrência'] = output_df['OcorrênciaAcc'] * output_df['Ocorrência']
taxid_count = output_df[['TaxId', 'Ocorrência']].drop_duplicates(subset='TaxId').reset_index(drop=True)


print(output_df)
print(taxid_count)

index_mapping = dict(zip(output_df['ACC'], output_df['TaxId']))
binding_results_df = pd.DataFrame.from_dict(binding_results, orient="index")
# binding_results_df = binding_results_df.rename(index=index_mapping)

primers_df = binding_results_df.explode(list(binding_results_df.columns)).reset_index()
primers_df = primers_df.rename(columns={"index": "TaxId"})
primers_df.to_csv('primers_df.tsv', sep='\t', index=False)
print(primers_df)

