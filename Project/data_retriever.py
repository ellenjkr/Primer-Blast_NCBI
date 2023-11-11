import requests
import time
import re
import pandas as pd
from bs4 import BeautifulSoup

from database_search import DataBaseSearch
import concurrent.futures
import threading

class DataRetriever():
	primertool_url = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi'
	# Classificações taxonômicas obtidas do html da página de taxonomia do NCBI
	default_tax_tags = [
		'species',
		'genus',
		'family',
		'order',
		'class',
		'phylum',
		'kingdom',
		'superkingdom'
	]
	ncbi_taxonomy_url = 'https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi/'


	def __init__(self, primer_params, left_primer, right_primer, saving_path, threads):
		super(DataRetriever, self).__init__()
		self.primer_params = primer_params
		self.left_primer = left_primer
		self.right_primer = right_primer
		self.saving_path = saving_path
		self.threads = threads
		self.database_search = DataBaseSearch(self.default_tax_tags)

	# Realiza a pesquisa e obtem o job key para redirecionar para a página de resultados
	def primerblast_search(self):
		request = requests.post(self.primertool_url, data=self.primer_params)
		parser = BeautifulSoup(request.text, 'html.parser')
		job_key = parser.find("input", {"name":"job_key"}).attrs['value']
		ctg_time = parser.find("input", {"name":"ctg_time"}).attrs['value']
		# ctg_time = None
		# job_key = '190JulZGW-580F7VU7V65ymua9UEvXDIBQ'
		return (job_key, ctg_time)

	def get_results_content(self, job_key):
		print('Aguardando os resultados do primer-blast\n\n')
		# Aguarda os dados carregarem, tenta novamente enquanto não estiverem carregados
		time.sleep(240)
		params = {'job_key': job_key}
		request = requests.get(self.primertool_url, params=params)
		parser = BeautifulSoup(request.text, 'html.parser')
		# classes prPairDtl representam as classificações dos resultados. Ex: Products on target templates
		primerblast_results = parser.find_all("div", class_="prPairDtl") 

		while primerblast_results == []:
			time.sleep(60)
			request = requests.get(self.primertool_url, params=params)
			parser = BeautifulSoup(request.text, 'html.parser')
			primerblast_results = parser.find_all("div", class_="prPairDtl")

		with open(f"{self.saving_path}.html", "w", encoding = 'utf-8') as file: 
		    # prettify the soup object and convert it into a string   
		    file.write(str(parser.prettify()))

		# with open("18_7.html") as fp:
		# 	parser = BeautifulSoup(fp, 'html.parser')

		# primerblast_results = parser.find_all("div", class_="prPairDtl") 
		return primerblast_results

	def get_mismatch_and_gaps_info(self, binded_primer, original_primer, primer_mask):
		# print(binded_primer)
		# print(original_primer)
		binded_primer_nuc = ''
		mismatches = 0
		conflicts_pos = []
		gaps = 0
		gaps_mask = []

		for pos, i in enumerate(primer_mask):
			if i == '-':
				gaps_mask.append(pos)


		for pos, i in enumerate(binded_primer):
			if pos in gaps_mask:
				conflicts_pos.append(pos)
				binded_primer_nuc += i
				mismatches += 1
			else:
				if i == '.':
					binded_primer_nuc += primer_mask[pos]
				elif i != '-':
					binded_primer_nuc += i
					mismatches += 1
					conflicts_pos.append(pos)
				elif i == '-':
					gaps += 1
					binded_primer_nuc += i
					conflicts_pos.append(pos)

		return (binded_primer_nuc, mismatches, conflicts_pos, gaps)



	def get_binding_info(self, organism_binding):

		# length = re.search(r'(product length = )(\d*)' , organism_binding).group(2)
		# # primers = re.findall(r'(Template \s*)(\d*\s*)([\.A-Z]*)' , organism_binding)

		# primers_f = re.findall(r'(Forward.*)\n(Template \s*)(\d*\s*)([\.A-Z-]*)' , organism_binding)

		# # Em alguns casos ele faz binding apenas do forward ou apenas do reverse
		# if primers_f == []:
		# 	primer_f = ""
		# 	template_f = ""
		# 	mismatched_primer_f, mismatches_f, mismatched_positions_f = ("", "", [])
		# else:
		# 	# Alguns primers tem traço no primer testado, nao no binding
		# 	forward_primer = primers_f[0][0].strip()
		# 	template_f = primers_f[0][2].strip()
		# 	primer_f = primers_f[0][3].strip()
		# 	# Alguns tem traços, precisa de tratamento, por enquanto ignoramos
		# 	if '-' in primer_f or '-' in forward_primer:
		# 		primer_f = ""
		# 		template_f = ""
		# 		mismatched_primer_f, mismatches_f, mismatched_positions_f = ("", "", [])
		# 	else:
		# 		mismatched_primer_f, mismatches_f, mismatched_positions_f = self.get_mismatch_info(primer_f, self.left_primer)

		# primers_r = re.findall(r'(Reverse.*)\n(Template \s*)(\d*\s*)([\.A-Z-]*)' , organism_binding)
		# if primers_r == []:
		# 	primer_r = ""
		# 	template_r = ""
		# 	mismatched_primer_r, mismatches_r, mismatched_positions_r = ("", "", [])
		# else:
		# 	reverse_primer = primers_r[0][0].strip()
		# 	template_r = primers_r[0][2].strip()
		# 	primer_r = primers_r[0][3].strip()
		# 	if '-' in primer_f or '-' in reverse_primer:
		# 		primer_r = ""
		# 		template_r = ""
		# 		mismatched_primer_r, mismatches_r, mismatched_positions_r = ("", "", [])
		# 	else:
		# 		mismatched_primer_r, mismatches_r, mismatched_positions_r = self.get_mismatch_info(primer_r, self.right_primer)

		
		

		# return {
		# 	'Tamanho': [length.strip()],
		# 	'Template F': [template_f],
		# 	'Primer F': [mismatched_primer_f],
		# 	'Mismatches F': [mismatches_f],
		# 	'Mismatches Pos F': [mismatched_positions_f],
		# 	'Template R': [template_r],
		# 	'Primer R': [mismatched_primer_r],
		# 	'Mismatches R': [mismatches_r],
		# 	'Mismatches Pos R': [mismatched_positions_r]
		# }

		valid_products = primers_f = re.findall(r'(product length = )(\d*)\n(Forward primer\s*\d*\s*)([\.A-Z-]*)\s*\d*\n(Template \s*)(\d*\s*)([\.A-Z-]*)(\s*\d*)\n\n(Reverse primer\s*\d*\s*)([\.A-Z-]*)\s*\d*\n(Template \s*)(\d*\s*)([\.A-Z-]*)' , organism_binding)
		if valid_products != []:
			valid_product = valid_products[0]
			length = valid_product[1]

			primers_f = re.findall(r'(Forward.*)\n(Template \s*)(\d*\s*)([\.A-Z-]*)' , organism_binding)

			forward_primer = valid_product[3]
			template_f = valid_product[5].strip()
			primer_f = valid_product[6]

			binded_primer_f, mismatches_f, conflicts_pos_f, gaps_f = self.get_mismatch_and_gaps_info(primer_f, self.left_primer, forward_primer)

			reverse_primer = valid_product[9]
			template_r = valid_product[11].strip()
			primer_r = valid_product[12]

			binded_primer_r, mismatches_r, conflicts_pos_r, gaps_r = self.get_mismatch_and_gaps_info(primer_r, self.right_primer, reverse_primer)

			if mismatches_r + gaps_r > 2 or mismatches_f + gaps_f > 2:
				return None
			return {
				'Tamanho': [length.strip()],
				'Template F': [template_f],
				'Primer F': [binded_primer_f],
				'Mismatches F': [mismatches_f],
				'Gaps F': [gaps_f],
				'Conflicts Pos F': [conflicts_pos_f],
				'Template R': [template_r],
				'Primer R': [binded_primer_r],
				'Mismatches R': [mismatches_r],
				'Gaps R': [gaps_r],
				'Conflicts Pos R': [conflicts_pos_r]
			}
		else:
			return None


	def get_results(self, primerblast_results):
		# Obtém organismos
		organisms = []
		for div in primerblast_results:
			if div.text.strip() != '':
				# Cada organismo é separado por um título que começa com ">"
				organisms += div.text.strip().split('>')
		
		if '' in organisms:
			organisms.remove('') # Remove espaços vazios
		# Lista de organismos
		unique_organisms = {'Nome': [], 'OcorrênciaAcc': []}

		# Lista de resultados para cada organismo
		binding_results = {}

		mismatched_positions = {'f': [], 'r': []}

		# for organism in organisms:
		# 	# Procura as duas primeiras palavras contendo apenas letras
		# 	# e que possuam pelo menos 3 e 2 letras, respectivamente
		# 	organism_name = re.search(r"[a-zA-Z]{3,} [a-zA-Z]{2,}", organism.strip()).group(0)
		# 	# Se a segunda palavra tiver menos de 3 letras iremos considerar que o organismo não possui um nome composto
		# 	organism_name = organism_name.split(" ")[0] if len(organism_name.split(" ")[1]) < 3 else organism_name
		# 	# Se for um nome novo, adiciona na lista
		# 	if organism_name not in unique_organisms['Nome']:
		# 		unique_organisms['Nome'].append(organism_name)
		# 		# Contabiliza uma ocorrência
		# 		unique_organisms['Ocorrência'].append(1)

		# 		binding_results[organism_name] = self.get_binding_info(organism)

		# 	# Se for um nome repetido contabiliza mais uma ocorrência
		# 	else:
		# 		unique_organisms['Ocorrência'][unique_organisms['Nome'].index(organism_name)] += 1

		# 		binding_info = self.get_binding_info(organism)
				
		# 		for key in binding_results[organism_name].keys():
		# 			binding_results[organism_name][key].extend(binding_info[key])

		for organism in organisms:
			organism_acc = organism.strip().split('\n')[0].split('.')[0]

			if organism_acc not in unique_organisms['Nome']:
				# unique_organisms['Nome'].append(organism_acc)
				# # Contabiliza uma ocorrência
				# unique_organisms['OcorrênciaAcc'].append(1)

				# binding_results[organism_acc] = self.get_binding_info(organism)
				binding_info =  self.get_binding_info(organism)
				if binding_info is not None:
					unique_organisms['Nome'].append(organism_acc)
					# Contabiliza uma ocorrência
					unique_organisms['OcorrênciaAcc'].append(1)

					binding_results[organism_acc] = binding_info

			# Se for um nome repetido contabiliza mais uma ocorrência
			else:
				binding_info = self.get_binding_info(organism)
				if binding_info is not None:
					unique_organisms['OcorrênciaAcc'][unique_organisms['Nome'].index(organism_acc)] += 1

					
					for key in binding_results[organism_acc].keys():
						binding_results[organism_acc][key].extend(binding_info[key])
		with concurrent.futures.ThreadPoolExecutor(max_workers=self.threads) as executor:
			results = list(executor.map(self.database_search.get_organism_taxid, zip(unique_organisms['Nome'], unique_organisms['OcorrênciaAcc'])))

		acc_report = pd.json_normalize(results)
		taxid_count = acc_report['TaxId'].value_counts().to_frame().reset_index()


		acc_report = acc_report.merge(taxid_count, on='TaxId', how='left')
		acc_report = acc_report.rename(columns={'count': 'Ocorrência'})
		acc_report['Ocorrência'] = acc_report['OcorrênciaAcc'] * acc_report['Ocorrência']
		taxid_count = acc_report[['TaxId', 'Ocorrência']].drop_duplicates(subset='TaxId').reset_index(drop=True)



		index_mapping = dict(zip(acc_report['Nome'], acc_report['TaxId']))

		binding_results_df = pd.DataFrame.from_dict(binding_results, orient="index")
		binding_results_df = binding_results_df.rename(index=index_mapping)

		return (acc_report, taxid_count, binding_results_df, mismatched_positions)
		# return (unique_organisms, binding_results, mismatched_positions)




	def get_taxonomy(self, organisms):
		print('Obtendo as taxonomias dos resultados')
		output_dict = {key: [] for key in self.default_tax_tags}
		output_dict['TaxId'] = []
		output_dict['Ocorrência'] = []
		taxonomies = []
		with concurrent.futures.ThreadPoolExecutor(max_workers=self.threads) as executor:
			results = list(executor.map(self.database_search.get_taxonomy, zip(organisms['TaxId'], organisms['Ocorrência'])))

		# Converte o dicionário para um dataframe
		output_df = pd.json_normalize(results)
		first_column = output_df.pop('TaxId')
		output_df.insert(0, 'TaxId', first_column)

		second_column = output_df.pop('Ocorrência')
		output_df.insert(1, 'Ocorrência', second_column)

		output_df.drop(columns='species', inplace=True)
		output_df.rename(columns={'tax_name': 'species'}, inplace=True)

		return output_df

	def build_primers_df(self, binding_results):
		# primers_df = pd.DataFrame.from_dict(binding_results, orient="index")
		primers_df = binding_results.explode(list(binding_results.columns)).reset_index()
		primers_df = primers_df.rename(columns={"index": "TaxId"})

		return primers_df


	def retrieve_data(self):
		job_key, ctg_time = self.primerblast_search()
		print(f'Os resultados podem ser visualizados em: https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?job_key={job_key}')
		primerblast_results = self.get_results_content(job_key)
		# unique_organisms, binding_results, mismatched_positions = self.get_results(primerblast_results)
		acc_report, taxid_count, binding_results_df, mismatched_positions = self.get_results(primerblast_results)

		# if len(unique_organisms['Nome']) > 0:
		# 	organisms_df = self.get_taxonomy(unique_organisms)
		if len(taxid_count['TaxId']) > 0:
			organisms_df = self.get_taxonomy(taxid_count)
			organisms_df.drop(['tax_id'], axis=1)
			print()
			print('Criando arquivo Excel')

			primers_df = self.build_primers_df(binding_results_df)
			conflicts_pos_f = primers_df.pop('Conflicts Pos F')
			conflicts_pos_r = primers_df.pop('Conflicts Pos R')

			writer = pd.ExcelWriter(f'{self.saving_path}.xlsx', engine='xlsxwriter')
	
			organisms_df.to_excel(writer, sheet_name='Espécies', startrow=1, header=False, index=False)
			primers_df.to_excel(writer, sheet_name='Primers', startrow=1, header=False, index=False)

			workbook = writer.book
			centered_cell = workbook.add_format({'align': 'center'})

			worksheet = writer.sheets['Espécies']
			(max_row, max_col) = organisms_df.shape
			column_settings = [{'header': column} for column in organisms_df.columns]
			worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings, 'style': 'Table Style Medium 9'})

			for pos, col in enumerate(organisms_df.columns):
				if col == 'Nome': 
					worksheet.set_column(pos, pos, 25)
				elif col == 'Ocorrência':
					worksheet.set_column(pos, pos, 18, centered_cell)
				else:
					worksheet.set_column(pos, pos, 20)



			worksheet = writer.sheets['Primers']
			(max_row, max_col) = primers_df.shape
			column_settings = [{'header': column} for column in primers_df.columns]
			worksheet.add_table(0, 0, max_row, max_col - 1, {'columns': column_settings, 'style': 'Table Style Medium 9'})


			red = workbook.add_format({'color': 'red'})
			green = workbook.add_format({'color': 'green'})

			for index, row in primers_df.iterrows():
				if row['Primer F'] != '':
					format_pairs = []

					# Get each DNA base character from the sequence.
					for pos, base in enumerate(row['Primer F']):
						if pos in conflicts_pos_f[index]:
							format_pairs.extend((red, base))
						else:
							format_pairs.extend((green, base))
					
					worksheet.write_rich_string(index + 1, 3, *format_pairs)

				if row['Primer R'] != '':
					format_pairs = []

					# Get each DNA base character from the sequence.
					for pos, base in enumerate(row['Primer R']):
						if pos in conflicts_pos_r[index]:
							format_pairs.extend((red, base))
						else:
							format_pairs.extend((green, base))
					
					worksheet.write_rich_string(index + 1, 7, *format_pairs)

			for pos, col in enumerate(primers_df.columns):
				if col in ['Primer F', 'Primer R']:
					worksheet.set_column(pos, pos, 40)
				else:
					worksheet.set_column(pos, pos, 18, centered_cell)
				

			writer.close()

		else:
			print('Nenhum resultado obtido')