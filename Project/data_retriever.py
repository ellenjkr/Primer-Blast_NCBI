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


	def __init__(self, primer_params, left_primer, right_primer, saving_path):
		super(DataRetriever, self).__init__()
		self.primer_params = primer_params
		self.left_primer = left_primer
		self.right_primer = right_primer
		self.saving_path = saving_path
		self.database_search = DataBaseSearch(self.default_tax_tags)

	# Realiza a pesquisa e obtem o job key para redirecionar para a página de resultados
	def primerblast_search(self):
		request = requests.post(self.primertool_url, data=self.primer_params)
		parser = BeautifulSoup(request.text, 'html.parser')
		job_key = parser.find("input", {"name":"job_key"}).attrs['value']
		ctg_time = parser.find("input", {"name":"ctg_time"}).attrs['value']
		# ctg_time = None
		# job_key = 'PTfik9Yz25v8pcGgzMDlkrbb9KCbyO-9mg'
		return (job_key, ctg_time)

	def get_results_content(self, job_key):
		print('Aguardando os resultados do primer-blast\n\n')
		# Aguarda os dados carregarem, tenta novamente enquanto não estiverem carregados
		time.sleep(60)
		params = {'job_key': job_key}
		request = requests.get(self.primertool_url, params=params)
		parser = BeautifulSoup(request.text, 'html.parser')
		# classes prPairDtl representam as classificações dos resultados. Ex: Products on target templates
		primerblast_results = parser.find_all("div", class_="prPairDtl") 

		while primerblast_results == []:
			time.sleep(40)
			request = requests.get(self.primertool_url, params=params)
			parser = BeautifulSoup(request.text, 'html.parser')
			primerblast_results = parser.find_all("div", class_="prPairDtl")

		return primerblast_results

	def get_mismatch_info(self, binded_primer, original_primer):
		print(binded_primer)
		print(original_primer)
		mismatched_primer = ''
		mismatches = 0
		mismatched_positions = []

		for pos, i in enumerate(binded_primer):
			if i == '.':
				mismatched_primer += original_primer[pos]
			else:
				mismatched_primer += i
				mismatches += 1
				mismatched_positions.append(pos)

		return (mismatched_primer, mismatches, mismatched_positions)


	def get_binding_info(self, organism_binding):
		print(organism_binding)
		length = re.search(r'(product length = )(\d*)' , organism_binding).group(2)
		# primers = re.findall(r'(Template \s*)(\d*\s*)([\.A-Z]*)' , organism_binding)

		primers_f = re.findall(r'(Forward.*)\n(Template \s*)(\d*\s*)([\.A-Z-]*)' , organism_binding)

		# Em alguns casos ele faz binding apenas do forward ou apenas do reverse
		if primers_f == []:
			primer_f = ""
			template_f = ""
			mismatched_primer_f, mismatches_f, mismatched_positions_f = ("", "", [])
		else:
			# Alguns primers tem traço no primer testado, nao no binding
			forward_primer = primers_f[0][0].strip()
			template_f = primers_f[0][2].strip()
			primer_f = primers_f[0][3].strip()
			# Alguns tem traços, precisa de tratamento, por enquanto ignoramos
			if '-' in primer_f or '-' in forward_primer:
				primer_f = ""
				template_f = ""
				mismatched_primer_f, mismatches_f, mismatched_positions_f = ("", "", [])
			else:
				mismatched_primer_f, mismatches_f, mismatched_positions_f = self.get_mismatch_info(primer_f, self.left_primer)

		primers_r = re.findall(r'(Reverse.*)\n(Template \s*)(\d*\s*)([\.A-Z-]*)' , organism_binding)
		if primers_r == []:
			primer_r = ""
			template_r = ""
			mismatched_primer_r, mismatches_r, mismatched_positions_r = ("", "", [])
		else:
			reverse_primer = primers_r[0][0].strip()
			template_r = primers_r[0][2].strip()
			primer_r = primers_r[0][3].strip()
			if '-' in primer_f or '-' in reverse_primer:
				primer_r = ""
				template_r = ""
				mismatched_primer_r, mismatches_r, mismatched_positions_r = ("", "", [])
			else:
				mismatched_primer_r, mismatches_r, mismatched_positions_r = self.get_mismatch_info(primer_r, self.right_primer)

		
		

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
		unique_organisms = {'Nome': [], 'Ocorrência': []}

		# Lista de resultados para cada organismo
		binding_results = {}

		mismatched_positions = {'f': [], 'r': []}

		for organism in organisms:
			# Procura as duas primeiras palavras contendo apenas letras
			# e que possuam pelo menos 3 e 2 letras, respectivamente
			organism_name = re.search(r"[a-zA-Z]{3,} [a-zA-Z]{2,}", organism.strip()).group(0)
			# Se a segunda palavra tiver menos de 3 letras iremos considerar que o organismo não possui um nome composto
			organism_name = organism_name.split(" ")[0] if len(organism_name.split(" ")[1]) < 3 else organism_name
			# Se for um nome novo, adiciona na lista
			if organism_name not in unique_organisms['Nome']:
				unique_organisms['Nome'].append(organism_name)
				# Contabiliza uma ocorrência
				unique_organisms['Ocorrência'].append(1)

				binding_results[organism_name] = self.get_binding_info(organism)

			# Se for um nome repetido contabiliza mais uma ocorrência
			else:
				unique_organisms['Ocorrência'][unique_organisms['Nome'].index(organism_name)] += 1

				binding_info = self.get_binding_info(organism)
				
				for key in binding_results[organism_name].keys():
					binding_results[organism_name][key].extend(binding_info[key])

		return (unique_organisms, binding_results, mismatched_positions)




	def get_taxonomy(self, organisms_names):
		output_dict = {key: [] for key in self.default_tax_tags}
		output_dict['Nome'] = []
		output_dict['Ocorrência'] = []
		taxonomies = []
		with concurrent.futures.ThreadPoolExecutor(max_workers=12) as executor:
			results = list(executor.map(self.database_search.get_taxonomy, zip(organisms_names['Nome'], organisms_names['Ocorrência'])))

		# Converte o dicionário para um dataframe
		output_df = pd.json_normalize(results)
		first_column = output_df.pop('Nome')
		output_df.insert(0, 'Nome', first_column)

		second_column = output_df.pop('Ocorrência')
		output_df.insert(1, 'Ocorrência', second_column)

		return output_df

	def build_primers_df(self, binding_results):
		primers_df = pd.DataFrame.from_dict(binding_results, orient="index")
		primers_df = primers_df.explode(list(primers_df.columns)).reset_index()
		primers_df = primers_df.rename(columns={"index": "Espécie"})

		return primers_df


	def retrieve_data(self):
		job_key, ctg_time = self.primerblast_search()
		print(f'Os resultados podem ser visualizados em: https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi?job_key={job_key}')
		primerblast_results = self.get_results_content(job_key)
		unique_organisms, binding_results, mismatched_positions = self.get_results(primerblast_results)
		if len(unique_organisms['Nome']) > 0:
			organisms_df = self.get_taxonomy(unique_organisms)
			
			print()
			print('Criando arquivo Excel')

			primers_df = self.build_primers_df(binding_results)
			mismatched_pos_f = primers_df.pop('Mismatches Pos F')
			mismatched_pos_r = primers_df.pop('Mismatches Pos R')

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
						if pos in mismatched_pos_f[index]:
							format_pairs.extend((red, base))
						else:
							format_pairs.extend((green, base))
					
					worksheet.write_rich_string(index + 1, 3, *format_pairs)

				if row['Primer R'] != '':
					format_pairs = []

					# Get each DNA base character from the sequence.
					for pos, base in enumerate(row['Primer R']):
						if pos in mismatched_pos_r[index]:
							format_pairs.extend((red, base))
						else:
							format_pairs.extend((green, base))
					
					worksheet.write_rich_string(index + 1, 6, *format_pairs)

			for pos, col in enumerate(primers_df.columns):
				if col == 'Espécie':
					worksheet.set_column(pos, pos, 25)
		
				elif col in ['Primer F', 'Primer R']:
					worksheet.set_column(pos, pos, 40)
				else:
					worksheet.set_column(pos, pos, 18, centered_cell)
				

			writer.save()

		else:
			print('Nenhum resultado obtido')