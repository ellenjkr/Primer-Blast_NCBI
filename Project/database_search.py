import pymysql.cursors
import subprocess
import re


class DataBaseSearch():
	def __init__(self, default_tax_tags):
		super(DataBaseSearch, self).__init__()
		self.default_tax_tags = default_tax_tags

	# def get_organism_taxid(self, cursor, organism):
	# 	cursor.execute(f"SELECT * FROM names WHERE name_txt = '{organism}' LIMIT 1")
	# 	result = cursor.fetchone()
	# 	if result is None:
	# 		cursor.execute(f"SELECT * FROM names WHERE MATCH(name_txt) AGAINST('{organism}' IN NATURAL LANGUAGE MODE)")
	# 		result = cursor.fetchone()

	# 	taxid = result['tax_id']
	# 	name_txt = result['name_txt']

	# 	return taxid, name_txt


	def get_organism_taxid(self, arguments):
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

					output_dict = {'Nome': acc, 'OcorrênciaAcc': occurrence}

					output_dict.update(result)

					return output_dict
				except Exception as e:
					print(f"Não foi possível encontrar a taxonomia de '{acc}'. Erro: {e}")


	def get_organism_data(self, cursor, organism):
		cursor.execute(f"SELECT parent_taxnodes_id, _rank FROM nodes WHERE tax_id = {organism} LIMIT 1")
		result = cursor.fetchone()

		return result

	def get_taxonomy(self, arguments):
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
					# organism_name = arguments[0]
					organism_taxid = arguments[0]
					occurrence = arguments[1]

					# print(f"Obtendo a taxonomia de '{organism_name}'")
					# organism_taxid, name_txt = self.get_organism_taxid(cursor, organism_name)
					cursor.execute(f"SELECT * FROM organismsTest WHERE tax_id = {organism_taxid} LIMIT 1")
					result = cursor.fetchone()
					
					output_dict = {'TaxId': organism_taxid, 'Ocorrência': occurrence}

					output_dict.update(result)
					# output_dict.update({taxonomy[key]: taxid2name[key] for key in taxonomy.keys() if taxonomy[key] in self.self.default_tax_tags})

					return output_dict
				except Exception as e:
					print(f"Não foi possível encontrar a taxonomia de '{organism_taxid}. Erro: {e}'")