import pymysql.cursors


class DataBaseSearch():
	def __init__(self, default_tax_tags):
		super(DataBaseSearch, self).__init__()
		self.default_tax_tags = default_tax_tags
		
	def get_organism_taxid(self, cursor, organism):
		cursor.execute(f"SELECT * FROM names WHERE name_txt LIKE '%{organism}%' AND name_class = 'scientific name' LIMIT 1")
		taxid = cursor.fetchone()['tax_id']

		return taxid


	def get_organism_data(self, cursor, organism):
		cursor.execute(f"SELECT parent_taxnodes_id, _rank FROM nodes WHERE tax_id = {organism} LIMIT 1")
		result = cursor.fetchone()

		return result

	def get_taxonomy(self, arguments):
		connection = pymysql.connect(
			host='localhost',
			user='root',
			password='batata2000',
			database='ncbi_data',
			cursorclass=pymysql.cursors.DictCursor
		)


		with connection:
			with connection.cursor() as cursor:
				try:
					organism_name = arguments[0]
					occurrence = arguments[1]

					print(f"Obtendo a taxonomia de '{organism_name}'")
					organism_taxid = self.get_organism_taxid(cursor, organism_name)

					taxids = []
					taxids.append(organism_taxid)

					taxonomy = {}
					result = self.get_organism_data(cursor, organism_taxid)

					rank = result['_rank']
					taxonomy[rank] = organism_taxid

					while int(result['parent_taxnodes_id']) != 1:
						organism_taxid = result['parent_taxnodes_id']
						taxids.append(organism_taxid)

						result = self.get_organism_data(cursor, organism_taxid)
						rank = result['_rank']
						if rank in self.default_tax_tags:
							taxonomy[rank] = organism_taxid



					cursor.execute(f"SELECT name_txt, tax_id FROM names WHERE tax_id IN {tuple(taxids)} AND name_class = 'scientific name'")
					names = cursor.fetchall()
					
					taxid2name = {list(i.values())[1]: list(i.values())[0] for i in names}
					output_dict = {'Nome': organism_name, 'Ocorrência': occurrence}
					output_dict.update({key: (taxid2name[taxonomy[key]] if key in taxonomy.keys() else None) for key in self.default_tax_tags})
					# output_dict.update({taxonomy[key]: taxid2name[key] for key in taxonomy.keys() if taxonomy[key] in self.self.default_tax_tags})

					return output_dict
				except Exception as e:
					print(f"Não foi possível encontrar a taxonomia de '{organism_name}. Erro: {e}'")