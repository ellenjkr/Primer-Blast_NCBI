import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from pathlib import Path


def get_degenerated_df(dfs):
	columns = dfs[0].columns.to_list()
	columns_operations = {}
	for col in columns:
		if col != 'TaxId' and col != 'Ocorrência':
			columns_operations[col] = 'first'
		elif col == 'Ocorrência':
			columns_operations[col] = 'mean'

	df = pd.concat(dfs).groupby('TaxId', as_index=False).agg(columns_operations)
	df['Ocorrência'] = df['Ocorrência'].astype(int)

	return df


def parse_data(main_path, tax_level):
	df = None
	# data = []
	# keys = []

	# for item in os.listdir('09112023'):
	# 	if os.path.isdir(f'09112023/{item}'):
	# 		for file in os.listdir(f'09112023/{item}'):
	# 			if 'xlsx' in file:
	# 				df = pd.read_excel(f'09112023/{item}/{file}', sheet_name='Espécies')
	# 				df_counts = df[['species', 'Ocorrência']].groupby('species').sum().reset_index()
	# 				# df_counts = df['species'].value_counts().reset_index(name='counts')
	# 				data.append(df_counts)
	# 				key = file.split('.')[0]
	# 				keys.append(key)




	# df = pd.concat(data, keys=keys)

	# df.reset_index(level=0, inplace=True)
	# df.rename(columns={'level_0': 'Primer'}, inplace=True)


	# pivot = pd.pivot_table(data=df, index=['Primer'], columns=['species'], values='Ocorrência')

	data = []
	keys = []

	degenerated_paths = []
	for path in Path(main_path).rglob('*.xlsx'):
		if 'degenerated' in str(path):
			if path.parent not in degenerated_paths:
				degenerated_paths.append(path.parent)
		else:
			df = pd.read_excel(path, sheet_name='Espécies')
			df_counts = df[[tax_level, 'Ocorrência']].groupby(tax_level).sum().reset_index()
			# df_counts = df['species'].value_counts().reset_index(name='counts')
			data.append(df_counts)
			key = path.name.split('.')[0]
			keys.append(key)

	for path in degenerated_paths:
		dfs = [pd.read_excel(path.joinpath(file), sheet_name='Espécies') for file in os.listdir(path) if 'xlsx' in file]
		df = get_degenerated_df(dfs)
		df_counts = df[[tax_level, 'Ocorrência']].groupby(tax_level).sum().reset_index()
		# df_counts = df['phylum'].value_counts().reset_index(name='counts')
		data.append(df_counts)
		keys.append(path.name.replace('degenerated_', ''))

	df = pd.concat(data, keys=keys)

	df.reset_index(level=0, inplace=True)
	df.rename(columns={'level_0': 'Primer'}, inplace=True)


	pivot = pd.pivot_table(data=df, index=['Primer'], columns=[tax_level], values='Ocorrência')

	return pivot
	

	
# Função para normalizar os valores em 100
def normalize_to_100(row):
	total = row.sum()  # Soma os valores na linha
	return (row / total) * 100


def comparison_chart(pivot, tax_level, normalized=False):
	df = deepcopy(pivot)
	if normalized:
		df = df.apply(normalize_to_100, axis=1)

	all_primers_top10 = []
	# print(df)

	for index, row in df.iterrows():
		top_10 = list(row.sort_values(ascending=False).head(10).index)
		
		for i in top_10:
			if i not in all_primers_top10:
				all_primers_top10.append(i)

	others = []
	for col in df.columns:
		if col not in all_primers_top10:
			others.append(col)
	

	if len(others) > 0:
		df['Outros'] = df[others].sum(axis=1)
		
	df = df.drop(columns=others)

	ax = df.plot.bar(stacked=True, legend=False, figsize=(10, 6), color=['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#FF0000', '#6AFF00', '#B9EDE0'], width=0.8)
	
	plt.xticks(rotation=45)
	if normalized:
		ax.set_ylim(top=100)

	ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
	plt.margins(0.2)
	plt.subplots_adjust(bottom=0.15, right=0.75)


	plt.title(f'Amplificações por Par de Primers: TOP 10 {tax_level}')
	plt.xlabel('Par de Primer')
	plt.ylabel('Organismos')


	plt.show()


def autolabel(ax):
	total = int(sum(rect.get_height() for rect in ax.patches))
	rect = ax.patches[0]
	height = rect.get_height()
	ax.annotate('{}'.format(total),
				xy=(rect.get_x() + rect.get_width() / 2, total),
				xytext=(0, 3),  # 3 points vertical offset
				textcoords="offset points",
				ha='center', va='bottom')


def generate_single_chart(row, tax_level):
	row = row.dropna()
	top_5 = list(row.sort_values(ascending=False).head(5).index)

	others = [idx for idx in row.index if idx not in top_5]

	if len(others) > 0:
		row.at['Outros'] = row.loc[others].sum()

	row = row.drop(others)
	row = row.sort_values(ascending=False)

	ax = pd.DataFrame(row).T.plot.bar(stacked=True, color=['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000', '#FF0000', '#FFFF00', '#6AFF00', '#B9EDE0'], width=0.8)

	plt.xticks(rotation=0)
	ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
	plt.margins(0.2)
	plt.subplots_adjust(bottom=0.15, right=0.5, left=0.20)
	
	plt.title(f'Amplificações {row.name}: TOP 5 {tax_level}', pad=15)
	plt.xlabel('Par de Primer', labelpad=6.0)
	plt.ylabel('Organismos', labelpad=6.0)

	autolabel(ax)

	plt.show()



def run(config):
	path = config['WORKING_PATH']
	tax_levels = [tax_level for tax_level in config['CHARTS'].keys() if config['CHARTS'][tax_level] is True]

	tax_levels_transl = {'species': 'Espécies', 'genus': 'Gêneros', 'family': 'Famílias', 'order': 'Ordens', 'class': 'Classes', 'phylum': 'Filos', 'kingdom': 'Reinos', 'superkingdom': 'Super-Reinos'}

	for tax_level in tax_levels:
		data = parse_data(path, tax_level)
		comparison_chart(data, tax_levels_transl[tax_level])
		comparison_chart(data, tax_levels_transl[tax_level], normalized=True)
		data.apply(lambda x: generate_single_chart(x, tax_levels_transl[tax_level]), axis=1)


if __name__ == '__main__':
	run({'WORKING_PATH': 'tcc_selection/', 'CHARTS': {'phylum': True, 'class': True}})