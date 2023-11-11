import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy


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



data = []
keys = []

os.chdir('InputFiles')
for item in os.listdir(os.getcwd()):
	if os.path.isdir(item):
		if item == 'bold_selection':
			for file in os.listdir(item):
				if 'xlsx' in file:
					df = pd.read_excel(f'bold_selection/{file}', sheet_name='Espécies')
					df_counts = df[['phylum', 'Ocorrência']].groupby('phylum').sum().reset_index()
		

					# df_counts = df['phylum'].value_counts().reset_index(name='Ocorrência')
					data.append(df_counts)
					key = file.split('.')[0]
					keys.append(key)

			# data.extend([pd.read_excel(f'bold_selection/{file}', sheet_name='Espécies') for file in os.listdir(item) if 'xlsx' in file])
		else:
			dfs = [pd.read_excel(f'{item}/{file}', sheet_name='Espécies') for file in os.listdir(item) if 'xlsx' in file]
			df = get_degenerated_df(dfs)
			df_counts = df[['phylum', 'Ocorrência']].groupby('phylum').sum().reset_index()
			# df_counts = df['phylum'].value_counts().reset_index(name='Ocorrência')
			data.append(df_counts)
			keys.append(item.replace('degenerated_', ''))



df = pd.concat(data, keys=keys)

df.reset_index(level=0, inplace=True)
df.rename(columns={'level_0': 'Primer'}, inplace=True)


pivot = pd.pivot_table(data=df, index=['Primer'], columns=['phylum'], values='Ocorrência')


# Função para normalizar os valores em 100
def normalize_to_100(row):
	total = row.sum()  # Soma os valores na linha
	return (row / total) * 100


def comparison_chart(pivot, normalized=False):
	df = deepcopy(pivot)
	if normalized:
		df = df.apply(normalize_to_100, axis=1)

	all_primers_top5 = []
	# print(df)

	for index, row in df.iterrows():
		top_5 = list(row.sort_values(ascending=False).head(5).index)
		
		for i in top_5:
			if i not in all_primers_top5:
				all_primers_top5.append(i)

	others = []
	for col in df.columns:
		if col not in all_primers_top5:
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

	plt.title('Amplificações por Par de Primers: TOP 5 Filos')
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


def generate_single_chart(row):
	row = row.dropna()
	# top_5 = list(row.sort_values(ascending=False).head(5).index)
	# others = [idx for idx in row.index if idx not in top_5]

	# if len(others) > 0:
	# 	row.at['Outros'] = row.loc[others].sum()

	# row = row.drop(others)
	row = row.sort_values(ascending=False)
	

	ax = pd.DataFrame(row).T.plot.bar(stacked=True, color=['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000', '#FF0000', '#FFFF00', '#6AFF00', '#B9EDE0'], width=0.8)

	plt.xticks(rotation=0)
	ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
	plt.margins(0.2)
	plt.subplots_adjust(bottom=0.15, right=0.5, left=0.20)
	
	plt.title(f'Amplificações {row.name}: TOP 5 Filos', pad=15)
	plt.xlabel('Par de Primer', labelpad=6.0)
	plt.ylabel('Organismos', labelpad=6.0)

	autolabel(ax)

	plt.show()


comparison_chart(pivot)
comparison_chart(pivot, normalized=True)

# pivot.apply(generate_single_chart, axis=1)

