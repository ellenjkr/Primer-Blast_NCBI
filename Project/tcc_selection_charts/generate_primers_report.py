import pandas as pd
import matplotlib.pyplot as plt
import os

from io import BytesIO



writer = pd.ExcelWriter(f'reports/tcc_bold_selection_report.xlsx', engine='xlsxwriter')
workbook = writer.book


classifications_keys = ['species', 'genus', 'family', '_order', 'class', 'phylum', 'kingdom', 'superkingdom']
classifications = {key: {'data': [], 'excel_chart': True} for key in classifications_keys}
species_data = []

os.chdir('01_xlsx_files')
for key in classifications.keys():
	for file in os.listdir(os.getcwd()):
		df = pd.read_excel(file, sheet_name='Espécies')
		counts = df[key].value_counts()
		classifications[key]['data'].append(counts)
		if df[key].unique().shape[0] > 255:
			classifications[key]['excel_chart'] = False


index = [f"{file.replace('.xlsx','')}" for file in sorted(os.listdir(os.getcwd()))]

for key in classifications.keys():
	if classifications[key]['excel_chart'] is True:
		chart = workbook.add_chart({'type': 'column', 'subtype': 'stacked'})
		df = pd.DataFrame(classifications[key]['data'], index=index)

		sheet_name = key.capitalize()
		df.to_excel(writer, sheet_name=f'{sheet_name} Report', startrow=0, header=True, index=True)

		# Configure the series of the chart from the species_dataframe data.
		for col_num in range(1, df.shape[1] + 1):
			chart.add_series({
				'name': [f'{sheet_name} Report', 0, col_num],
				'categories': [f'{sheet_name} Report', 1, 0, df.shape[0], 0],
				'values': [f'{sheet_name} Report', 1, col_num, df.shape[0], col_num],
				'gap': 100,
			})

		chart.style = 2
		chart.set_y_axis({'major_gridlines': {'visible': False}})
		chart.set_legend(
			{'position': 'bottom'}
		)
		# pd.DataFrame().to_excel(writer, sheet_name=f'{sheet_name} Chart', startrow=0, header=False, index=False)
		worksheet = writer.sheets[f'{sheet_name} Report']
		worksheet.insert_chart('B7', chart)

	else:
		df = pd.DataFrame(classifications[key]['data'], index=index)
		sheet_name = key.capitalize()
		# df.to_excel(writer, sheet_name=f'{sheet_name} Report', startrow=0, header=True, index=True)

		classifications[key]['data'] = [classification_df.rename_axis('Espécie').reset_index(name='counts') for classification_df in classifications[key]['data']]

		concat_df = pd.concat(classifications[key]['data'], keys=index)

		concat_df.reset_index(level=0, inplace=True)
		concat_df.rename(columns={'level_0': 'Primer'}, inplace=True)

		pivot = pd.pivot_table(data=concat_df, index=['Primer'], columns=['Espécie'], values='counts')

		ax = pivot.plot.bar(stacked=True, legend=False, figsize=(10, 6))
		plt.xticks(rotation=0)

		imgdata = BytesIO()
		fig = ax.get_figure()
		fig.savefig(imgdata, format="png")
		imgdata.seek(0)

		# pd.DataFrame().to_excel(writer, sheet_name=f'{sheet_name} Chart', startrow=0, header=False, index=False)
		worksheet = writer.sheets[f'{sheet_name} Report']
		worksheet.insert_image(
			6, 1, "",
			{'image_data': imgdata}
		)

writer.close()

# https://holypython.com/python-visualization-tutorial/colors-with-python/
# https://medium.com/@jb.ranchana/easy-way-to-create-stacked-bar-graphs-from-dataframe-19cc97c86fe3