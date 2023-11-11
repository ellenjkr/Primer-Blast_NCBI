import xlsxwriter
import pandas as pd

df = pd.read_csv('2.csv', sep=';', encoding='latin-1')

class_counts = df['class'].value_counts().rename_axis('Class').reset_index(name='Qtd')

writer = pd.ExcelWriter(f'test.xlsx', engine='xlsxwriter')
class_counts.to_excel(writer, sheet_name='Teste', startrow=0, header=True, index=False)


workbook = writer.book
worksheet = writer.sheets['Teste']

chart = workbook.add_chart({'type': 'column', 'subtype': 'stacked'})

for _class in range(1, class_counts.shape[0] + 1):  # Percorre cada linha do dataframe
    chart.add_series({
        'name':       ['Teste', _class, 0],  # O nome está na coluna zero
        'values':     ['Teste', _class, 1, _class, 1],  # Os valores estão na linha percorrida, na coluna 1
        'gap':        2,
    })

    # chart.add_series({
    #     'name':       ['Teste', 0, col_num],
    #     'categories': ['Teste', 1, 0, 1, 0],
    #     'values':     ['Teste', 1, col_num, 1, col_num],
    #     'gap':        2,
    # })


chart.set_y_axis({'major_gridlines': {'visible': False}})

worksheet.insert_chart('D2', chart)

writer.save()

# https://pandas-xlsxwriter-charts.readthedocs.io/chart_stacked_column_farms.html