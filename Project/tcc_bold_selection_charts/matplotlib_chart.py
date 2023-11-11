import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

from io import BytesIO


dataframes = []
keys = []
for file in os.listdir('Primers_Plantas'):
    df = pd.read_excel(f'Primers_Plantas/{file}', sheet_name='Espécies')
    df_counts = df['genus'].value_counts().rename_axis('Espécie').reset_index(name='counts')
    dataframes.append(df_counts)
    key = file.split('.')[0]
    keys.append(key)

df = pd.concat(dataframes, keys=keys)

df.reset_index(level=0, inplace=True)
df.rename(columns={'level_0': 'Primer'}, inplace=True)


pivot = pd.pivot_table(data=df, index=['Primer'], columns=['Espécie'], values='counts')


ax = pivot.plot.bar(stacked=True, legend=False, figsize=(10, 6))
plt.xticks(rotation=0)

imgdata = BytesIO()
fig = ax.get_figure()
fig.savefig(imgdata, format="png")
imgdata.seek(0)
plt.show()


writer = pd.ExcelWriter(f'teste_matplotlib.xlsx', engine='xlsxwriter')
workbook = writer.book
pd.DataFrame().to_excel(writer, sheet_name='chart')
worksheet = writer.sheets['chart']
worksheet.insert_image(
    0, 0, "",
    {'image_data': imgdata}
)

writer.save()

# from io import BytesIO
# import matplotlib.pyplot as plt

# imgdata = BytesIO()
# fig, ax = plt.subplots()
# results.resid.hist(ax=ax)
# fig.savefig(imgdata, format="png")
# imgdata.seek(0)
