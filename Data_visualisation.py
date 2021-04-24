import pandas as pd
import matplotlib.pyplot as plt

features = pd.read_csv('features.csv')

print(features)

surf = features[['protIDs', 'surfaces']]
surf.sort_values(by=['surfaces'], inplace=True)

print(surf)

#plt.plot(surf, bins=100)  # density=False would make counts
plt.rcParams.update({'font.size': 5})
surf.plot(kind='bar', x='protIDs', y='surfaces')
plt.title('Surface area of proteins in squared Angstrom')
plt.xlabel('Protein ID')
plt.savefig('surface.png', bbox_inches='tight')





