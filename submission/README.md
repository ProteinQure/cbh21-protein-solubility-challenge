#Transcendent Flavour   
This is our submission for the CBH2021 hackathon.
###Welcome to Flavour Town 
[Github](https://github.com/MatthewParkinGit/cbh21-protein-solubility-challenge)
### The Team  
#### K. Armstrong 
- PhD student in Computer Science, Bioinformatics(M.Sc.), Biochemistry(BSc) - Machine learning, deep learning and Bioinformatics

 ##### G. Grey 
- PhD student in Biochemistry, Biochemistry (M.Bio.), Biochemistry(BSc) - structural and functional protein biochemistry

##### M. Parkin 
- Post-Apprenticeship - General Computing, Software Engineering and Networking

####Aims  
1 - Learn from both sequence and structure  
2 - Ensure features are representative of all the data  
3 - Machine Learning -> predictions.py  

###Machine Learning Outline
- Preprocessor -> ColumnTransformer with KBins discretiser and onehot encoding
- Selector -> Recursive Feature Elimination
- Regressor -> Random Forest Regressor
- Parameter Tuning -> GridSearchCV with a verbose of 10 and cv of 5

###Statistical potentials/Features used in the data frame  
Counts for each type of amino acid in each protein  
Counts for the percentage of the amino acids in each protein  
Length of protein  
Gravy  
Aromaticity  
Turns  
Helices  
Sheets  
Solvent Accessible Surface  
Polar and Apolar Surface Area   
Iso Electric Point 
Instability Index  
Charge at PH (1,7,14)

<t> [659 rows x 55 columns] </t>


###Improvements  
More data points
Different machine learning model
    - Try different preprocessing
    - Try deep learning approach