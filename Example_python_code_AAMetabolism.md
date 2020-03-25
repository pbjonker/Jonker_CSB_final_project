```python
##This is an example piece of code used to convert data obtained from link:
##http://amp.pharm.mssm.edu/Harmonizome/gene_set/Metabolism+of+amino+acids+and+derivatives/Reactome+Pathways
##This step uses regular expressions to isolate the genes of interest from a raw gene list. The raw list for this data
##set is in the Jonker_CSB_final_project repo, under the name "AAMetabolism.txt". This stands for Amino Acid Metabolism.
AAMetabolism = open("AAMetabolism.txt", 'r')
import re
AAMetabolism_raw = re.findall(r':.\w+.,', AAMetabolism.read())
##Since the output of the code is a list, i had to reconvert it back into strings, which i did below.
listToStr = ' '.join([str(elem) for elem in AAMetabolism_raw])
##Then, again, used regular expressions to reisolate the word-like objects (gene names)
AAMetabolism = re.findall(r'\w+',listToStr)
##Export to a CSV
import csv
with open('AAMetabolism.csv', 'w') as AAMetabolism1:
    wr = csv.writer(AAMetabolism1, quoting=csv.QUOTE_ALL)
    wr.writerow(AAMetabolism)
```


```python

```
