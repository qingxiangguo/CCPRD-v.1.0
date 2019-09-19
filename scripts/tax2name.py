#!/opt/python2.7/bin/python2.7

f =open('2','r')
l_line = f.readlines()
f.close()

from ete3 import NCBITaxa
ncbi = NCBITaxa()

for line in l_line:
	taxid2name = ncbi.get_taxid_translator([line])
	print taxid2name,
