

from scipy import stats

import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt

import urllib


def get_gene_name(entry_tag,attribute):
	# <name type="ordered locus">YLR420W</name>
	gene_names=entry_tag.find('.//{http://uniprot.org/uniprot}gene').getchildren()
	name_attribute='|'.join([name.text for name in gene_names if name.attrib['type']==attribute])
	return name_attribute

def get_name(entry_tag):
	return entry_tag.find('.//{http://uniprot.org/uniprot}fullName').text


def get_entry_tag(entry):
	for entry_tag in entries:
		if entry_tag.getchildren()[0].text==entry:
			return entry_tag



url='https://www.uniprot.org/uniprot/?query=proteome:UP000002311&format=xml&cols=id'

page=urllib.urlopen(url).read()

filename_xml='uniprot_proteome_sce.xml'
open(filename_xml,'w').write(page)


tree = ET.parse(filename_xml)
root = tree.getroot()
entries=root.findall('.//{http://uniprot.org/uniprot}entry')



proteins=pd.read_excel('2016 Rabin Supp Data.xlsx',sheet_name='Protein log2 RA').fillna('')
dilutions=['0.05','0.11','0.16','0.22','0.30']
# tuple of limiting substrate and marker color
substrates=[('C','k'),('N','b'),('P','g'),('U','r'),('L','y')]

aybrah=pd.read_excel('https://github.com/kcorreia/aybrah/raw/master/aybrah.xlsx',sheet_name='aybrah').set_index('FOG').fillna('')
aybraham=pd.read_excel('https://github.com/kcorreia/aybraham/raw/master/aybraham.xlsx',sheet_name='reactions').fillna('').set_index('id')

# cytosolic ribosome FOGs
fogs=aybraham['FOG']['BIOMASS_PROTEIN_iMM904_TRNA_1G']

# clean up fog expression
for char in list('()?& '):
	fogs=fogs.replace(char,'|')


for ax in axs.reshape(-1): 
  ax.set_ylabel(str(i))

completed_fogs=[]

for fog in sorted(filter(None,set(fogs.split('|')))):
#for fog in ['FOG00007']:
	print fog
	if os.path.isfile(fog+'.png'):
		continue
	seqids=aybrah['sce'][fog].split(';')
	# check for paralogs/ohnologs
	children=aybrah[(aybrah.Parent.str.contains(fog))]
	if len(children)>0:
		fogs_added=children['sce'].index.tolist()
		seqids.extend(    ';'.join(children['sce'].tolist()).split(';')   )
	seqids=filter(None,seqids)
	if len(seqids)==0:
		continue
	fig, axes= plt.subplots(1, len(seqids))
	if len(seqids)>1:
		axes_iterate=axes.reshape(-1)
	else:
		axes_iterate=[axes]
	for seqid,ax in zip(seqids,axes_iterate):
		entry=seqid.split('|')[1]
		print('get gene information')
		entry_tag=get_entry_tag(entry)
		locus_ordered=get_gene_name(entry_tag,'ordered locus')
		name_gene=get_gene_name(entry_tag,'primary')
		name_full=get_full_name(entry_tag)
		# find relevant row in Rabin supplemental data
		df_temp=proteins[proteins.Gene.str.contains(locus_ordered)]
		# not all genes have data
		ax.set_title(name_gene)
		if len(df_temp)==0:
			print('no expression data')
			continue
		else:
			index=df_temp.index.item()
		substrates=[('C','k'),('N','b'),('P','g'),('U','r'),('L','y')]
		filename=name_gene+' '+name_full
		print('create plot')
		for substrate,color in substrates:
			dilutions_float=[float(D)for D in dilutions]
			expression=[proteins[substrate+D][index] for D in dilutions]
			slope, intercept, r_value, p_value, std_err = stats.linregress(dilutions_float,expression)
			expression_predicted=[slope*D+intercept for D in dilutions_float]
			ax.plot( dilutions_float, expression, color+'s',label=substrate+'-lim: r'+r'$^{2}$'+str(round(r_value**2,2))+', slope '+str(round(slope,2)))
			ax.plot( dilutions_float, expression_predicted, color+':')
		print('edit axis plot etc')
		ax.set_xlabel('Dilution rate '+r'$h^{-1}$')
		ax.set_ylabel('Protein log2 RA')
		ax.set_ylim(top=3.0)
		ax.legend(loc='lower right', prop={'size': 8})
	plt.suptitle(fog)
	#plt.title(name_gene+' '+locus_ordered+': '+name_full)
	plt.savefig(fog+'.png')
	plt.close()



