import json
import matplotlib.pyplot as plt

def read_json(inJSON):
	with open(inJSON,'rt') as cj:
		dataJSON = json.load(cj)
	return dataJSON

def box_plot(inArray,Labels,Title,Filename):
	fig,ax = plt.subplots()
	ax.boxplot(inArray,labels=Labels)
	ax.set_title(Title,fontsize='xx-large')
	#ax.tick_params(labelsize='xx-large')
	ax.set_xticklabels(Labels,rotation=45,ha='right',rotation_mode='anchor')
	plt.tight_layout()
	plt.savefig(Filename,format='svg')
	plt.close()

parser = argparse.ArgumentParser(description='Boxplot entrypoint script.')
parser.add_argument('input_dictionary', help='Input dictionary to plot.')
parser.add_argument('output_svgfile', help='Filename and path for output svg.')
args = parser.parse_args()

dict1 = read_json(args.input_dictionary)

labels, data = [*zip(*dict1.items())]
box_plot(data,labels,'Cosine Similarity',args.output_svgfile)