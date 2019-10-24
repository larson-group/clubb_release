#!/home/roehls/python_plots/my_py/bin/python

###############################################################################
# analyze_results.py
#
# This python script will take the tuned results and create a scatter plot,
# a box plot, and a Plotgen plot making it visually easy to compare the tuned
# results.  The scatter plot, plots all of the results, while the box plot only
# plots the top 25% and Plotgen only plots the top 10 of each case.
#
# File History:
#  v1.0: Initial Release
#
# Notes:
#  To use this program you MUST have pylab and numpy installed.  To see instructions on 
#  how to do this see the TWiki page: 
#  https://github.com/larson-group/sys_admin/wiki/Pylab
#  You must then also change the location of python to the virtual environment 
#  python location.
# Required Arguments
#  -The folders that contain the ens_tune_* folders.
###############################################################################

import os
import glob
import re
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages

# Set up our colors that we will use, they are: orange, blue, green
# purple, black
ourColors =[(1,.65,0),(0,0,1),(0,1,0),(.63,.13,.94),(0,0,0)];


# Set up some global variables
name = []
values = []
costFunc = []
numCases = len(sys.argv)-1

# Set up all of the arrays that we will need
files = []
folders = []
allCost = zeros((numCases,2000))
allNames = zeros((numCases,2000))
allParam = zeros((numCases,4000,53))
topResults = ones((numCases,500))
topResultsName = ones((numCases,500))


##############################################################
# readParameterFile
# Function that will read the parameter file
# Fills the name array with parameters names and
# Fills the values array with parameters values
# Arguments: String, file to be read
##############################################################
def readParameterFile(file):
	input = []
	list([values.pop() for z in xrange(len(values))])	#Clear the array
	currFile = open(file)

	for line in currFile:
		input = re.split('=',line)
		if len(input) == 2:
			name.append(input[0].replace(' ','')) 
			values.append(float(input[1]))
	currFile.close()


##############################################################
# readCostFunction
# Function that will read the cost function of a given ensemble run
# It must open the tune.log file and get the cost function
# Arguments: String, folder where tune.log is located
##############################################################
def readCostFunction(folder):
	input = []
	currFile = open(folder+'/tune.log')
	for line in currFile:
		input = re.split('\$\$',line)
		if len(input) > 1:
			costFunc.append(float(input[1].replace('\n','').replace(' ','')))
			break


##############################################################
# writeBoxHtml
# A function that will write the html file the box plots.
# It looks through all the box plots created and makes a html
# page for them so they can be easily viewed
# Arguments: None
##############################################################
def writeBoxHtml():
	boxFile = open('plots/box/box.html', 'w')
	boxFile.write("<html>\n<head>\n<title>Box Plots</title>\n</head>\n<body>\n")
	boxFile.write("<table border=\"1\" ALIGN=\"center\"><tr><th>Legend</th></tr>")
	for i in range(1,len(sys.argv)):
		boxFile.write("<tr><td>" + str(i) +"</td><td>"+ sys.argv[i]+"</td></tr>")
	boxFile.write("</table>\n")
	for pngFiles in glob.glob("plots/box/*.png"):
		pngFiles = re.split('/',pngFiles)
		pngFiles = pngFiles[2]
		boxFile.write("<img width=\"324\" height=\"312\" align=\"BOTTOM\" border=\"0\" style=\"padding: 5px;\" src=\"" + pngFiles + "\" alt=\"" + pngFiles + "\" />\n")
	boxFile.write("</div>\n</body>\n</html>")
	boxFile.close()



##############################################################
# writeScatterHtml
# A function that will write the html file the scatter plots.
# It looks through all the scatter plots created and makes a html
# page for them so they can be easily viewed
# Arguments: None
##############################################################
def writeScatterHtml():
	scatterFile = open('plots/scatter/scatter.html', 'w')
	scatterFile.write("<html>\n<head>\n<title>Scatter Plots</title>\n</head>\n<body>\n")
	for pngFiles in glob.glob("plots/scatter/*.png"):
		pngFiles = re.split('/',pngFiles)
		pngFiles = pngFiles[2]
		scatterFile.write("<img width=\"324\" height=\"312\" align=\"BOTTOM\" border=\"0\" style=\"padding: 5px;\" src=\"" + pngFiles + "\" alt=\"" + pngFiles + "\" />\n")
	scatterFile.write("</div>\n</body>\n</html>")
	scatterFile.close()


##############################################################
# writeIndexHtml
# A function that will write the index html file.  This simply
# creates a page with links to correct page
# Arguments: None
##############################################################
def writeIndexHtml():
	indexFile = open('plots/index.html', 'w')
	indexFile.write("<html>\n<head>\n<title>index</title>\n</head>\n<body>\n")
	indexFile.write("<a href=\"scatter/scatter.html\">Scatter Plots</a><br>\n")
	indexFile.write("<a href=\"box/box.html\">Box Plots</a><br>\n")
	indexFile.write("<a href=\"plotgen/index.html\">Plotgen Plots</a><br>\n")
	indexFile.write("</body>\n</html>")
	indexFile.close()


##############################################################
# fileReading
# A function that will write the index html file.  This simply
# creates a page with links to correct page
# Arguments: None
##############################################################
def fileReading():
	index = 0
	for runs in range(0,numCases):
		list([costFunc.pop() for z in xrange(len(costFunc))])		#Clear the array
		tempFolder=glob.glob(sys.argv[runs+1]+'/ens_tune_*')		#Get all the tuning folders for each case
		folders.append([])
		for p in range(0,len(tempFolder)-1):
			folders[runs].append(tempFolder[p])			#append them to the folder array
		folders[runs].sort()	
		files.append([])
		for fold in folders[runs]:					#For each folder find:
			fileNum=re.sub('.+/.+_','',fold)			
			allNames[runs][index]=fileNum				#the name
			files[runs].append(glob.glob(fold+'/tunable_parameters*.in')[0]) #the locatin of tunable_parameters.in
			readCostFunction(fold)					#the cost function
			for p in range(0,(len(costFunc))):			#Append cost function to all cost functions
				allCost[runs][p]=costFunc[p]
			index = index + 1
		temp = list(topResults[runs])					#Set up temp arrays
		tempName = list(topResultsName[runs])
		for numFile in range(0,len(costFunc)):				#Search for the lowest cost functions
			for numResult in range(0,(len(costFunc))):
				if(allCost[runs][numFile] < temp[numResult]):
					temp.insert(numResult,allCost[runs][numFile])
					tempName.insert(numResult,allNames[runs][numFile])
					break
		for numResult in range(0,len(costFunc)):			#Set the global arrays equal to the temp arrays
			topResults[runs][numResult] = temp[numResult]
			topResultsName[runs][numResult] = tempName[numResult]
		index = 0 

	for k in range(0,numCases):						#Reads all of the files parameters
		i = 0
		for paraFile in files[k]:
			readParameterFile(paraFile)
			for p in range(0,(len(values)-1)):
				allParam[k][i][p] = values[p]
			i = i +1
			if(values[p] == 0):
				print "Error at i=",i," and k=",k
				exit()





##############################################################
# plotScatterPlots
# Creates the scatter plots for all of the points
# It saves each scatter plot as a .png
##############################################################	
def plotScatterPlots():
	for j in range(0,len(values)):
	#for j in range(0,1):
		if (allParam[0][0][j] == allParam[0][1][j] == allParam[0][2][j] == allParam[0][3][j]) :	
			print "Skipping Plot ", name[j]
		else:
			#Creates the scatter plots
			maxX=allParam[0][0][j]
			minX=allParam[0][0][j]
			maxY=allCost[0][0]
			print "Creating Plot ", name[j]
			for k in range(0,numCases):
				for i in range(0,len(files[k])):
					scatter(allParam[k][i][j],allCost[k][i],s=0,color='white')
					plt.text(allParam[k][i][j],allCost[k][i],int(allNames[k][i]),size='xx-small',color=ourColors[k%5])
					if(allParam[k][i][j] > maxX):
						maxX = allParam[k][i][j]
					if(allParam[k][i][j] < minX):
						minX = allParam[k][i][j]
					if(allCost[k][i]>maxY):
						maxY = allCost[k][i]
			for t in range(0,numCases):
				plt.text(maxX,(maxY-((maxY/100)*(t*2))),sys.argv[t+1].rstrip('/'),color=ourColors[t%5],size='medium')
			title(name[j]+ ' Parameter vs Cost Function')
			xlabel(name[j]+ ' Value')
			ylabel('Cost Function')	
			grid(True)
			draw()
			savefig("plots/scatter/"+str(j)+"temp.png")
			#pdfScatter.savefig()
			clf()


##############################################################
# plotBoxPlots
# Creates the box plots for the top 25% of the points
# It saves each box plot as a .png
##############################################################	
def plotBoxPlots():
	for j in range(0,len(values)):
	#for j in range(0,1):
		if (not(allParam[0][1][j] == allParam[0][1][j] == allParam[0][2][j] == allParam[0][3][j])):

			limit = zeros(len(sys.argv))
			data = []
			highestLimit = 0
			legendStr = ""
			for k in range(0,numCases):
					limit[k] = len(files[k])/4
					if(limit[k] > highestLimit):
						highestLimit = limit[k]
			top25Results = zeros(((numCases),highestLimit))
			top25ResultsName = zeros(((numCases),highestLimit))
			for k in range(0,numCases):
				for i in range(0,int(limit[k])):
					top25Results[k][i]=allParam[k][i][j]
					top25ResultsName[k][i]=topResultsName[k][i]
			box = []
			allBox = []
			#Creates the Box and Wisker Plots
			for k in range(0,numCases):
				data.append(top25Results[k])
			boxplot(data)
			ylabel('Value')
			title(name[j])
			grid(True)
			savefig("plots/box/"+str(j)+"temp.png")
			clf()
				
				
##############################################################
# plotPlotgen
# Runs plotgen for all of the cases. 
##############################################################	
def plotPlotgen():
	plotgenString = "plotgen -c -l -r -ensemble "
	for k in range(0,numCases):
		for i in range(0,10):
			plotgenString += str(sys.argv[k+1]) + '/ens_tune_' + str(int(topResultsName[k][i])) + '/ '
	plotgenString += "./plots/plotgen"
	os.system(plotgenString)


# Clear out the folders and make sure they exsist
os.system("rm -rf plots && mkdir plots")
os.system("mkdir plots/scatter")
os.system("mkdir plots/box")


fileReading()
plotScatterPlots()
plotBoxPlots()
plotPlotgen()

writeBoxHtml()
writeScatterHtml()
writeIndexHtml()

os.system("rm -rf tuner_results.maff && zip -r tuner_results.maff plots/")
