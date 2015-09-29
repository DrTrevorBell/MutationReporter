#!/usr/bin/python

# 17 June 2010

import datetime

starttime = datetime.datetime.now()

import climb, cgi, os, re, sys, socket

# import cgitb
# cgitb.enable()

# sys.stderr = sys.stdout

form = cgi.FieldStorage()

version = "<h1>Mutation Reporter Tool</h1>"

def closePage():
	print headerLinks
	print '<hr><input type="button" value="Go Back" onclick="goBack()" /><br>'
	# <br><a href="/mrt/index.html">Submit another</a>'
	# print '<a href="/virologyjournalreview/index.html">Submit another</a>'
	print '<hr>'
	print '<tt><p align=right>Version 1.3 (July 2012)<br>'
	endtime = datetime.datetime.now()
	print 'Request from %s<br>' % (cgi.os.environ['REMOTE_ADDR'])
	print 'Served by %s at %s<br>' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0])
	print 'Run completed %s<br>' % (endtime.strftime("%Y-%m-%d %H:%M"))
 	print '%03.6f seconds</tt></p>' % ((endtime - starttime).microseconds / float(1000000))
	print '</body></html>'

	f = open('/var/www/tmp/fm.log', 'a')
	f.write( 'Mutation Reporter\t' )
	f.write( 'Request from %s\t' % (cgi.os.environ['REMOTE_ADDR']) )
	f.write( 'Served by %s at %s\t' % (socket.gethostname(), socket.getaddrinfo(socket.gethostname(), None)[0][4][0]) )
	f.write( 'Run completed %s\t' % (endtime.strftime("%Y-%m-%d %H:%M")) )
	f.write( '%03.6f seconds\t' % ((endtime - starttime).microseconds / float(1000000)) )
	f.write( '\n')
	f.close()

def badFile(errorText, explainText):
	print version
	print "<h2>%s</h2>" % (errorText)
	print explainText
	closePage()
	# file.close()
	sys.exit()

print 'Content-Type: text/html'
print

print '<!DOCTYPE html'
print '	PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"'
print '	 "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">'
print '<html xmlns="http://www.w3.org/1999/xhtml" lang="en-US" xml:lang="en-US">'
print '<head><title>HVDR Mutation Reporter Tool</title>'
print '<link rel="stylesheet" type="text/css" href="/res/hvdr.css"/>'
print '<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"/>'
print '</head>'
print '<script type="text/javascript">'
print 'function goBack()'
print '  {'
print '  window.history.back()'
print '  }'
print '</script>'
print '<body>'

# email = form["emailaddress"].value
# if len(email) <= 6 or email.find('@') == -1 or email.find('.') == -1:
# 	badFile("Bad email address!", "")

headerLinks = ''	# required when error triggered before processing file

infile = form["infile"].file
filename = form["infile"].filename
if len(filename) == 0:
	badFile('Error: Specify the input file!', '')

loci = form['loci'].value
if len(loci) == 0:
	badFile('Error: Bad loci position!', '')
else:
	loci = loci.replace(' ', '')

anchorMotif = form['anchormotif'].value
anchorMotif = anchorMotif.strip()	# remove leading and trailing whitespace
anchorMotif = anchorMotif.upper()	# prevents searching for lowercase motifs but better from a user-error point of view
if len(anchorMotif) == 0:
	badFile('Error: Bad anchor motif!', '')

anchorPosition = form['anchorposition'].value
if len(anchorPosition) == 0:
	badFile('Error: Bad anchor position!', '')

truncation = form['truncation'].value

aa = form['sequencetype'].value == 'aa'

if form.has_key('suppresszeros'):
	sZ = True
else:
	sZ = False

if form.has_key('percentages'):
	perc = True
else:
	perc = False

if form.has_key('aminoacidcoloroutput'):
	highlightAA = True
else:
	highlightAA = False

if form.has_key('sortsequences'):
	sortSeq = True
else:
	sortSeq = False

if form.has_key('outputgrouping'):
	grouping = form["outputgrouping"].value
	grouping = grouping.replace(' ', '')
else:
	grouping = ''

referenceMotif = form['referencemotif'].value
graphWidth = form['graphwidth'].value
graphHeight = form['graphheight'].value
graphTitle = form['graphtitle'].value
if form.has_key('suppressconserved'):
	suppressConserved = True
else:
	suppressConserved = False
if form.has_key('stacked'):
	positionMode = 'stack'
else:
	positionMode = 'dodge'
if form.has_key('graphscale'):
	graphScale = '+ ylim(0,100)'
else:
	graphScale = ''

if form.has_key('discoverymode'):
	discoveryMode = True
	grouping = ''		# just to be sure
	referenceMotif = ''	# just to be sure
else:
	discoveryMode = False

# os.system('beep -f 8192 -r 16 -l 32 -d 64')	# make a noise

randomToken = climb.randomStamp()
climb.TEMPFOLDER = '/var/www/tmp/'
climb.TEMPFOLDER += 'MRT-' + randomToken + '/'

os.makedirs(climb.TEMPFOLDER)

tempF = climb.TEMPFOLDER + filename
temp = open(tempF, 'wb')
for line in infile:
	temp.write(line)
temp.close()
infile.close()

S = climb.Sequence()
S.load(tempF, filterPattern=form["pattern"].value, filterInclude=(form["filter"].value == "include"))

S.seqCase()

if sortSeq:
	S.seq.sort()

if len(truncation) != 0:
	S.truncate = int(truncation)

# motifPos = S.find(anchorMotif)[0][2]	# extract the starting nucleotide position from the first result
# motifPos = re.search(anchorMotif, S.seq[0]['seq']).start() + 1	# S.find adds one

motifPos = re.search(anchorMotif, S.seq[0]['seq'])
if motifPos != None:
	motifPos = motifPos.start() + 1
	numMotifs = re.findall(anchorMotif, S.seq[0]['seq'])
else:
	badFile('Anchor motif (%s) not found' % anchorMotif, '')

cropRange = str(motifPos) + '-' + str(S.seqLength()[0][1])	# length of first sequence
S.nucCopy(cropRange)

# Error checking to ensure that motif matches with grouping

if aa:
	setOne = climb.AMINOACIDS[:]
	setTwo = [climb.GAP, '?', 'X', '*']
else:
	setOne = climb.BASES[:]
	setTwo = climb.NONBASES[:]

if discoveryMode:
	tempLoci = []
	for i in climb.parseList(loci):
		for j in range(i[0], i[1]+1):
			tempLoci.append(str(j))		# convert loci to a list of strings for use later

	tempColumn = 0
	nonConservedLoci = ''
	discoveryDistribution = S.baseDistribution(loci, set1=setOne, set2=setTwo, distributionMapping=int(anchorPosition))
	for i in discoveryDistribution:
		discoveryCount = 0
		for j in i:
			if i[j] > 0:
				discoveryCount += 1
		if discoveryCount > 1:			# more than one residue found in column
			nonConservedLoci = nonConservedLoci + tempLoci[tempColumn] + ','
		tempColumn += 1
	loci = nonConservedLoci[:-1]	# remove final comma
	grouping = ''			# disable grouping
	if len(loci) == 0:		# all positions are conserved
		badFile('Discovery Mode: All positions are conserved', '')

# if no grouping is specified, use grouping = length of motif
if grouping == '':
	doGrouping = False
	cc = 0
	for i in climb.parseList(loci):
		cc += i[1] - i[0] + 1
	grouping = str(cc)
else:
	doGrouping = True

groupList = grouping.split(',')

# create a list of each actual output position
colList = []
for i in climb.parseList(loci):
	for j in range(i[0], i[1]+1):
		colList.append(str(j))

temp = 0
for i in groupList:
	temp += int(i)

if temp != len(colList):
	badFile("Grouping does not match loci!", "")

if len(referenceMotif) > 0:
	if len(referenceMotif) != temp:
		badFile("Reference motif length is incorrect!", "")

print '<h1>Results</h1>'
print '<br>'

if referenceMotif != '':
	graphLink = '&bull; <a href="#graph" style="color:green;">Mutation Distribution Graph</a> '
else:
	graphLink = ''
# button to jump to a target anchor; instead use a link for "Back"
# headerLinks = '<center><input type="button" value="Go Back" onClick="javascript:location.href = \'#graph\';" /> &bull; <a href="#inputparameters">Input Parameters</a> &bull; <a href="#locidistribution">Loci Distribution</a> &bull; <a href="#residuedistributionsummary">Residue Distribution Summary</a> %s&bull; <a href="#motifdistribution">Motif Distribution</a> &bull;</center><br>' % (graphLink)

headerLinks = '<tt><center>&bull; <a href="javascript:javascript:history.go(-1)" style="color:green;">&larr; Back</a> &bull; <a href="#inputparameters" style="color:green;">Input Parameters</a> &bull; <a href="#locidistribution" style="color:green;">Loci Distribution</a> &bull; <a href="#residuedistributionsummary" style="color:green;">Residue Distribution Summary</a> %s&bull; <a href="#motifdistribution" style="color:green;">Motif Distribution/Graph</a> &bull;</p></center></tt><br>' % (graphLink)

######################
## Input Parameters ##
######################

print '<a name="inputparameters"><h2>Input Parameters</h2></a><br>'
print headerLinks
print '<table border="1" style="font-family:monoscape; font-size:16px">'

def printTableRow(p1, p2):
	print '<tr><td>%s</td><td><b>%s</b></td></tr>' % (p1, p2)

printTableRow('Input file', filename)
printTableRow('Number of Sequences', str(len(S.seq)))
printTableRow('Anchor Motif', anchorMotif)
printTableRow('Anchor Position', str(anchorPosition))
printTableRow('Loci', loci)
if doGrouping:
	printTableRow('Output grouping', grouping)
printTableRow('Anchor Motif Found at Position', str(motifPos))
outNote = 'Anchor Motif Occurrences'
if len(numMotifs) > 1:
	outNote += '<br>(Only the first occurrence is used)'
printTableRow(outNote, str(len(numMotifs)))
if len(referenceMotif) > 0:
	printTableRow('Reference motif', referenceMotif)
if discoveryMode:
	printTableRow('Discovery Mode', 'On')
print '</table><br>'


####################################################################################################
### Table 1: Sequence IDs, loci in columns are per "output grouping" and residues in each column ###
####################################################################################################
print '<a name="locidistribution"><h2>Loci Distribution</h2></a><br>'
print headerLinks
print '<table border="1" style="font-family:monospace; font-size:16px">'

# output the correct number of positions for each column
print '<tr><td><b>ID</b></td>'
previous = 0
for j in groupList:
	print '<td><b>'
	for k in colList[0:int(j)]:
		print k, '<br>'	# extra <br> after last item or 1 <br> after only one item, but seems not to affect HTML output
	print '</b></td>'
	for k in range(int(j)):
		del colList[0]		# remove 'k' entries from colList
print '</tr>'

for i in S.extract(loci, mapping=int(anchorPosition)):
	print '<tr><td>%s</td>' % (i[0])

	start = 0
	end = 0
	for j in groupList:
		end += int(j)
		print '<td align="center">%s</td>' % (i[1][start:end])		# center the content of the cell (one base in a column with 4 digit header for example)
		start = end
	print '</tr>'

	# print '<tr><td>%s</td><td>%s</td></tr>' % (i[0], i[1])
print '</table>'
print '<br>'

##########################################################################
### Table 2: Summary table -- all residues as rows and loci as columns ###
##########################################################################
print '<a name="residuedistributionsummary"><h2>Residue Distribution Summary</h2></a><br>'
print headerLinks
print '<table border="1" style="font-family:monospace; font-size:16px">'
lociDistribution = S.baseDistribution(loci, set1=setOne, set2=setTwo, distributionMapping=int(anchorPosition))
print climb.distributionShow(lociDistribution, loci, len(S.seq), rowList=setOne+setTwo, html=True, suppressZeros=sZ, percentage=perc, firstCell="&nbsp;")   # or <b>Residue</b>
print '</table>'
print '<br>'

#######################
### Plot mutationts ###
#######################
# lociDistribution is a list of dictionaries; each list entry represents one locus; the keys are the residues and the values are the raw counts

if len(referenceMotif) > 0:
	lociList = []
	for i in climb.parseList(loci):
		for j in range(i[0], i[1]+1):
			lociList.append(j)

	referencePosition = 0
	series = [['Locus', 'Residue', 'Percentage', 'Order']]		# Order used later to enforce factor/level order for plot when locus preceeded by WT residue
	for i in lociDistribution:
		for j in i:	# j is the key
			# construct Name and Value lists of the residues and counts to include in the graph by including
			# only those with counts which are not the reference residue and counts which are greater than zero (option?)
			# Also include the locus
			if (j != referenceMotif[referencePosition]) and (i[j] != 0):
				series.append([referenceMotif[referencePosition]+str(lociList[referencePosition]), j, i[j]/float(len(S.seq)) * 100, lociList[referencePosition]])
			elif j == referenceMotif[referencePosition] and (i[j]/float(len(S.seq))) == 1 and not suppressConserved:
				series.append([referenceMotif[referencePosition]+str(lociList[referencePosition]), j, 0, lociList[referencePosition]])
		referencePosition += 1

	# Write data out as a csv file and then import this into R
	# Easier than trying to send the data over as a data frame

	import csv
	graphFile = open(climb.TEMPFOLDER + 'graphFile.csv', 'wb')
	graphOut = csv.writer(graphFile)
	graphOut.writerows(series)
	graphFile.close()

	# creating the graph via rpy resulted in no output
	# this may be related to passing in the variable 'x'
	# doing this via a script file instead
	# import rpy

	# sanity checks on the graph parameters
	if not graphWidth.isdigit():
		graphWidth = '480';
	if not graphHeight.isdigit():
		graphHeight = '480';
	graphTitle = graphTitle.replace('\\', ' ')

	fontSize1 = int(int(graphWidth) / 54)
	fontSize2 = int(int(graphWidth) / 60)

	graphFile = open(climb.TEMPFOLDER + 'graphFile.R', 'w')
	graphFile.write("x1 = read.csv('graphFile.csv')\n")
	graphFile.write("library('ggplot2')\n")
	graphFile.write("png('graphFile.png', width=%s, height=%s, units='px')\n" % (graphWidth, graphHeight))
	# http://stackoverflow.com/a/11364250/580010
	graphFile.write("residues = c('C', 'T', 'G', 'A', '-', '?', 'X', '*', 'N', 'W', 'S', 'M', 'K', 'R', 'Y', 'B', 'D', 'H', 'V', 'E' ,'F', 'I', 'J', 'L', 'O', 'P' ,'Q', 'U', 'Z')\n")
	graphFile.write("color_key = c('blue','red','black','green','gray', 'lightgray', 'dimgray', 'maroon', 'darkgrey', 'yellow', 'darkblue', 'cyan', 'darkred', 'darkgreen', 'magenta', 'orange', 'darkgoldenrod4','darkolivegreen1', 'palevioletred1', 'lightblue', 'plum', 'lightgreen', 'limegreen', 'pink', 'moccasin', 'sandybrown', 'navajowhite3', 'mistyrose1', 'rosybrown3')\n")
	graphFile.write("names(color_key) = unique(c(as.character(residues)))\n")
        graphFile.write("x1 = transform(x1, Locus=reorder(Locus, Order))\n")	# re-order by Locus (http://stackoverflow.com/a/3744432/580010)
	# http://stackoverflow.com/a/1331400/580010
	graphFile.write("p = ggplot(x1,aes(x = factor(Locus))) %s + geom_bar(aes(fill = Residue, weight = Percentage), position='%s') + labs(x = 'Locus',y = 'Percentage') + opts(title = '%s', axis.text.x = theme_text(size = %s, angle=90, hjust=1), axis.title.x = theme_text(size = %s), axis.text.y = theme_text(size = %s), axis.title.y = theme_text(size = %s), legend.text = theme_text(size = %s), legend.title = theme_text(size = %s), panel.grid.major = theme_line(size = 0.5, colour = 'darkgray'), panel.grid.minor = theme_blank(), panel.background = theme_blank(), axis.ticks = theme_blank())\n" % (graphScale, positionMode, graphTitle, fontSize1, fontSize2, fontSize1, fontSize2, fontSize2, fontSize2))
	graphFile.write("p + scale_fill_manual(values = color_key[names(color_key) %in% x1$Re])\n")
	graphFile.write("dev.off()\n")
	graphFile.close()

	import subprocess
	# http://stackoverflow.com/a/3339266/580010
	process = subprocess.Popen(['R', 'CMD', 'BATCH', 'graphFile.R'], cwd=climb.TEMPFOLDER)
	process.wait()
	print '<a name="graph"><h2>Mutation Distribution Graph</h2></a><br>'
	print headerLinks
	print '<br><img src="%s">' % ('/tmp/MRT-' + randomToken + '/graphFile.png')
	print '<br><a href="%s">Download</a> raw graph data in CSV format<br><br>' % ('/tmp/MRT-' + randomToken + '/graphFile.csv')

########################################
### Table 3: Summary table -- motifs ###
########################################
print '<a name="motifdistribution"><h2>Motif Distribution</h2></a><br>'
print headerLinks
print '<table border="1" style="font-family:monospace; font-size:16px">'
M = S.extract(loci, mapping=int(anchorPosition))
motifList = {}	# stores the counts of each motif; key is motif
motifSeq = {}	# stores the sequence IDs for each motif in a list; key is motif
motifLabel = {}	# stores the short label for each motif; key is motif
motifLabelCount = 1	# start labels from 1
for i in M:
	if motifList.has_key(i[1]):
		motifList[i[1]] += 1
		motifSeq[i[1]].append(i[0])
	else:
		motifList[i[1]] = 1
		motifSeq[i[1]] = [i[0]]

# pass in the string conversion of the first returned value of the loc -- str(climb.parseList(loci)[0][0])
# passing '0' instead
# update: passing in link to details instead

# 16 January 2013
# The climb.distributionShow method is used to display various output tables
# and can output both plain text and HTML; however, modifying this routine to
# include a new left-most column containing the motif labels is not a feasible
# solution; instead, simply output the motif distribution directly

######################################################################
######################################################################
######################################################################

##### print climb.distributionShow([motifList], '0',  len(M),  motifList.keys(),                        html=True, suppressZeros=sZ, percentage=perc, firstCell="<b>Motif</b>")
##### def         distributionShow(D,           loci, numSeqs, rowList=BASES+NONBASES, percentage=True, html=False, suppressZeros=False, firstCell=''):

print '<tr><td><b>Label</b></td><td><b>Motif</b></td><td><b>Distribution</b></td></tr>'

tempSort = []
D = [motifList]	# to match variables from original method
for kk, vv in D[0].iteritems():	# extract the dictionary from the list
	tempSort.append([kk, vv])
tempSort = sorted(tempSort, key=lambda element: (-element[1], element[0]))	# ascending order of "value, key"; - for descending sort
rowList = []
for tempLoop in tempSort:
	rowList.append(tempLoop[0])

out2 = ''
for i in rowList:	# loop through each base; force the order
	motifLabel[i] = 'L%04i' % motifLabelCount	# if key not present, add key with label
	motifLabelCount += 1
	tt1 = ''
	for tt2 in i:
		tt1 = tt1 + tt2 # + '\t'
	out2 += '<tr><td><b>%s</b></td><td>%s</b>' % (motifLabel[i], tt1) 	# row heading
	for j in D:
		out2 += '</td><td>'
		if (not sZ) or (j[i] != 0):	# populate table if value is not zero or if suppressZeros is false
			if perc:
				out2 += '%06.2f' % (j[i] / float(len(M)) * 100)	# too many leading zeros? accommodating '100%'
			else:
				out2 += '%3i' % (j[i])	# leading zeros untidy

		else:
			out2 += '&nbsp;'

	out2 += '</tr>\n'

print out2

######################################################################
######################################################################
######################################################################

print '</table>'

# copied from code above for optionally plotting the mutation distribution
# plots percentages or counts of motifs, in descending order
import csv
graphFile = open(climb.TEMPFOLDER + 'motifDistributionGraphFile.csv', 'wb')
graphOut = csv.writer(graphFile)

if perc:
	divisor = len(M) / float(100.0)
	yLabel = 'Percentage'
else:
	divisor = 1
	yLabel = 'Count'

for  i in motifList:
	graphOut.writerow([motifLabel[i], i, motifList[i]/float(divisor)])
graphFile.close()




# need parameters for this graph?
# sanity checks on the graph parameters
if not graphWidth.isdigit():
	graphWidth = '480';
if not graphHeight.isdigit():
	graphHeight = '480';
graphTitle = graphTitle.replace('\\', ' ')

fontSize1 = int(int(graphWidth) / 54)
fontSize2 = int(int(graphWidth) / 60)





graphFile = open(climb.TEMPFOLDER + 'motifDistributionGraphFile.R', 'w')
graphFile.write("library('ggplot2')\n")
graphFile.write("v = read.csv('motifDistributionGraphFile.csv', header=F)\n")	# V1 is label, V2 is motif, V3 is value
graphFile.write("v$V4 = reorder(as.character(v$V1), -v$V3)\n")	# -v$V2 for reverse (descending) sorting
graphFile.write("png('motifDistributionGraphFile.png', width=%s, height=%s, units='px')\n" % (graphWidth, graphHeight))
graphFile.write("ggplot(v, aes(x=v$V4, y=v$V3), stat='identity') + geom_bar(fill='dark red', stat='identity', alpha=1.0) + labs(x = 'Motif Label',y = '%s') + opts(title = 'Motif Distribution', plot.title = theme_text(size = '24'))\n" % yLabel)
# graphFile.write("p = barplot(as.vector(x1$V2), names.arg=as.vector(x1$V1))\n")
graphFile.write("dev.off()\n")
graphFile.close()

import subprocess
# http://stackoverflow.com/a/3339266/580010
process = subprocess.Popen(['R', 'CMD', 'BATCH', 'motifDistributionGraphFile.R'], cwd=climb.TEMPFOLDER)
process.wait()
# print '<br><a name="graph"><h2>Motif Distribution Graph</h2></a><br>'
# print headerLinks
print '<br><img src="%s">' % ('/tmp/MRT-' + randomToken + '/motifDistributionGraphFile.png')
print '<br><a href="%s">Download</a> raw graph data in CSV format<br><br>' % ('/tmp/MRT-' + randomToken + '/motifDistributionGraphFile.csv')




######################################
### Table 4: Motif summary details ###
######################################
detailsFilename = '/tmp/MRT-' + randomToken + '/details.html'
detailsFile = open('/var/www' + detailsFilename, 'w')
detailsLink = '<a href="%s" TARGET="_blank">View</a>' % detailsFilename

print '<br>%s detailed output<br>' % (detailsLink)
print '<br>'

detailsFile.write('<html><head><title>Detailed output for motif summary</title></head><body>')

detailsFile.write('<table border="1" style="font-family:monospace; font-szie:16px">')

tempID = ''
tempBG1 = 'lightblue'
tempBG2 = 'lightgreen'
tempBG = tempBG1
for i in motifSeq:	# i is the key
	for j in motifSeq[i]:
		if tempID != i:
			if tempBG == tempBG1:
				tempBG = tempBG2
			else:
				tempBG = tempBG1
			tempID = i
		tt1 = ''
		for tt2 in i:
			tt1 = tt1 + tt2 + '\t'
		detailsFile.write('<tr><td bgcolor="%s">%s</td><td bgcolor="%s">%s</td></tr>' % (tempBG, tt1, tempBG, j))
detailsFile.write('</table></body></html>')
detailsFile.close()

# if not aa:
if False:
	print '<table border="1" style="font-family:monospace; font-size:16px">'
	TT = S.translateLoci(loci, translateMapping = int(anchorPosition))
	print climb.translateLociShow(TT, html=True, highlightAminoAcids=highlightAA)
	print '</table>'
	print '<br>'
S.unload(override=True)

closePage()

