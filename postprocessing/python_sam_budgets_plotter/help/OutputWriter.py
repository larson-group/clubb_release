import os
import logging
from datetime import datetime as dt

#-------------------------------------------------------------------------------
#   L O G G E R
#-------------------------------------------------------------------------------
logger = logging.getLogger('plotgen.help.ow')
#logger.setLevel(logging.INFO)
logger.setLevel(logging.DEBUG)
#logger.setLevel(logging.CRITICAL)

#-------------------------------------------------------------------------------
#    F U N C T I O N S
#-------------------------------------------------------------------------------
###############################################################################
# Writes the index page that contains the two frames.
###############################################################################
def writeIndex(FILE, mode):
    logger.info('writeIndex')
    with open(FILE, "w+") as myfile:
        myfile.write("""    <html>
        <title>""" + mode + """</title>
        <frameset cols="180,*">
            <frame src="navigation.html" frameborder="0" name="nav">
            <frame src="plots.html" frameborder="0" name="plots">
        </frameset>
    </html>
""")

###############################################################################
# Writes the beginning of the navigation page HTML. This should be called 
# before the cases are written to the navigation page.
###############################################################################
def writeNavPageStart(FILE):
    logger.info('writeNavPageStart')
    with open(FILE, "w+") as myfile:
        myfile.write("""    <html>
        <h4>Cases:</h4>
""")

###############################################################################
# Writes a case to the navigation page
###############################################################################
def writeNavPageCase(FILE, case, link):
    logger.info('writeNavPageCase')
    with open(FILE, "a") as myfile:
        myfile.write("""       <a href="plots.html#""" + link + """" target="plots">""" + case + """</a><br/>
""")
    
###############################################################################
# Writes HTML to the end of the navigation page. This should be called after
# all cases were added.
###############################################################################
def writeNavPageClose(FILE):
    logger.info('writeNavPageClose')
    with open(FILE, "a") as myfile:
        myfile.write("""    </html>
""")
    
###############################################################################
# Writes a case title to the HTML file
###############################################################################
def writeCaseTitle(FILE, CASE):
    logger.info('writeCaseTitle')
    with open(FILE, "a") as myfile:
        myfile.write("	<a name=\"" +  CASE + """" ></a>
	<div align="CENTER">
		<font size="+2">
		<font color="#0000ff"> <a href="#"""+ CASE + """">""" + CASE + """</a> </font> </font>
	</div>
""")

###############################################################################
# Writes HTML under the case header
###############################################################################
def writeSubHtml(FILE, HTML):
    logger.info('writeSubHtml')
    with open(FILE, "a") as myfile:
        myfile.write("""	<div align="CENTER">
		""" + HTML + """
	</div>
""")

###############################################################################
# Writes the HTML header information for CLUBB
###############################################################################
def writeHeader(FILE, mode):
    logger.info('writeHeader')
    from time import strftime
    with open(FILE, "w+") as myfile:
        myfile.write("""<html>
<head>
	<title>""" + mode + """</title>
	<p><div align="CENTER">
		<font size="+3" color="#811212">""" + mode + """</font>
		<font size="-1"><br/>""" + strftime("%m/%d/%Y") + """</font></p>
	</div></p>

</head>
<body>
""")

###############################################################################
# Writes the SAM_CLUBB Budgets Variable Equivalence Table
###############################################################################
def writeSamBudgetSubHeader(FILE):
    logger.info('writeSamBudgetSubHeader')
    text = """    <br />
    <DIV ALIGN="CENTER"><TABLE CELLPADDING=3 BORDER="1">
    <TR>
        <TD ALIGN="CENTER" COLSPAN=11><B>Budget Variable Equivalence Table</B></TD>
    </TR>
    <TR>
	<TD ALIGN="CENTER" COLSPAN=2><B>W2 (wp2)</B></TD>
	<TD ALIGN="CENTER" COLSPAN=2><B>TW/QW (wpthlp/wprtp)</B></TD>
	<TD ALIGN="CENTER" COLSPAN=2><B>T2/Q2 (thlp2/rtp2)</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER"><B>SAM</B></TD>
        <TD ALIGN="CENTER"><B>CLUBB</B></TD>
        <TD ALIGN="CENTER"><B>SAM</B></TD>
        <TD ALIGN="CENTER"><B>CLUBB</B></TD>
        <TD ALIGN="CENTER"><B>SAM</B></TD>
        <TD ALIGN="CENTER"><B>CLUBB</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">ADV (T)</TD>
        <TD ALIGN="CENTER">ta</TD>
        <TD ALIGN="CENTER">ADV (T)</TD>
        <TD ALIGN="CENTER">ta</TD>
        <TD ALIGN="CENTER">ADVTR (T)</TD>
        <TD ALIGN="CENTER">ta</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">PRES (P)</TD>
        <TD ALIGN="CENTER">pr2+pr3</TD>
        <TD ALIGN="CENTER">PRES (P)</TD>
        <TD ALIGN="CENTER">pr1+pr2+pr3</TD>
        <TD ALIGN="CENTER">GRAD (G)</TD>
        <TD ALIGN="CENTER">tp</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">REDIS (R)</TD>
        <TD ALIGN="CENTER">pr1</TD>
        <TD ALIGN="CENTER">GRAD (G)</TD>
        <TD ALIGN="CENTER">tp</TD>
        <TD ALIGN="CENTER">DISSIP (D)</TD>
        <TD ALIGN="CENTER">dp1</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">BUOY (B)</TD>
        <TD ALIGN="CENTER">bp</TD>
        <TD ALIGN="CENTER">BUOY (B)</TD>
        <TD ALIGN="CENTER">bp</TD>
        <TD ALIGN="CENTER">DIFFTR</TD>
        <TD ALIGN="CENTER">dp2</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">DIFF (D)</TD>
        <TD ALIGN="CENTER">dp1+dp2</TD>
        <TD ALIGN="CENTER">DIFF (D)</TD>
        <TD ALIGN="CENTER">dp1</TD>
        <TD ALIGN="CENTER">PREC (PC)</TD>
        <TD ALIGN="CENTER">&nbsp;</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">_RES (Residual)</TD>
        <TD ALIGN="CENTER">&nbsp;</TD>
        <TD ALIGN="CENTER">PREC (PC)</TD>
        <TD ALIGN="CENTER">&nbsp;</TD>
        <TD ALIGN="CENTER">_RES (Residual)</TD>
        <TD ALIGN="CENTER">&nbsp;</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">&nbsp;</TD>
        <TD ALIGN="CENTER">&nbsp;</TD>
        <TD ALIGN="CENTER">_RES (Residual)</TD>
        <TD ALIGN="CENTER">&nbsp;</TD>
        <TD ALIGN="CENTER">&nbsp;</TD>
        <TD ALIGN="CENTER">&nbsp;</TD>
    </TR>
</TABLE>
<br />
<B>Equivalencies above are derived from equations in <a href="http://journals.ametsoc.org/doi/pdf/10.1175/1520-0469%282002%29059%3C2550%3ASODCCC%3E2.0.CO%3B2"> this paper.</a></B>
<br />
<br />
</DIV>"""
    writeSubHtml(FILE, text)

###############################################################################
# Writes the footer of a page
###############################################################################
# NOTE: (Steffen Domke) Replaced static year of 2016 with datetime.datetime.now().year
def writeFooter(FILE):
    logger.info('writeFooter')
    with open(FILE, "a") as myfile:
        myfile.write("""	<br /> <br /> <br /> <br />
	<hr noshade size=5 width=70%>
	<div align="CENTER">
		<font size="-2">
		Copyright &#169; {} Larson Group. All rights reserved. 
		</font>
	</div>
</body>
</html>
""".format(dt.now().year))

###############################################################################
# Inserts an image
###############################################################################
def placeImage(FILE, img):
    logger.info('placeImage')
    imgWidth = '49%'
    imgHeight = 'auto'
    with open(FILE, "a") as myfile:
        myfile.write("		<img width=\"" + imgWidth + "\" height=\"" + imgHeight + "\" align=\"BOTTOM\" border=\"0\" style=\"padding: 5px;\" src=\"" + img + ".png\" alt=\"" + img + "\" />\n")

###############################################################################
# Gets a list of all images in a specific directory
###############################################################################
def getListOfImages(directory, imageNames):
    logger.info('getListOfImages')
    fileList = os.listdir(directory)
    fileList.sort()
    return fileList

###############################################################################
# Calls all navigation page methods
###############################################################################
def writeNavPage(directory, case):
    nav = directory + 'navigation.html'
    writeNavPageStart(nav)
    writeNavPageCase(nav, case, case)
    writeNavPageClose(nav)

###############################################################################
# Calls all plot page methods which are needed
###############################################################################
def writePlotsPage(directory, case, mode, imageNames):
    plots = directory + 'plots.html'
    
    writeHeader(plots, mode)
    writeSamBudgetSubHeader(plots)
    writeCaseTitle(plots, case)
    for image in imageNames:
        placeImage(plots, './jpg/' + image)
    writeFooter(plots)
