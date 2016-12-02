import os
import logging

#-------------------------------------------------------------------------------
#   L O G G E R
#-------------------------------------------------------------------------------
FORMAT='%(asctime)s:%(levelname)s:%(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger('OutputWriter')
#logger.setLevel(logging.DEBUG)

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
# Writes sub text under the header
###############################################################################
def writeSubHeader(FILE, TEXT):
    logger.info('writeSubHeader')
    with open(FILE, "a") as myfile:
        myfile.write("""	<div align="CENTER">
	        <b><font size="-1" color="#430e9a"> """ + TEXT + """</font></b>
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
# Writes the SAM_CLUBB Variable Equivalence Table
###############################################################################
def writeSamSubHeader(FILE):
    logger.info('writeSamSubHeader')
    writeSubHeader(FILE, "2D SAM_CLUBB runs use a 64-km horizontal domain and a 10-s timestep, with CLUBB called every SAM timestep. 3D SAM_CLUBB runs use 4x4 columns with a 16-km horizontal grid spacing and a 10-s timestep with CLUBB called every SAM timestep. All SAM_CLUBB runs except LBA use Morrison microphysics. CLUBB standalone runs use a 10-s timestep and the Morrison microphysics.")
    writeSubHeader(FILE, "When two variables are listed, the first variable is the SAM-CLUBB variable and the second is the SAM-Standalone variable.")
    text = """    <br />
    <DIV ALIGN="CENTER"><TABLE CELLPADDING=3 BORDER="1">
    <TR>
        <TD ALIGN="CENTER" COLSPAN=11><B>Variable Equivalence Table</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER"><B>CLUBB</B></TD>
        <TD ALIGN="CENTER"><B>SAM CLUBB</B></TD>
        <TD ALIGN="CENTER"><B>SAM Standalone</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">wpthlp</TD>
        <TD ALIGN="CENTER">wpthlp+tlflux</TD>
        <TD ALIGN="CENTER">tlflux</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">wprtp</TD>
        <TD ALIGN="CENTER">wprtp+qtflux</TD>
        <TD ALIGN="CENTER">qtflux</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">thlp2</TD>
        <TD ALIGN="CENTER">thlp2+tl2</TD>
        <TD ALIGN="CENTER">tl2</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">rtp2</TD>
        <TD ALIGN="CENTER">rtp2 + qt2</TD>
        <TD ALIGN="CENTER">qt2</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">upwp</TD>
        <TD ALIGN="CENTER">uw + upwp</TD>
        <TD ALIGN="CENTER">uw</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">vpwp</TD>
        <TD ALIGN="CENTER">vw + vpwp</TD>
        <TD ALIGN="CENTER">vw</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">up2</TD>
        <TD ALIGN="CENTER">up2 + u2</TD>
        <TD ALIGN="CENTER">u2</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">vp2</TD>
        <TD ALIGN="CENTER">vp2 + v2</TD>
        <TD ALIGN="CENTER">v2</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">wp2</TD>
        <TD ALIGN="CENTER">wp2 + w2</TD>
        <TD ALIGN="CENTER">w2</TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER">wp3</TD>
        <TD ALIGN="CENTER">wp3 + w3</TD>
        <TD ALIGN="CENTER">w3</TD>
    </TR>
</TABLE>
<br />
</DIV>"""
    writeSubHtml(FILE, text)
  
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
    
############
def writeMorrBudgetSubHeader(FILE):
    logger.info('writeMorrBudgetSubHeader')
    text = """    <br />
    <br />
    <br />
    <DIV ALIGN="CENTER"><TABLE CELLPADDING=3 BORDER="1">
    <TR>
        <TD ALIGN="CENTER" COLSPAN=11><B>Microphysical Processes</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=11><B>Processes Affecting Mixing Ratios (kg/kg/s)</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>EPRD</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>sublimation of cloud ice</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>PRACS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of rain by snow to form snow</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>EPRDG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>sublimation of graupel</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>PRAI</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of cloud ice by snow</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>EPRDS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>sublimation of snow</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>     </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>PRC</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>autoconversion</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>EVPMG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>evaporation of melted graupel</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>     </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>PRCI</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>autoconversion of cloud ice to snow</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>EVPMS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>evaporation of melted snow</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>PRD</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>depositional growth of cloud ice</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>MNUCCC</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>contact freezing of cloud droplets</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>PRDG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>depositional growth of graupel</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>MNUCCD</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>freezing of aerosol</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>PRDS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>depositional growth of snow</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>MNUCCR</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>contact freezing of rain droplets</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>PRE</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>rain evaporation</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>PCC</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>saturation adjustment</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>PSACR</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of snow by rain to form graupel</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>PGMLT</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>freezing of rain to form graupel</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>PSACWG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of cloud water by graupel</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>PGRACS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of rain by snow to form graupel</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>PSACWI</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of cloud water by cloud ice</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>PGSACW</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>riming of cloud water by snow to form graupel</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>PSACWS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of cloud water by snow</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>PIACR</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of cloud ice by rain to form graupel</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>PSMLT</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>freezing of rain to form snow</B></TD>
    </TR>
    <TR>writeIndex
        <TD ALIGN="CENTER" COLSPAN=2><B>PIACRS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of cloud ice by rain to form snow</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>QMULTG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>splintering from droplets accreted onto graupel</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>PRA</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>accretion</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>QMULTR</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>splintering from rain droplets accreted onto snow</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>PRACG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of rain by graupel</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>QMULTRG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>splintering from rain droplets accreted onto graupel</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>PRACI</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of cloud ice by rain to form graupel</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>QMULTS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>splintering from cloud droplets accreted onto snow</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>PRACIS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of cloud ice by rain to form snow</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>QC_INST</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>changes due to instantaneous processes</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>QG_INST</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>changes due to instantaneous processes</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>QI_INST</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>changes due to instantaneous processes</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>QR_INST</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>changes due to instantaneous processes</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>QS_INST</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>changes due to instantaneous processes</B></TD>
    </TR>

    <TR>
	<TD ALIGN="CENTER" COLSPAN=11><B>                 </B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=11><B>Processes Affecting Number Concentrations (#/kg/s)</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NACT</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>cloud drop formation by aerosols</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NPRC1</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>change in rain due to autoconversion</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NCSTEN</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>cloud water sedimentation</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NPRCI</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>autoconversion of cloud ice to snow</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NEGFIX_NC</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>removal of negative cloud water concentration</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NPSACWG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of cloud drops by graupel</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NEXFIX_NG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>removal of negative graupel concentration</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NPSACWI</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>cloud droplet accretion by cloud ice</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NEGFIX_NI</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>removal of negative ice concentration</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NPSACWS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>cloud droplet accretion by snow</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NEGFIX_NR</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>removal of negative rain concentration</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NRAGG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>self collection of rain drops</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NEGFIX_NS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>removal of negative snow concentration</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NRSTEN</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>rain sedimentation</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NGMLTG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>melting of graupel</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NSAGG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>aggregation of snow</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NGMLTR</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>melting of graupel to form rain</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NSCNG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>conversion of snow to graupel</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NGRACS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of rain by snow</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NSMLTR</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>melting of snow to form rain</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NGSTEN</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>graupel sedimentation</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NSMLTS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>melting of snow</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NIACR</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of ice by rain to form graupel</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NSSTEN</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>snow sedimentation</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NIACRS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of ice by rain to form snow</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NSUBG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>sublimation of graupel</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NIM_MORR_CL</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>clipping of large ice concentration</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NSUBI</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>sublimation of cloud ice</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NISTEN</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>cloud ice sedimentation</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NSUBR</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>evaporation of rain</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NMULTG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>splintering due to accretion of droplets by graupel</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NSUBS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>sublimation of snow</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NMULTR</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>splintering due to rain riming snow</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NC_INST</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>changes due to instantaneous processes</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NMULTRG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>splintering due to rain riming snow</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NG_INST</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>changes due to instantaneous processes</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NMULTS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>splintering due to accretion of droplets by snow</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NI_INST</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>changes due to instantaneous processes</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NNUCCC</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>contact freezing</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NR_INST</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>changes due to instantaneous processes</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NNUCCD</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>ice nucleation by freezing aerosol</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>NS_INST</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>changes due to instantaneous processes</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NNUCCR</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>contact freezing of rain</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>SIZEFIX_NC</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>adjust cloud number when size is too large/small</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NPRA</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>accretion</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>SIZEFIX_NG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>adjust graupel number when size is too large/small</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NPRACG</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of rain by graupel</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>SIZEFIX_NI</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>adjust graupel number when size is too large/small</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NPRACS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>collection of rain by snow</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>SIZEFIX_NR</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>adjust rain number when size is too large/small</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NPRAI</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>accretion of cloud ice by snow</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>SIZEFIX_NS</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>adjust snow number when size is too large/small</B></TD>
    </TR>
    <TR>
        <TD ALIGN="CENTER" COLSPAN=2><B>NPRC</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>autoconversion</B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B>    </B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B></B></TD>
        <TD ALIGN="CENTER" COLSPAN=2><B></B></TD>
    </TR>


</TABLE>
<BR/>
<BR/>
</DIV>"""
    writeSubHtml(FILE, text)
    
############
def writeWrfHeader(FILE):
    logger.info('writeWrfHeader')
    writeSubHeader(FILE, "WRF_CLUBB runs use a 100-km horizontal grid spacing, 3x3 grid columns, and a 60-s timestep.")
    
    text = """<br />
<br />
"""
    
    writeSubHtml(FILE, text)
    
############
def writeGfdlHeader(FILE):
    logger.info('writeGfdlHeader')
    writeSubHeader(FILE, "GFDL-SCM-CLUBB and CLUBB-standalone runs use a 10-min timestep with 48 levels")
    
    text = """<br />
<br />
"""
    
    writeSubHtml(FILE, text)

############
def writeFooter(FILE):
    logger.info('writeFooter')
    with open(FILE, "a") as myfile:
        myfile.write("""	<br /> <br /> <br /> <br />
	<hr noshade size=5 width=70%>
	<div align="CENTER">
		<font size="-2">
		Copyright &#169; 2016 Larson Group. All rights reserved. 
		</font>
	</div>
</body>
</html>
""")

###############################################################################
# Inserts an image
###############################################################################
def placeImage(FILE, img):
    logger.info('placeImage')
    imgWidth = '49%'
    imgHeight = 'auto'
    with open(FILE, "a") as myfile:
        myfile.write("		<img width=\"" + imgWidth + "\" height=\"" + imgHeight + "\" align=\"BOTTOM\" border=\"0\" style=\"padding: 5px;\" src=\"" + img + "\" alt=\"" + img + "\" />\n")

############
def printDivCenter(FILE):
    logger.info('printDivCenter')
    with open(FILE, "a") as myfile:
        myfile.write("""	<div align="CENTER">
""")

############
def printCloseDivCenter(FILE):
    logger.info('printCloseDivCenter')
    with open(FILE, "a") as myfile:
        myfile.write("""	</div>
""")

##################################
def getListOfImages(directory):
    logger.info('getListOfImages')
    fileList = os.listdir(directory)
    fileList.sort()
    return fileList

def writeNavPage(directory, case):
    nav = directory + 'navigation.html'
    writeNavPageStart(nav)
    writeNavPageCase(nav, case, case)
    writeNavPageClose(nav)

def writePlotsPage(directory, case, mode):
    plots = directory + 'plots.html'
    images = directory + 'jpg/'
    
    writeHeader(plots, mode)
    writeSamBudgetSubHeader(plots)
    writeCaseTitle(plots, case)
    imageList = getListOfImages(images)
    logger.debug(imageList)
    for image in imageList:
        placeImage(plots, './jpg/' + image)
    writeFooter(plots)
