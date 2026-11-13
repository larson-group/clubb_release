"""Variables for gallery.py."""
import os

root = os.path.dirname(os.path.realpath(__file__)) +'/../output'  # path to jpgs or folders of jpgs and output root
tmp = '/tmp'              # temporary folder to move corrupt files to
index = 'index.html'      # filename for html files
plots_filename = 'plots.html'
navigation = 'navigation.html'

header = ("""<!doctype html>
<html>
<head>
  <title>PyPlotgen Output</title>
  <meta charset="utf-8" />
  <meta http-equiv="Content-type" content="text/html; charset=utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <style type="text/css">
    body {
      background-color: white;
      font-family: "Open Sans", "Helvetica Neue", Helvetica, Arial, sans-serif;
      margin: 0;
      padding: 0;
    }
    div {
      border-radius: 0.25em;
      margin: 1em auto;
      padding: 2em;
      text-align: center;
    }
    h1{
      font-size: 48px;
    }
    p {
      font-size: 16px;
      padding-bottom: 1.5em;
    }
    a:link, a:visited {
      font-size: 24px;
      text-decoration: underline;
    }
    img {
      padding: 0.1em;
      border-radius: 0.25em;
    }
  </style>
</head>
<body>
<div align="CENTER"><h1>PyPlotgen Output</h1></div>
<div>
""")

nav_header = '<html>\n<h4>Cases:</h4>'
nav_footer = '</html>'
nav_a_tag = '<a href="plots.html#%s" target="plots">%s</a><br/>'

idx_page = '<html>\n\t' \
            '<title>PyPlotgen</title>\n\t' \
            '<frameset cols="180,*">\n\t' \
            '<frame src="navigation.html" frameborder="0" name="nav">\n\t' \
            '<frame src="plots.html" frameborder="0" name="plots">\n\t'\
            '</frameset>\n\t'\
            '</html>'

a_tag = '<div align="CENTER"><a name="%s"><font size="+6"><b>%s</b></font></a></div>'
br = '\n<br>'

setup_file_link = '<a href="%s">Case setup information for clubb folder %s</a>' + br

br = '\n<br>'
footer = '\n</div></body></html>'
nav = '%s</div>\n<div>'
img_src = '\n<img src="%s">'
case_title = '\n<h1>%s</h1>'
case_description='\n<p>%s</p>'
timestamp = '\n<p>Page created on %s</p>'
url_dir = '\n<p><a href="%s">%s</a></p>'
url_img = '\n<a href="%s"><img title="%s" src="%s"></a>'
url_mov = '\n<video controls><source src="%s"></video>'
