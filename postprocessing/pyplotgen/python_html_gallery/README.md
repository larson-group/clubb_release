# _CLUBB Users_

You should be able to ignore this folder!

This folder is included as a dependency for generating a webpage to display pyplotgen output images.  
Please see their documentation from the git repo below. 

----
# python-html-gallery

Python script to generate basic HTML image gallery pages from jpgs.

# Installation

Requires Python Imaging Library - `PIL` or `Pillow`

Install with `pip install PIL` or `easy_install PIL` or build and install it from [source](http://www.pythonware.com/products/pil/).

Then `git clone https://github.com/drduh/python-html-gallery`

# Use

Put jpgs or folders of jpgs into a directory specified in `static.py`

Then run `python gallery.py` to create thumbnails, gallery pages and an index page.

`pcap.sh` is an example shell script for generating jpgs from `tcpflow`/`foremost` for gallery content.
