#!/usr/bin/env python
#
# Generates thumbnails and HTML gallery pages for jpgs in a directory.

import datetime
import fnmatch
import glob
import multiprocessing
import os
import random
import re
import shutil
import sys
from multiprocessing import Pool
from multiprocessing import freeze_support

from config import Case_definitions
from python_html_gallery import static

try:
    from PIL import Image
except ImportError:
    try:
        import Image
    except ImportError:
        print('Requires Python Imaging Library. See README.md.')
        sys.exit(1)


def ListFiles(regex, path):
    """Returns list of matching files in path."""
    rule = re.compile(fnmatch.translate(regex), re.IGNORECASE)
    return [name for name in os.listdir(path) if rule.match(name)] or None


def ListDirs(path):
    """Returns list of directories in path."""
    return [d for d in os.listdir(path) if os.path.isdir(
        os.path.join(path, d))]


def Now(time=True):
    """Returns formatted current time."""
    if time:
        return datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    else:
        return datetime.datetime.now().strftime('%Y%m%d')


def RandomThumb(page,image_extension=".png"):
    """Returns path to random thumbnail for a given page."""
    return random.choice(
        glob.glob(os.path.join(page.split('/')[0], '*_thumb'+image_extension)))


def OrganizeRoot(output_dir,image_extension=".png"):
    """Creates directories for images in root directory."""

    static.root = output_dir
    try:
        os.chdir(static.root)
    except OSError:
        print('Could not cd into %s' % static.root)
        sys.exit(1)

    fs = ListFiles('*'+image_extension, '.')
    if fs:
        for jpg in fs:
            datehour = Now(time=False)
            if not os.path.exists(datehour):
                print('Creating directory: %s' % datehour)
                os.makedirs(datehour)
            if not os.path.exists(os.path.join(datehour, jpg)):
                shutil.move(jpg, datehour)
            else:
                print('%s already exists' % os.path.join(datehour, jpg))


def GenerateHtmlImages(page, jpgs):
    """
    Creates HTML image elements

  Args:
    page: str, name of page for thumbnails.
    jpgs: list, jpg files to create thumbnails for.
  Returns:
    url_imgs: list, image links to write.
  """
    url_imgs = []
    for jpg in jpgs:
        jpg = page + '/' + jpg
        url_imgs.append(static.url_img % (jpg, jpg, jpg))
    return url_imgs


def get_start_end_minutes(casename):
    """
    Get the start and end time minutes for a case
    as defined in Case_definitions.py
    :param casename: Name of the case as defined by the 'name' parameter of it's entry in CASES_TO_PLOT in Case_definitions.py
    :return: tuple of the order (start_time,end_time)
    """
    for case in Case_definitions.CASES_TO_PLOT:
        if case['name'] == casename:
            return case['start_time'], case['end_time']

def get_description(casename):
    """
    Get the description from the case definintion
    :param casename: Name of the case as defined by the 'name' parameter of it's entry in CASES_TO_PLOT in Case_definitions.py
    :return: The case's description text
    """
    for case in Case_definitions.CASES_TO_PLOT:
        if case['name'] == casename:
            return case['description']


def WriteGalleryPage(page,image_extension=".png"):
    """Writes a gallery page for jpgs in path.

  Args:
    page: str, name of page under root directory.
  """
    # os.chdir(static.root)

    with open(static.plots, 'a') as plots_file:
        # plots_file.write(static.header % page)
        # plots_file.write(static.case_title % page)
        start_time,end_time = get_start_end_minutes(page)
        case_title = page + " minutes " + str(start_time) + "-" + str(end_time)
        plots_file.write(static.a_tag % (page, case_title))
        case_description = get_description(page)
        plots_file.write(static.case_description % (case_description))

        # Write case_setup.txt links
        for setup_file in glob.glob(page+'/*.txt'):
            setup_file_tail_len = 10
            # Find the number of characters that occur in the setupfile's name before the name of the clubb folder.
            # The format for setup filenames is casename/casename_inputfoldername_setup.txt
            casename_prefix_len = len(page) * 2 + 2 # length is num characters in page/page_  including the '/' and '_'
            setup_file_src_folder = setup_file[casename_prefix_len: -setup_file_tail_len]
            plots_file.write(static.setup_file_link % (setup_file, setup_file_src_folder) + "\n")

        plots_file.write(static.timestamp % Now())

        try:
            img_paths = '*'+image_extension
            case_images = ListFiles(img_paths, page)
            jpgs = sorted(case_images, reverse=True)[::-1]
        except TypeError:
            print('%s: No images found' % page)
            return

        for e in GenerateHtmlImages(page, jpgs):
            plots_file.write(e)

        # plots_file.write(static.footer)


def WriteGalleryPages(multithreaded=False,image_extension=".png"):
    """Write gallery pages for directories in root path."""
    with open(static.plots, 'w') as index_file:
        index_file.write(static.header)

    all_pages = sorted(ListDirs(static.root))

    if multithreaded:
        freeze_support()  # Required for multithreading
        n_processors = multiprocessing.cpu_count()
        with Pool(processes=n_processors) as pool:
            pool.map(WriteGalleryPage, all_pages)
    else:
        for page in all_pages:
            WriteGalleryPage(page,image_extension)

    with open(static.plots, 'a') as index_file:
        index_file.write(static.footer)


def WriteNavigation():
    """Write navigation file with gallery links and thumbnails in root path."""
    os.chdir(static.root)

    with open(static.navigation, 'w') as nav_file:
        nav_file.write(static.nav_header)
        # nav_file.write(static.timestamp % Now())

        page_count = 0
        for page in sorted(ListDirs(static.root)):
            page_count += 1
            try:
                nav_file.write(static.nav_a_tag % (page, page))
            except IndexError:
                print('%s: No thumbnails found, removing' % page)
                os.unlink(page)

        nav_file.write(static.nav_footer)

    print('Wrote %s with %s gallery link(s)' % (
        os.path.join(static.root, static.navigation), page_count))

def WriteIndex():
    """

    """
    os.chdir(static.root)

    with open(static.index, 'w') as index_file:
        index_file.write(static.idx_page)

    print("Wrote index.html")

def main(output_dir, multithreaded=False, image_extension=".png"):
    """Main function."""
    OrganizeRoot(output_dir,image_extension)
    WriteGalleryPages(multithreaded=multithreaded,image_extension=image_extension)
    WriteNavigation()
    WriteIndex()

if __name__ == '__main__':
    main()
