#!/usr/bin/env bash
#
# Cycles tcpflow output through foremost to recover jpg files,
# then executes python script to generate gallery pages and index.
#
# For example, capture HTTP traffic and pass it to tcpflow,
#   cd ~/pcap && tcpdump -U tcp port 80 -w - | tcpflow -r -
#
# Then in another shell, run this script,
#   ./pcap.sh
#
# You may need to use sudo if writing to a protected directory such as /var/www.

set -o errtrace
set -o nounset
set -o pipefail

foremost=$(command -v foremost)  # path to foremost
python=$(command -v python)      # path to python
gallery=/home/drduh/gallery.py   # path to python gallery script
output=/home/drduh/www           # path to write jpgs
tcpflow=/home/drduh/pcap         # path to where tcpflow is writing
pause=60                         # seconds between each run


sanity_check() {
  # Ensure foremost and required directories exist and are sane.

  if [[ -z ${foremost} && ! -x ${foremost} ]] ; then
    echo "foremost is not available" ; exit 1
  fi

  if [[ ! -d ${tcpflow} && "${tcpflow}" != "/" ]] ; then
    echo "${tcpflow} does not exist" ; exit 1
  fi

  if [[ ! -d ${output} ]] ; then
    echo "${output} does not exist" ; exit 1
  fi

  if [[ ! -f ${gallery} ]] ; then
    echo "${gallery} does not exist" ; exit 1
  fi
}


sanity_check


while : ; do
  date

  # Clear out any previous foremost output
  rm -r "${tcpflow}/output" 2>/dev/null

  if [[ "$(ls -A ${tcpflow})" ]] ; then
    cd ${tcpflow} && \
      ${foremost} -i * &>/dev/null

    if [[ -d ${tcpflow}/output/jpg ]] ; then
      echo "Found $(ls -l ${tcpflow} | wc -l) jpgs from foremost"

      for filename in $(find ${tcpflow}/output/jpg/ -type f -name "*.jpg") ; do
        # add random number to filename as it can repeat in foremost output
        cp $filename ${output}/$RANDOM-$(basename $filename)
      done

      echo "Updating gallery pages ..."
      ${python} ${gallery}

    else
      echo "No new jpgs from foremost"
    fi

    # Clear out any tcpflow for next run
    rm -r ${tcpflow}/* &>/dev/null
  else
    echo "No new clubb from tcpflow"
  fi

  echo "Sleeping for ${pause} seconds" ; sleep ${pause}
  echo ""
done

