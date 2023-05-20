import subprocess
import sys
import os
import re

"""
Author: Tyler Cernik
Date: 1/16/2023
This program will be used to create autocommit messages of all the commits that have been made to clubb since the last autocommit.
This is done over a few steps:
    1. First we will read in arguments from the command line where this script is run from:
        - The path to the host repo that we will be accessing to view the last autocommit in it
        - The path to the internal directory of clubb that we will be looking for commits in
    
    2. Then we will put all commit information into a dictionary where the key is the commit hash and the value is the entire commit
    3. Then we will get the commit hash of the newest clubb commit in the host repo's latest autocommit message
    4. Then we will get all the commits in clubb that are newer than the commit hash we got in step 3


TODO:
    - set up program to be able to take in an argument for which repo we are getting the autocommit messages from
       (wrf, e3sm, cam) (or maybe just the path to the repo)
    - have the program find the latest autocommit in that repo and get the latest clubb commit from it.
    - Only grab clubb commits newer than the clubb commit in the latest autocommit
    - make sure to only add commits that have files changed that we would actually have change(driver code, etc.)
    - create the commmit message for the new autocommit and return it to standard output.
"""

external_repo = "" # This will be the path to the repo that we are getting the previous autocommit messages from
clubb_internal_directory = "" # We will look for commits that update files in this repo inside of clubb
log_dict = {} # This will be a dictionary that will hold the commit hash as the key and the entire commit as the value

# clubb_internal_directory = "src/CLUBB_core"

def get_arguments():
    """
    This method will get the arguments from the command line and set the values of the variables
    external_repo and clubb_interal_directory
    """

    # if there is no argv[2], then we will assume that the internal directory is src/CLUBB_core

    # if there is no argv[1], we will throw an error and exit the program, returning an empty string as the message

    if(len(sys.argv) > 1 and os.path.isdir(sys.argv[1])):
        # if sys argv[1] is not empty and is a valid directory, then we will set the external_repo to that directory
        external_repo = sys.argv[1]
        if(sys.argv[2]):
            clubb_internal_directory = sys.argv[2]
        else:
            print("Path to internal directory is not valid, defaulting to src/CLUBB_core")
            clubb_internal_directory = "src/CLUBB_core"
    else:
        print("the path to the repo you want to get the autocommit messages from is incorrect, exiting")
        sys.exit(1)

    # print("external_repo: " + external_repo)
    # print("clubb_internal_directory: " + clubb_internal_directory)

    return external_repo, clubb_internal_directory


def run(external_repo, clubb_internal_directory, log_dict):
    """
    This method will get the git log from the clubb repo and put it into a dictionary where the key is the commit hash and the value is the entire commit
    """
    # This method will be called from the base of the clubb repo, otherwise the clubb_internal_directory will need to be changed
    # to reflect as such.

    # Git log pretty format that will get the git hash (%H) and the commit message (%b)
    format = "--pretty=\"Start_of_Commit\nCommit %H\nAuthor: %an\nDate: %ad\n%s\n%b\""

    # get the git log from the clubb repo
    git_log = subprocess.run(["git", "log", format, "--", clubb_internal_directory], stdout=subprocess.PIPE)

    # decode the git log from bytes to a string, and split it into a list of strings
    # git_log = list(map(lambda x: x.strip("\n") ,str(git_log.stdout.decode("utf-8")).split("commit_hash:")))

    git_log = list(str(git_log.stdout.decode("utf-8")).strip("\"").split("Start_of_Commit\n"))

    # clean all empty strings from the list
    git_log = list(filter(None, git_log))

    for commit in git_log:
        hash = re.search('Commit (.+)', commit).group(1)
        log_dict[hash] = commit

    # Testing code - writes all items in the git_log list to a file
    # out_file = open("/home/cernikt/Documents/autocommit_message_maker/clubb/utilities/git_log.txt", "w")
    # # index = git_log.index()
    # for line in git_log:
    #     out_file.write(line)
    #     out_file.write("\n")
    # out_file.close()


def get_autommit_message_commit(external_repo, clubb_internal):
    """
    This method will get the latest autocommit message from the external repo and return the commit hash of the latest clubb commit
    If it cannot find a clubb commit in the latest autocommit message, it will return -1
    """
    commit_format = "--pretty=\"commit_hash:%H %s Body: %b\""
    # get the latest autocommit message
    # Gets it in the format commit_hash:hash commit_title
    # get the latest clubb commit from that autocommit

    #print(clubb_internal.split("/"))

    clubb_internal_src = clubb_internal.split("/")[1].strip()

    #print(external_repo)

    os.chdir(external_repo)

    
    repo_git_log = subprocess.run(["git", "log", commit_format], stdout=subprocess.PIPE)
    repo_git_log = str(repo_git_log.stdout.decode("utf-8")).split("commit_hash:")
    repo_git_log = list(filter(None, repo_git_log))


    #repo_out_file = open(("/home/cernikt/Documents/autocommit_message_maker/clubb/utilities/repo_git_log.txt"), "w")
    # index = git_log.index()
    #for line in repo_git_log:
      #  repo_out_file.write(line)
    # repo_out_file.close()

    hash = ""
    for line in repo_git_log:
        if ("Autoupdated " + clubb_internal_src) in line:

            # this method looks through the git log of the host repo and finds the latest autocommit message


            # print(line)
            # print(len(line))
            # this offset is here to account for the length of the "Autoupdated src/CLUBB_core" string
            if(len(line) <= (len("Autoupdated " + clubb_internal_src)+80)):
                return -1
            offset = line.index("Autoupdated " + clubb_internal_src) + len("Autoupdated " + clubb_internal_src)

            # 54 is the length of the string "Commit " plus the length of the hash
            hash = re.search('Commit (.+)', line[:offset + 54]).group(1).strip()
            # print(hash)
            # once we find the latest autocommit, we can break out of the loop
            return hash
    
    # do a git show on the hash to get the commit message,
    # then get the latest clubb commit from that autocommit message

    return -1


def get_new_clubb_commits(log_dict, external_repo, clubb_internal_directory):
    """
    This method will run previously defined methods to get the latest autocommit message from the external repo
    This is done by first calling teh get autommit message commit method, which will return the commit hash of the latest clubb commit in the host repo
    Then, based on what that is, we will return a string that contains the commit message for the newest auto commit message,
    which contains the commit messages that make up that autocommit.
    """
    # get the latest clubb commit from the latest autocommit message in the external repo

    internal = clubb_internal_directory.split("/")[1].strip()

    run(external_repo, clubb_internal_directory, log_dict)

    hash = get_autommit_message_commit(external_repo, clubb_internal_directory)
    # print(hash)
    if(hash == -1):
        # this case occurs when there are no clubb commit message in the latest autocommit message in the host repo
        # This occurs when either:
        # 1. The external repo is not properly specified in the arguments (This likely won't happen, since this script is run in another script)
        # 2. There are no autocommit messages in the external repo (This can occur if it is the first time this script is run during the autocommit process)

        # In this case, we will just return a message that contains the latest clubb commit message that was made.
        # Because of this, once the script is run again, it will be able to find the latest autocommit message, and will be able to get the latest clubb commit message from that
        # it also adds a message that says that there are no autocommit messages found, and that the latest clubb commit message was printed instead

        key = list(log_dict.keys())[0]
        temp = (("Autoupdated " + internal + "\n") + ("").join(list(log_dict[key])).replace('"', ''))
        
        temp += ("No autocommit messages found, printed latest commit from " + internal + " instead. \n Please make sure the that there are previous autocommit messages, or that the external repo is properly configured.")
        return temp
    else:
        # We did find a clubb commit in the latest autocommit message in the host repo, so we can proceed as normal

        # get all of the clubb commits in the internal repo
    
        # filter the list of clubb commits to only include commits that are after the latest autocommit clubb hash
        new_clubb_commits = []

        # Gettting all of the clubb commits that are after the latest autocommit clubb commit
        for key in list(log_dict.keys()):

            if(key in hash):
                return ("Autoupdated " + internal + "\n") + ("\n").join(new_clubb_commits).replace('"', '')
            new_clubb_commits.append(log_dict[key])

        # return the list of new clubb commits as a string
        return ("Autoupdated " + internal + "\n") + ("\n").join(new_clubb_commits).replace('"', '')


if __name__=="__main__":

    try:
        external_repo, clubb_internal_directory = get_arguments()
        out = get_new_clubb_commits(log_dict, external_repo, clubb_internal_directory)

        sys.stdout.flush()

        print(out)
    except Exception as e:
        # print("Caught exception: {}".format(e))
        print("Autoupdated " + clubb_internal_directory.split("/")[1].strip())

    # final_out_file = open(("/home/cernikt/Documents/autocommit_message_maker/clubb/utilities/final_git_log.txt"), "w")
    # index = git_log.index()
    # for line in out:
    #     final_out_file.write(line)
    #     # final_out_file.write("\n")
    # final_out_file.close()
