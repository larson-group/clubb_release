import subprocess
import sys
import os
import re

"""
Author: Tyler Cernik
Date: 1/16/2023
This program will be used to create autocommit messages of all the commits that have been made to clubb since the last autocommit.
The program will be called from the autocommit script and will return the autocommit message to standard output.
The autocommit script will then take that message and use it to create the autocommit.
This method will also take in the path to the repo that we are getting the autocommit messages from.
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
    # This method will get the arguments from the command line and set the values of the variables
    # external_repo and clubb_interal_directory

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
    #global log_dict
    # This method will be called from the base of the clubb repo, otherwise the clubb_internal_directory will need to be changed
    # to reflect as such.

    # Git log pretty format that will get the git hash (%H) and the commit message (%b)
    format = "--pretty=\"Start_of_Commit\nCommit %H\nAuthor: %an\nDate: %ad\n%s\n%b\""

    # get the git log from the clubb repo
    git_log = subprocess.run(["git", "log", format, "--", clubb_internal_directory], stdout=subprocess.PIPE)
    #"--pretty[format:Commit %H %nMessage: <full-commit-message>], stdout=subprocess.PIPE)
    # print(type(git_log))
    #git_log = git_log.stdout.decode("utf-8")
    #print(type(git_log.stdout))
    #print(str(git_log.stdout).split("\n\n"))


    # decode the git log from bytes to a string, and split it into a list of strings
    # git_log = list(map(lambda x: x.strip("\n") ,str(git_log.stdout.decode("utf-8")).split("commit_hash:")))

    git_log = list(str(git_log.stdout.decode("utf-8")).strip("\"").split("Start_of_Commit\n"))




    # print(len(git_log))


    # clean all empty strings from the list
    git_log = list(filter(None, git_log))


    for commit in git_log:
        # print(commit)
        # print(commit[7:48])
        hash = re.search('Commit (.+)', commit).group(1)
        # commit[7:48].strip("\n")
        log_dict[hash] = commit




    # Testing code - writes all items in the git_log list to a file
    # out_file = open("/home/cernikt/Documents/autocommit_message_maker/clubb/utilities/git_log.txt", "w")
    # # index = git_log.index()
    # for line in git_log:
    #     out_file.write(line)
    #     out_file.write("\n")
    # out_file.close()

def get_autommit_message_commit(external_repo, clubb_internal):
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

            # print(line)
            # print(len(line))
            if(len(line) <= (len("Autoupdated " + clubb_internal_src)+80)):
                return -1
            offset = line.index("Autoupdated " + clubb_internal_src) + len("Autoupdated " + clubb_internal_src)

            hash = re.search('Commit (.+)', line[:offset + 54]).group(1).strip()
            # print(hash)
            # once we find the latest autocommit, we can break out of the loop
            return hash
    
    # do a git show on the hash to get the commit message,
    # then get the latest clubb commit from that autocommit message

    return -1


def get_new_clubb_commits(log_dict, external_repo, clubb_internal_directory):
    # get the latest clubb commit from the latest autocommit message in the external repo

    #global log_dict

    # external_repo, clubb_internal_directory = get_arguments()

    internal = clubb_internal_directory.split("/")[1].strip()

    run(external_repo, clubb_internal_directory, log_dict)

    hash = get_autommit_message_commit(external_repo, clubb_internal_directory)
    # print(hash)
    if(hash == -1):
        key = list(log_dict.keys())[0]
        #print(list(log_dict[key]))
        temp = (("Autoupdated " + internal + "\n") + ("").join(list(log_dict[key])).replace('"', ''))
        
        temp += ("No autocommit messages found, printed latest commit from " + internal + " instead. \n Please make sure the that there are previous autocommit messages, or that the external repo is properly configured.")
        return temp
    else:
        
        # get all of the clubb commits in the internal repo
    
        # filter the list of clubb commits to only include commits that are after the latest autocommit clubb hash
        new_clubb_commits = []

        # print(log_dict.keys())

        # print(log_dict[hash])

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
