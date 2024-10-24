
'''
Creates easily readable output compatible with multiprocessing.
Separates console printing from error log file printing.

:author: Benjamin A. Stephens
:date: November 2020
'''

import os
import logging
import shutil

def logToFile(message):
    """Writes a message to the log file."""
    logging.info("proc_id: "+str(os.getpid())+"\t"+message)


def logToFileAndConsole(message):
    """Writes a message to the console along with log file."""
    print(message)
    logging.info("proc_id: "+str(os.getpid())+"\t"+message)

def initializeProgress(image_extension, movie_extension):
    """Prints the initial progress indicator, updated later"""
    print("Running...\n" + "----------\n" +
          "NOTE: If processing multiple cases (ARM, BOMEX, etc.) with multithreading,\n"+
          "the total number of panels to be processed may increase as threads complete\n"+
          "and begin processing new cases.\n" +
          "----------")
    if movie_extension is not None:
        movie_extension = "."+movie_extension
        print(" --> ANIMATION NOTE: Time slices must sometimes be filtered from animations\n" +
              "     due to mismatched time stepping. Affected cases are noted in error.log.")
        if movie_extension == '.mp4':
            if shutil.which('ffmpeg') is not None:
                print(" --> FFmpeg is present and will be used (see README for more info).")
            else:
                print(" --> FFmpeg not found, videos may not play on all browsers (see README for more info).")
        print("----------")
        print("\rProgress: {:4d} of {:>4} total {} panels complete".format(0,'?',movie_extension), end="")
    else:
        print("\rProgress: {:4d} of {:>4} total {} panels complete".format(0,'?',image_extension), end="")

def updateProgress(total_progress_counter, image_extension, movie_extension):
    """Updates consolve progress continuously to indicate progress to user."""
    percent_complete=int(total_progress_counter[1]/total_progress_counter[0]*100)
    if movie_extension is not None:
        movie_extension = "."+movie_extension
        print("\rProgress: {:4d} of {:4d} total {} panels complete".format(
              total_progress_counter[1],total_progress_counter[0],movie_extension),end="")
    else:
        print("\rProgress: {:4d} of {:4d} total {} panels complete".format(
              total_progress_counter[1],total_progress_counter[0],image_extension),end="")

def writeFinalErrorLog(errorlog,finalerrorlog):
    """Rewrites output file in order if multiple processes were active"""
    #open/read file
    f1 = open(errorlog,"r")
    Lines = f1.readlines()
    f1.close()

    #find length of process ID string
    for i in range(31,len(Lines[0])):
        if Lines[0][i] == "\t":
            proc_end=i
            name_start=proc_end+13
            break

    # extract times and process ids
    times = []
    procs = []
    table = []
    for i in range(len(Lines)):
        times.append(Lines[i][0:20])
        procs.append(Lines[i][31:proc_end])
        if "Processing: " in Lines[i]:
            table.append([Lines[i][name_start:-1],Lines[i][0:20],Lines[i][31:proc_end]])

    table=sorted(table)
    proc_nums=[]
    for i in range(len(table)):
        proc_nums.append(table[i][2])

    #for i in range(len(table)):
    f2 = open(finalerrorlog,'w')

    # write new output file, organized alphabetically according
    # to test case (ARM, etc.) and then time stamp
    written_lines=0
    while procs[written_lines]==procs[0]:
        f2.write(Lines[written_lines])
        written_lines+=1
    for i in range(0,len(proc_nums)):
        if proc_nums.count(proc_nums[i]) > 1:
            other_times=[]
            for j in range(i+1,len(proc_nums)):
                if table[j][2]==table[i][2]:
                    other_times.append(table[j][1])
            if len(other_times) > 0:
                nextmintime=min(other_times)
                for k in range(len(Lines)):
                    if procs[k]==proc_nums[i] and times[k] >= table[i][1] and times[k] < nextmintime:
                        f2.write(Lines[k])
                        written_lines+=1
            else:
                for k in range(len(Lines)):
                    if procs[k]==proc_nums[i] and times[k] >= table[i][1]:
                        f2.write(Lines[k])
                        written_lines+=1
        else:
            for k in range(len(Lines)):
                if procs[k]==proc_nums[i]:
                    f2.write(Lines[k])
                    written_lines+=1
    for i in range(written_lines,len(Lines)):
            f2.write(Lines[i])

    #close new file and remove temp error log
    f2.close()
    os.remove(errorlog)

def warnUser(finalerrorlog):
    """Gives the user a sense of warnings and errors in the output file."""
    f = open(finalerrorlog,"r")
    Lines = f.readlines()
    f.close()

    error=False
    for i in range(0,len(Lines)):
        if "Error" in Lines[i]:
            print("error.log messages include:\n")
            print(Lines[i],end=" ")
            error=True
            break

    for i in range(0,len(Lines)):
        if ("Warning" in Lines[i]) or ("warning" in Lines[i]):
            if error == False:
                print("error.log messages include:\n")
                print(Lines[i])
                break
            else:
                print(Lines[i])
                break