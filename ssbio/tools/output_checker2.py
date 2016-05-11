from __future__ import absolute_import
import os
import glob
import argparse
import re
import subprocess
import ast
from io import open

class OutputChecker(object):
    u"""
    Ever run a program or script that generates output files, and wonder if they
    completed successfully? output_checker.py will help you out by:
    1) Asking you what the extension of the input files are (ex: "*.in")
    1) Asking you what the extension of the output files are (ex: "*.out")
    2) Asking you what should be contained in those files which indicate
        completeness (ex: "SCRIPT DONE")
    3) Looping over these files and searching for this phrase
    4) Reporting progress
    """

    def __init__(self, working_dir=None, input_file_extension=None, output_file_extension=None, complete_phrase=None, display_files=None, fast=None):
        if not working_dir:
            self.working_dir = raw_input(u"Working directory: ")
        else:
            self.working_dir = working_dir

        if not input_file_extension:
            self.input_file_extension = raw_input(u"Input files: ")
        else:
            self.input_file_extension = input_file_extension

        if not output_file_extension:
            self.output_file_extension = raw_input(u"Output files: ")
        else:
            self.output_file_extension = output_file_extension

        if not complete_phrase:
            self.complete_phrase = raw_input(u"Complete phrase: ")
        else:
            self.complete_phrase = complete_phrase

        self.display_files = display_files
        self.fast = fast

    def file_searcher(self, ext):
        u"""
        Returns a list of all files which have provided extension in the working directory
        Output:
        files - list of files
        """
        wd = self.working_dir

        searchstring = os.path.join(wd, ext)
        # print("Searching for files like so: {}".format(searchstring))
        files = glob.glob(searchstring)
        return files

    def grep_text(self, text, pattern):
        u"""
        Simple code to grep for a pattern in text
        """
        grepper = re.compile(pattern)
        if grepper.search(unicode(text)):
            return True
        return False

    def grep_file(self, pattern, file_obj):
        u"""
        From: http://code.activestate.com/recipes/577069-access-grep-from-python
        """
        for line_num, line in enumerate(file_obj):
            if self.grep_text(line, pattern):
                return True
        return False

    def tail(self, f, n):
        process = subprocess.Popen(
            [u'tail', u'-n {0}'.format(n), u'{0}'.format(f)], stdout=subprocess.PIPE)
        stdout = process.communicate()[0]
        return stdout

    def outfiles_done(self):
        u"""
        Returns the output files that have the complete phrase in them
        """
        complete = []

        outfiles = self.file_searcher(self.output_file_extension)
        for f in outfiles:
            if self.fast:
                tail_of_file = self.tail(f, 10)
                if self.grep_text(tail_of_file, self.complete_phrase):
                    complete.append(f)
            else:
                outfile = open(f)
                if self.grep_file(self.complete_phrase, outfile):
                    complete.append(f)

        if self.display_files:
            print complete

        return complete


if __name__ == u'__main__':
    p = argparse.ArgumentParser(description=u'Output Parser v0.1')
    p.add_argument(u'--working_dir', u'-wd')
    p.add_argument(u'--input', u'-i')
    p.add_argument(u'--output', u'-o')
    p.add_argument(u'--complete', u'-c')
    p.add_argument(u'--display', u'-d', action=u'store_true', default=False)
    p.add_argument(u'--fast', u'-f', action=u'store_true', default=False)
    args = p.parse_args()

    oc = OutputChecker(args.working_dir, args.input,
                       args.output, args.complete, args.display, args.fast)
    # print('\nWorking directory: {}'.format(oc.working_dir))
    # print('Input files extension: {}'.format(oc.input_file_extension))
    # print('Output files extension: {}'.format(oc.output_file_extension))
    # print('Complete phrase to look for: {}'.format(oc.complete_phrase))
    # print('Display completed files: {}'.format(oc.display_files))
    # print('Fast searcher (search from tail): {}'.format(oc.fast))

    infiles = oc.file_searcher(oc.input_file_extension)
    outfiles = oc.file_searcher(oc.output_file_extension)

    print u'Number of input files: {0}'.format(len(infiles))
    print u'Number of output files: {0}'.format(len(outfiles))
    completed_files = oc.outfiles_done()
    print u'Number of completed output files: {0}'.format(len(completed_files))

    # print(outfiles)
    # print(completed_files)
