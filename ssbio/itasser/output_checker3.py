import os
import glob
import argparse
import re
import subprocess
import ast


class OutputChecker():
    """
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
            self.working_dir = input("Working directory: ")
        else:
            self.working_dir = working_dir

        if not input_file_extension:
            self.input_file_extension = input("Input files: ")
        else:
            self.input_file_extension = input_file_extension

        if not output_file_extension:
            self.output_file_extension = input("Output files: ")
        else:
            self.output_file_extension = output_file_extension

        if not complete_phrase:
            self.complete_phrase = input("Complete phrase: ")
        else:
            self.complete_phrase = complete_phrase

        self.display_files = display_files
        self.fast = fast

    def file_searcher(self, ext):
        """
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
        """
        Simple code to grep for a pattern in text
        """
        grepper = re.compile(pattern)
        if grepper.search(str(text)):
            return True
        return False

    def grep_file(self, pattern, file_obj):
        """
        From: http://code.activestate.com/recipes/577069-access-grep-from-python
        """
        for line_num, line in enumerate(file_obj):
            if self.grep_text(line, pattern):
                return True
        return False

    def tail(self, f, n):
        process = subprocess.Popen(
            ['tail', '-n {0}'.format(n), '{0}'.format(f)], stdout=subprocess.PIPE)
        stdout = process.communicate()[0]
        return stdout

    def outfiles_done(self):
        """
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
            print(complete)

        return complete


if __name__ == '__main__':
    p = argparse.ArgumentParser(description='Output Parser v0.1')
    p.add_argument('--working_dir', '-wd')
    p.add_argument('--input', '-i')
    p.add_argument('--output', '-o')
    p.add_argument('--complete', '-c')
    p.add_argument('--display', '-d', action='store_true', default=False)
    p.add_argument('--fast', '-f', action='store_true', default=False)
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

    print('Number of input files: {0}'.format(len(infiles)))
    print('Number of output files: {0}'.format(len(outfiles)))
    completed_files = oc.outfiles_done()
    print('Number of completed output files: {0}'.format(len(completed_files)))

    # print(outfiles)
    # print(completed_files)
