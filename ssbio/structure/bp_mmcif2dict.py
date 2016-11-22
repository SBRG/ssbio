import shlex
from Bio.PDB import MMCIF2Dict


class MMCIF2DictFix(MMCIF2Dict):

    def __init__(self, filename):
        super(MMCIF2DictFix, self).__init__()

    def _tokenize(self, handle):
        for line in handle:
            if line.startswith("#"):
                continue
            elif line.startswith(";"):
                token = line[1:].strip()
                for line in handle:
                    line = line.strip()
                    if line == ';':
                        break
                    token += line
                yield token
            else:
                try:
                    tokens = shlex.split(line)
                except ValueError:
                    # error "No closing quotation"
                    line = line.replace("'",'"')
                    # if odd - add a closing " to that line
                    if not line.count('"') % 2 == 0:
                        line = '{}"'.format(line)
                    tokens = shlex.split(line)
                for token in tokens:
                    yield token