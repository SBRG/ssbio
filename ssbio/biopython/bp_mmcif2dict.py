import shlex
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


class MMCIF2DictFix(MMCIF2Dict):
    """Fixes for MMCIF2Dict according to biopython#481 and biopython#523. Methods override parent.
    """

    def __init__(self, filename):
        super(MMCIF2DictFix, self).__init__(filename=filename)

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