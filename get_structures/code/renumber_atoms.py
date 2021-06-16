#!/usr/bin/env python3
import os
import sys
def check_input(args):
    start = int(args[0][0:])
    fh = open(args[1], 'r')
    return (start, fh)

def renumber_atoms(fh, start):
    terline='TER'+' '*77
    atom = start
    for line in fh:
        if line.startswith('ATOM'):
            yield line[:6] + str(atom).rjust(5) + line[11:]
            atom += 1

        elif line.startswith('TER'):
            yield terline+'\n'

def main():
    starting_atom, pdbfh = check_input(sys.argv[1:])
    new_pdb = renumber_atoms(pdbfh, starting_atom)

    try:
        _buffer = []
        _buffer_size = 5000  # write N lines at a time
        for lineno, line in enumerate(new_pdb):
            if not (lineno % _buffer_size):
                sys.stdout.write(''.join(_buffer))
                _buffer = []
            _buffer.append(line)

        sys.stdout.write(''.join(_buffer))
        sys.stdout.flush()
    except IOError:
        pass
    pdbfh.close()
    sys.exit(0)


if __name__ == '__main__':
    main()
