#!/usr/bin/env python3
import os
import sys
def check_input(args):
    start = int(args[0][0:])
    fh = open(args[1], 'r')
    return (start, fh)

def pad_line(line):
    """Helper function to pad line to 80 characters in case it is shorter"""
    size_of_line = len(line)
    if size_of_line < 80:
        padding = 80 - size_of_line + 1
        line = line.strip('\n') + ' ' * padding + '\n'
    return line[:81]  # 80 + newline character

def renumber_residues(fh, start):
    terline='TER'+' '*77
    serial = start
    prev_line_res = None  
    res = start - 1 
    _pad_line=pad_line
    for line in fh:
        line = _pad_line(line)
        if line.startswith('ATOM'):
            line_res = line[17:27]
            if line_res != prev_line_res:
                prev_line_res = line_res
                res += 1
            yield line[:22] + str(res).rjust(4) + line[26:]

        elif line.startswith('TER'):
            yield terline+'\n'

def main():
    starting_res, pdbfh = check_input(sys.argv[1:])
    new_pdb = renumber_residues(pdbfh, starting_res)

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
