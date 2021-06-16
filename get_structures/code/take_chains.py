#!/usr/bin/env python3
import os
import sys
def check_input(args):
    chains = args[0][1:]
    chainset = set([c.strip() for c in chains.split(',') if c.strip().isalnum()])
    fh = open(args[1], 'r')
    return (chains, fh)

def take_chains(fh, chainset):
    ter=True
    terline='TER'+' '*77
    for line in fh:
        if line.startswith('ATOM'):
            if line[21] in chainset:
                yield line
                ter=False
        elif line.startswith('TER'):
            if not ter:
                yield terline+'\n'
                ter=True

def main():
    chainset, pdbfh = check_input(sys.argv[1:])
    new_pdb = take_chains(pdbfh, chainset)

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
