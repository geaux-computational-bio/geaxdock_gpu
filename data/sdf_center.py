#!/usr/bin/env python

"""
parse a sdf file and calculate the geometric center of the ligand,
does not count the Hydrogen atom,
add the centers to the first line of ff file in the form of
"CENTER <xxx> <yyy> <zzz>"
"""
import shutil

# x, y, z are in the 0th, 1st, and 2nd collumns
options = {'x': 0,
           'y': 1,
           'z': 2}


def getNonHydrogenLines(sdf_ifn):
    lines = []
    for line in file(sdf_ifn):
        if 'M  END' not in line:
            lines.append(line)
        else:
            break

    lna = int(lines[3].split()[0])  # number of atoms
    atom_lines = lines[4: (4 + lna)]
    non_H_lines = [l for l in atom_lines if 'H' not in l]
    return non_H_lines


def getCoordVals(lines, coord='x'):
    idx = options[coord]
    vals = [line.split()[idx] for line in lines]
    vals = [float(val) for val in vals]
    return vals


def getCenters(non_H_lines):
    centers = []
    for coord in ['x', 'y', 'z']:
        vals = getCoordVals(non_H_lines, coord=coord)
        center = sum(vals) / float(len(vals))
        centers.append(center)

    return centers


def ligCenter(sdf_ifn):
    non_H_lines = getNonHydrogenLines(sdf_ifn)
    return getCenters(non_H_lines)


def addCenter(sdf_ifn, ff_ifn):
    # back up the ff file
    bk = ff_ifn + '.bk'
    shutil.copyfile(ff_ifn, bk)

    lines = file(ff_ifn).readlines()
    if "CENTER" in lines[0]:
        pass
    else:
        centers = ligCenter(sdf_ifn)
        centers = map(str, centers)
        center_line = "CENTER " + ' '.join(centers) + "\n"
        lines.insert(0, center_line)

        with open(ff_ifn, 'w') as f:
            f.writelines(lines)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="add ligand center to the \
    force field file")
    parser.add_argument("-l", '--sdf', type=str,
                        help="ligand sdf file")
    parser.add_argument("-f", '--ff', type=str,
                        help="force field file")
    args = parser.parse_args()

    sdf_ifn = args.sdf
    ff_ifn = args.ff

    addCenter(sdf_ifn, ff_ifn)
