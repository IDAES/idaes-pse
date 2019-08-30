#!/usr/bin/python


def main():
    import alamopy.examples
    infile = 'input.txt'
    outfile = 'output.txt'
    fin = open(infile, 'r')
    fout = open(outfile, 'w')
    newline = fin.readline()
    newlist = newline.split()
    n = int(newlist[0])
    for p in range(0, n):
        newline = fin.readline()
        newlist = newline.split()
        ninputs = len(newlist)
        x = [0] * (ninputs + 1)
        for k in range(0, ninputs):
            x[k] = float(newlist[k])
        x[ninputs] = alamopy.examples.sixcamel(*x[:-1])
        for k in range(0, len(x)):
            fout.write(str(x[k]) + ' ')
        fout.write(' \n')


if __name__ == '__main__':
    main()
